"""100 particiones pareadas, selección interna y evaluación externa de modelos."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
from xgboost import XGBClassifier

from .data import ROOT, write_json

METRICS = ["accuracy", "f1", "sensitivity", "specificity"]
CONFIGS = {
    "LogisticRegression": [{"C": 0.1}, {"C": 1.0}],
    "RandomForest": [{"n_estimators": 100, "max_depth": 12, "min_samples_leaf": 2},
                     {"n_estimators": 150, "max_depth": None, "min_samples_leaf": 1}],
    "XGBoost": [{"n_estimators": 100, "max_depth": 3, "learning_rate": 0.1},
                {"n_estimators": 150, "max_depth": 5, "learning_rate": 0.05}],
}


def metrics(y_true, y_pred):
    """Positivo=activo (1). Matriz con orden [inactivo, activo]."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    if tn + fp == 0 or tp + fn == 0:
        raise ValueError("La evaluación requiere ambas clases.")
    return dict(accuracy=float(accuracy_score(y_true, y_pred)),
                f1=float(f1_score(y_true, y_pred, zero_division=0)),
                sensitivity=float(recall_score(y_true, y_pred, zero_division=0)),
                specificity=float(tn / (tn + fp)),
                precision=float(precision_score(y_true, y_pred, zero_division=0)),
                tn=int(tn), fp=int(fp), fn=int(fn), tp=int(tp))


def preprocessor(representation):
    if representation == "SMILES":
        # Respetar mayúsculas: C y c tienen significado químico diferente.
        return TfidfVectorizer(analyzer="char", ngram_range=(2, 4), lowercase=False,
                               min_df=2, max_features=2048, dtype=np.float32)
    return Pipeline([
        ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
        ("variance", VarianceThreshold()),
        ("scale", StandardScaler()),
    ])


def model(name, params, seed):
    if name == "LogisticRegression":
        return LogisticRegression(**params, solver="liblinear", max_iter=2000,
                                  class_weight="balanced", random_state=seed)
    if name == "RandomForest":
        return RandomForestClassifier(**params, class_weight="balanced", n_jobs=1,
                                      random_state=seed)
    return XGBClassifier(**params, objective="binary:logistic", eval_metric="logloss",
                         tree_method="hist", n_jobs=1, random_state=seed,
                         subsample=0.8, colsample_bytree=0.8, reg_lambda=1)


def fit_model(estimator, x, y, name):
    if name == "XGBoost":
        # Ponderación calculada solo en la partición usada para ajustar.
        estimator.set_params(scale_pos_weight=float((y == 0).sum() / (y == 1).sum()))
    estimator.fit(x, y)
    return estimator


def matrix(x):
    # Ceros explícitos para XGBoost: ausencia de n-grama no significa dato ausente.
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float32)


def split_indices(y, seed):
    train, test = train_test_split(np.arange(len(y)), test_size=0.2, stratify=y, random_state=seed)
    fit, valid = train_test_split(train, test_size=0.25, stratify=y[train], random_state=seed + 10000)
    return train, test, fit, valid


def run_seed(seed, molecules, padel, out):
    shard = out / "runs" / f"seed_{seed:03d}.json"
    if shard.exists():
        return json.loads(shard.read_text(encoding="utf-8"))
    y = molecules.label.to_numpy()
    train, test, fit, valid = split_indices(y, seed)
    split = {"seed": seed, **{k: molecules.molecule_id.iloc[v].tolist()
             for k, v in {"train": train, "test": test, "inner_fit": fit, "inner_validation": valid}.items()}}
    write_json(out / "splits" / f"seed_{seed:03d}.json", split)
    rows, trials, predictions = [], [], []
    for representation, raw in [("SMILES", molecules.smiles), ("PaDEL", padel)]:
        inner_prep = preprocessor(representation)
        x_fit = matrix(inner_prep.fit_transform(raw.iloc[fit]))
        x_valid = matrix(inner_prep.transform(raw.iloc[valid]))
        final_prep = preprocessor(representation)
        x_train = matrix(final_prep.fit_transform(raw.iloc[train]))
        x_test = matrix(final_prep.transform(raw.iloc[test]))
        for name, candidates in CONFIGS.items():
            scores = []
            for candidate, params in enumerate(candidates):
                estimator = fit_model(model(name, params, seed), x_fit, y[fit], name)
                score = f1_score(y[valid], estimator.predict(x_valid), zero_division=0)
                scores.append(score)
                trials.append(dict(seed=seed, representation=representation, model=name,
                                   candidate=candidate, params=params, validation_f1=float(score)))
            best = int(np.argmax(scores))
            estimator = fit_model(model(name, candidates[best], seed), x_train, y[train], name)
            probability = estimator.predict_proba(x_test)[:, 1]
            predicted = (probability >= 0.5).astype(int)
            row = dict(seed=seed, representation=representation, model=name,
                       candidate=best, validation_f1=float(scores[best]),
                       n_train=len(train), n_test=len(test), n_features=x_train.shape[1],
                       **metrics(y[test], predicted))
            rows.append(row)
            predictions.extend(dict(seed=seed, representation=representation, model=name,
                                    molecule_id=molecules.molecule_id.iloc[idx],
                                    y_true=int(y[idx]), y_pred=int(pred), probability=float(prob))
                               for idx, pred, prob in zip(test, predicted, probability))
            if seed == 0:
                (out / "models").mkdir(exist_ok=True)
                joblib.dump({"preprocessor": final_prep, "model": estimator,
                             "representation": representation, "positive_class": "active",
                             "descriptor_columns": padel.columns.tolist() if representation == "PaDEL" else None},
                            out / "models" / f"{representation}_{name}_seed0.joblib", compress=3)
                if name == "XGBoost":
                    estimator.save_model(out / "models" / f"{representation}_XGBoost_seed0.ubj")
    majority = int(np.bincount(y[train]).argmax())
    baseline = dict(seed=seed, **metrics(y[test], np.full(len(test), majority)))
    result = dict(rows=rows, trials=trials, predictions=predictions, baseline=baseline)
    write_json(shard, result)
    print(f"Completada semilla {seed}: 6 evaluaciones externas", flush=True)
    return result


def run(repeats=100, workers=3, output="results/main"):
    if repeats < 1:
        raise ValueError("repeats debe ser positivo")
    out = ROOT / output
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT / "data/processed/mmp9_binary.csv"
    fp = ROOT / "data/processed/padel_descriptors.csv.gz"
    molecules = pd.read_csv(source)
    padel = pd.read_csv(fp, index_col="molecule_id")
    assert padel.index.is_unique and set(padel.index) == set(molecules.molecule_id)
    padel = padel.reindex(molecules.molecule_id).reset_index(drop=True)
    assert molecules.smiles.is_unique and set(molecules.label) == {0, 1}
    signature = {
        "data_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "features_sha256": hashlib.sha256(fp.read_bytes()).hexdigest(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "configs": CONFIGS, "seeds": list(range(repeats)),
        "split": "stratified 80/20 outer, 75/25 inner on outer training",
        "selection_metric": "inner validation F1, active=1; ties choose first candidate",
        "threshold": 0.5,
        "versions": {p: importlib.metadata.version(p) for p in
                     ["numpy", "pandas", "scikit-learn", "xgboost", "rdkit", "padelpy", "joblib"]},
    }
    manifest = out / "manifest.json"
    if manifest.exists() and json.loads(manifest.read_text(encoding="utf-8")) != signature:
        raise RuntimeError("La configuración cambió. Usar otra carpeta --output para evitar mezclar resultados.")
    write_json(manifest, signature)
    results = []
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(run_seed, seed, molecules, padel, out) for seed in range(repeats)]
        for future in as_completed(futures):
            results.append(future.result())
    scores = pd.DataFrame([row for result in results for row in result["rows"]])
    scores = scores.sort_values(["representation", "model", "seed"])
    scores.to_csv(out / "metrics.csv", index=False)
    summary = scores.groupby(["representation", "model"])[METRICS].agg(["mean", "std"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary.reset_index().to_csv(out / "summary.csv", index=False)
    pd.DataFrame([row for r in results for row in r["predictions"]]).sort_values(
        ["seed", "representation", "model", "molecule_id"]).to_csv(out / "predictions.csv.gz", index=False, compression="gzip")
    pd.DataFrame([row for r in results for row in r["trials"]]).sort_values(
        ["seed", "representation", "model", "candidate"]).to_csv(out / "tuning.csv", index=False)
    pd.DataFrame([r["baseline"] for r in results]).sort_values("seed").to_csv(out / "baseline.csv", index=False)
    print(summary.to_string(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--output", default="results/main")
    args = parser.parse_args()
    run(args.repeats, args.workers, args.output)
