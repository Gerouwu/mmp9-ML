"""Audita resultados guardados sin volver a entrenar los modelos."""
import argparse
import hashlib
import json

import joblib
import numpy as np
import pandas as pd

from .data import ROOT, write_json
from .experiment import METRICS, metrics, split_indices, matrix


def verify(output="results/main"):
    folder = ROOT / output
    df = pd.read_csv(ROOT / "data/processed/mmp9_binary.csv")
    scores = pd.read_csv(folder / "metrics.csv")
    predictions = pd.read_csv(folder / "predictions.csv.gz")
    summary = pd.read_csv(folder / "summary.csv")
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for key, path in [("data_sha256", ROOT / "data/processed/mmp9_binary.csv"),
                      ("features_sha256", ROOT / "data/processed/padel_descriptors.csv.gz"),
                      ("code_sha256", ROOT / "mmp9/experiment.py")]:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == manifest[key]
    trials = pd.read_csv(folder / "tuning.csv")
    baseline = pd.read_csv(folder / "baseline.csv")
    seeds = manifest["seeds"]
    assert len(scores) == 6 * len(seeds)
    assert not scores.duplicated(["seed", "representation", "model"]).any()
    y = df.set_index("molecule_id").label
    for seed in seeds:
        split = json.loads((folder / "splits" / f"seed_{seed:03d}.json").read_text(encoding="utf-8"))
        assert not set(split["train"]) & set(split["test"])
        assert not set(split["inner_fit"]) & set(split["inner_validation"])
        assert set(split["inner_fit"]) | set(split["inner_validation"]) == set(split["train"])
        assert set(split["train"]) | set(split["test"]) == set(df.molecule_id)
        indices = split_indices(df.label.to_numpy(), seed)
        for key, idx in zip(["train", "test", "inner_fit", "inner_validation"], indices):
            assert df.molecule_id.iloc[idx].tolist() == split[key]
        base = baseline.loc[baseline.seed.eq(seed)].iloc[0]
        majority = int(y.loc[split["train"]].value_counts().idxmax())
        recalculated = metrics(y.loc[split["test"]], [majority] * len(split["test"]))
        for key, value in recalculated.items():
            assert np.isclose(base[key], value)
        for (representation, name), group in predictions.loc[predictions.seed.eq(seed)].groupby(["representation", "model"]):
            assert group.molecule_id.is_unique
            assert set(group.molecule_id) == set(split["test"])
            assert np.array_equal(group.y_true, y.loc[group.molecule_id])
            assert group.probability.between(0, 1).all()
            assert np.array_equal(group.y_pred, (group.probability >= 0.5).astype(int))
            row = scores[(scores.seed == seed) & (scores.representation == representation) & (scores.model == name)].iloc[0]
            recalculated = metrics(group.y_true, group.y_pred)
            for key, value in recalculated.items():
                assert np.isclose(row[key], value), (seed, representation, name, key)
            candidates = trials[(trials.seed == seed) & (trials.representation == representation) & (trials.model == name)].sort_values("candidate")
            assert len(candidates) == 2
            assert int(row.candidate) == int(candidates.loc[candidates.validation_f1.idxmax(), "candidate"])
    for _, row in summary.iterrows():
        subset = scores[(scores.representation == row.representation) & (scores.model == row.model)]
        assert set(subset.seed) == set(seeds)
        for key in METRICS:
            assert np.isclose(row[key + "_mean"], subset[key].mean())
            assert np.isclose(row[key + "_std"], subset[key].std(ddof=1))
    if 0 in seeds:
        descriptors = pd.read_csv(ROOT / "data/processed/padel_descriptors.csv.gz", index_col="molecule_id")
        indexed = df.set_index("molecule_id")
        for (representation, name), group in predictions.loc[predictions.seed.eq(0)].groupby(["representation", "model"]):
            saved = joblib.load(folder / "models" / f"{representation}_{name}_seed0.joblib")
            raw = indexed.loc[group.molecule_id, "smiles"] if representation == "SMILES" else descriptors.loc[
                group.molecule_id, saved["descriptor_columns"]]
            transformed = matrix(saved["preprocessor"].transform(raw))
            probability = saved["model"].predict_proba(transformed)[:, 1]
            assert np.allclose(group.probability, probability, rtol=1e-6, atol=1e-7)
    result = {"status": "passed", "evaluations": len(scores), "prediction_rows": len(predictions),
              "seeds": len(seeds), "checks": ["input_hashes", "split_disjointness", "split_reproduction", "labels", "threshold",
              "confusion_counts", "all_metrics", "sample_standard_deviation", "hyperparameter_selection", "baseline", "model_reloading"]}
    write_json(folder / "verification.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/main")
    verify(parser.parse_args().output)
