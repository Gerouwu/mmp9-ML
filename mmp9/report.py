"""Genera informe sencillo (Markdown y HTML) a partir de resultados reales."""
import argparse
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

from .data import ROOT
from .experiment import METRICS

LABELS = {"accuracy": "Accuracy", "f1": "F1-score", "sensitivity": "Sensibilidad", "specificity": "Especificidad"}


def markdown_table(frame):
    header = "| " + " | ".join(frame.columns) + " |\n"
    separator = "| " + " | ".join(["---"] * len(frame.columns)) + " |\n"
    return header + separator + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in frame.itertuples(index=False, name=None))


def generate(results="results/main"):
    folder = ROOT / results
    out = ROOT / "reports"
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(folder / "metrics.csv")
    summary = pd.read_csv(folder / "summary.csv")
    baseline = pd.read_csv(folder / "baseline.csv")
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    curated = json.loads((ROOT / "data/processed/curation.json").read_text(encoding="utf-8"))
    features = json.loads((ROOT / "data/processed/features.json").read_text(encoding="utf-8"))
    provenance = json.loads((ROOT / "data/raw/provenance.json").read_text(encoding="utf-8"))
    status = json.loads((ROOT / "data/raw/chembl_status.json").read_text(encoding="utf-8"))
    molecules = pd.read_csv(ROOT / "data/processed/mmp9_binary.csv")
    repeats = len(manifest["seeds"])
    assert scores.groupby(["representation", "model"]).size().eq(repeats).all()
    assert len(scores) == 6 * repeats
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    table = summary[["representation", "model"]].rename(columns={"representation": "Representación", "model": "Modelo"})
    for metric in METRICS:
        table[LABELS[metric] + " (%)"] = [f"{100*m:.2f} ± {100*s:.2f}" for m, s in
                                           zip(summary[metric + "_mean"], summary[metric + "_std"])]
    table.to_csv(out / "rendimiento_resumen.csv", index=False)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ticks = summary.representation + "\n" + summary.model
    colors = ["#087e8b" if model == "XGBoost" else "#607b96" for model in summary.model]
    for ax, metric in zip(axes.flat, METRICS):
        ax.bar(np.arange(len(summary)), 100 * summary[metric + "_mean"], color=colors,
               yerr=100 * summary[metric + "_std"], capsize=4)
        ax.set(ylim=(0, 105), ylabel="Porcentaje", title=LABELS[metric])
        ax.set_xticks(np.arange(len(summary)), ticks, rotation=22, ha="right", fontsize=8)
        ax.axhline(100 * baseline[metric].mean(), color="#b04c38", ls="--", label="Clase mayoritaria")
        ax.legend(fontsize=8)
    fig.suptitle(f"MMP-9: rendimiento externo, {repeats} repeticiones (media ± DE)")
    fig.savefig(figures / "rendimiento.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), constrained_layout=True)
    for ax, representation in zip(axes, ["SMILES", "PaDEL"]):
        subset = scores[(scores.model == "XGBoost") & (scores.representation == representation)]
        # Media de conteos por ejecución; no son moléculas únicas agregadas.
        cm = np.array([[subset.tn.mean(), subset.fp.mean()], [subset.fn.mean(), subset.tp.mean()]])
        ax.imshow(cm, cmap="Blues")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i,j]:.2f}", ha="center", va="center",
                        color="white" if cm[i,j] > cm.max()/2 else "black", fontsize=14)
        ax.set(xticks=[0, 1], yticks=[0, 1], xticklabels=["Inactiva", "Activa"],
               yticklabels=["Inactiva", "Activa"], xlabel="Predicción", ylabel="Clase real",
               title=f"XGBoost · {representation}")
    fig.suptitle("Matriz de confusión: conteo medio por partición de prueba")
    fig.savefig(figures / "xgboost_confusion.png", dpi=160)
    plt.close(fig)

    # EDA descriptivo: estos valores no son entradas adicionales de los modelos.
    eda = molecules[["molecule_id", "class", "pic50"]].copy()
    mols = [Chem.MolFromSmiles(s) for s in molecules.smiles]
    for key, fn in [("MW", Descriptors.MolWt), ("LogP", Descriptors.MolLogP),
                    ("HDonors", Lipinski.NumHDonors), ("HAcceptors", Lipinski.NumHAcceptors)]:
        eda[key] = [fn(mol) for mol in mols]
    eda.to_csv(out / "eda_lipinski.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    counts = molecules["class"].value_counts().reindex(["active", "inactive"])
    bars = axes[0].bar(["Activas", "Inactivas"], counts, color=["#087e8b", "#b04c38"])
    axes[0].bar_label(bars)
    axes[0].set(title="Conjunto de clasificación", ylabel="Moléculas")
    for category, color in [("active", "#087e8b"), ("inactive", "#b04c38")]:
        axes[1].hist(molecules.loc[molecules["class"] == category, "pic50"], bins=18,
                     alpha=0.7, color=color, label=category)
    axes[1].set(title="Distribución de pIC50", xlabel="pIC50 = 9 − log10(IC50 en nM)", ylabel="Moléculas")
    axes[1].legend()
    fig.savefig(figures / "datos.png", dpi=160)
    plt.close(fig)

    xgb = summary.loc[summary.model.eq("XGBoost")].sort_values("f1_mean", ascending=False)
    best = xgb.iloc[0]
    leader = summary.sort_values("f1_mean", ascending=False).iloc[0]
    exclusions = pd.DataFrame(list(curated["excluded_by_reason"].items()), columns=["Motivo (orden secuencial)", "Registros"])
    selection = pd.read_csv(folder / "tuning.csv")
    chosen = scores.groupby(["representation", "model", "candidate"]).size().reset_index(name="Selecciones")
    chosen.columns = ["Representación", "Modelo", "Candidato", "Selecciones"]
    sections = []

    def section(title, text):
        sections.append((title, text))

    section("Objetivo y resultado principal", f"Se clasificaron moléculas frente a MMP-9 humana (ChEMBL CHEMBL321, UniProt P14780) "
            f"como activas o inactivas. Se ejecutaron {len(scores)} evaluaciones externas: dos representaciones, tres modelos y "
            f"{repeats} repeticiones por combinación. Dentro de XGBoost, {best.representation} obtuvo el mayor F1 medio "
            f"({100*best.f1_mean:.2f}%). El mayor F1 observado entre las seis combinaciones correspondió a "
            f"{leader.model} con {leader.representation} ({100*leader.f1_mean:.2f}%). Esta comparación es descriptiva, "
            "no una demostración de superioridad estadística ni una validación experimental de nuevos inhibidores.")
    section("Datos y curación — data.py [1]", f"Fuente: {status['chembl_db_version']}, descarga UTC {provenance['downloaded_utc']}. "
            f"Se descargaron {curated['raw_activities']} registros IC50. Se conservaron mediciones exactas (=), en nM, positivas "
            "y finitas, sin alertas de validez ni duplicados potenciales, de ensayos de unión (B) con confianza de asignación 9. "
            "Los filtros se aplicaron secuencialmente; un registro descartado se cuenta solo en su primer motivo.\n\n" + markdown_table(exclusions) +
            f"\n\nQuedaron {curated['eligible_measurements']} mediciones elegibles. Se normalizó el fragmento padre y se agruparon "
            "estructuras canónicas, conservando estereoquímica; se tomó la mediana de IC50 por estructura. "
            f"Se excluyeron {curated['conflicting_molecules']} molécula(s) con mediciones tanto activas como inactivas. "
            f"El conjunto curado contiene {curated['curated_molecules']} moléculas: {curated['class_counts']['active']} activas "
            f"(IC50 ≤ 1000 nM), {curated['class_counts']['inactive']} inactivas (IC50 ≥ 10000 nM) y "
            f"{curated['class_counts']['intermediate']} intermedias. Las intermedias no se entrenaron ni evaluaron. "
            f"Se usaron {len(molecules)} moléculas en ambas representaciones. Los identificadores, mediciones y exclusiones "
            "se conservan en data/raw y data/processed.\n\n![Datos](figures/datos.png)")
    section("Representaciones — features.py [2] y experiment.py [3]",
            "SMILES: TF-IDF de n-gramas de caracteres de longitud 2 a 4, respetando mayúsculas, frecuencia mínima de dos "
            "moléculas y máximo 2048 características. Vocabulario y ponderaciones se ajustan solo con entrenamiento. "
            "Es una representación textual; los n-gramas no equivalen necesariamente a subestructuras químicas.\n\n"
            f"PaDEL: {features['n_descriptors']} descriptores 1D/2D calculados con PaDEL-Descriptor real, sin huellas ni "
            f"descriptores 3D. Fracción de entradas ausentes: {100*features['missing_fraction']:.4f}%. "
            "Se imputan medianas, se eliminan características constantes y se estandarizan valores únicamente con entrenamiento. "
            "Las columnas completamente ausentes en entrenamiento se rellenan con cero y se eliminan por varianza. "
            "La alineación se comprueba por identificador molecular. IC50, pIC50, clase e identificadores no se incluyen como entradas.")
    config_text = "\n\n".join(f"{name}: " + "; ".join(f"candidato {i}: `{json.dumps(p)}`" for i, p in enumerate(params))
                                for name, params in manifest["configs"].items())
    section("Diseño experimental e hiperparámetros — experiment.py [3]",
            f"Semillas 0–{repeats-1}. En cada repetición: partición externa estratificada 80/20, con "
            f"{int(scores.n_train.iloc[0])} moléculas de entrenamiento y {int(scores.n_test.iloc[0])} de prueba. "
            "Dentro del 80% se separa 75/25 para ajustar dos candidatos y seleccionar el de mayor F1 de validación "
            "(empates: primer candidato). Después se reajustan preprocesamiento y modelo elegido sobre todo el 80% "
            "y se evalúa una vez sobre el 20%. Se reutilizan las mismas particiones entre representaciones y modelos. "
            "La selección interna usa una partición de validación, no validación cruzada de múltiples pliegues.\n\n" + config_text +
            "\n\nRegresión logística: solver liblinear y pesos balanceados. Random Forest: pesos balanceados. "
            "XGBoost: hist, subsample=0.8, colsample_bytree=0.8, reg_lambda=1 y scale_pos_weight="
            "n_inactivas/n_activas calculado solo en los datos de ajuste. Umbral de decisión fijo: 0.5. "
            "Se emplean ceros explícitos en la matriz TF-IDF para XGBoost. Los detalles y elecciones se guardan en "
            "manifest.json, tuning.csv y metrics.csv.\n\n" + markdown_table(chosen))
    section("Métricas y rendimiento", "La clase positiva es activa (1); inactiva es 0. Accuracy=(TP+TN)/N; "
            "F1=2TP/(2TP+FP+FN); sensibilidad=TP/(TP+FN); especificidad=TN/(TN+FP). "
            "La media y desviación estándar muestral (ddof=1) se calculan sobre las métricas de cada repetición; "
            "no a partir de una matriz de confusión global. Los valores siguientes son porcentajes (DE en puntos porcentuales).\n\n" +
            markdown_table(table) + "\n\n![Rendimiento](figures/rendimiento.png)\n\n"
            "![Matrices XGBoost](figures/xgboost_confusion.png)\n\n"
            "Las matrices muestran conteos medios por repetición; por ello contienen decimales. "
            "Una molécula puede aparecer en prueba en varias repeticiones.")
    section("Interpretación y limitaciones",
            f"La referencia que siempre predice la clase mayoritaria obtuvo accuracy={100*baseline.accuracy.mean():.2f}%, "
            f"F1={100*baseline.f1.mean():.2f}%, sensibilidad={100*baseline.sensitivity.mean():.2f}% y "
            f"especificidad={100*baseline.specificity.mean():.2f}%. El desbalance exige interpretar las cuatro métricas juntas. "
            f"Para XGBoost con {best.representation}, la sensibilidad media fue {100*best.sensitivity_mean:.2f}% "
            f"y la especificidad {100*best.specificity_mean:.2f}%.\n\n"
            "Las particiones aleatorias pueden compartir familias químicas entre entrenamiento y prueba. No se ha realizado "
            "validación externa, temporal ni separación por esqueletos moleculares: el rendimiento puede ser optimista para "
            "familias nuevas. Las 100 repeticiones se solapan, no son 100 muestras independientes; la DE representa "
            "variabilidad entre particiones y no un intervalo de confianza. La heterogeneidad de los ensayos persiste "
            "incluso tras filtrar calidad. Excluir valores censurados e intermedios limita el dominio evaluado. "
            "La búsqueda de hiperparámetros es pequeña y no garantiza el óptimo. Los modelos guardados corresponden "
            "a la semilla 0, sin seleccionarla por rendimiento; son artefactos reproducibles del experimento, "
            "no modelos clínicamente validados. El informe IEEE se preparará posteriormente.")
    section("Archivos Python citados y evidencia",
            "[1] [mmp9/data.py](../mmp9/data.py): descarga, verificación de objetivo, filtros, agregación y etiquetas.\n\n"
            "[2] [mmp9/features.py](../mmp9/features.py): cálculo PaDEL por lotes y verificación de alineación.\n\n"
            "[3] [mmp9/experiment.py](../mmp9/experiment.py): particiones, candidatos, entrenamiento, métricas y modelos.\n\n"
            "[4] [mmp9/report.py](../mmp9/report.py): este informe y sus figuras; EDA de Lipinski.\n\n"
            "[5] [tests/test_pipeline.py](../tests/test_pipeline.py): pruebas de umbrales, métricas, particiones y preprocesamiento.\n\n"
            "[6] [mmp9/verify.py](../mmp9/verify.py): auditoría de métricas contra predicciones, separación e integridad.\n\n"
            f"Resultados: [métricas individuales](../{results}/metrics.csv), [resumen](../{results}/summary.csv), "
            f"[hiperparámetros](../{results}/tuning.csv), [configuración y versiones](../{results}/manifest.json).\n\n"
            "El archivo docente original se conserva sin alteraciones en la raíz y no se ejecuta en el flujo nuevo.")
    section("Fuentes de datos y documentación",
            "[ChEMBL: objetivo MMP-9](https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL321.json), "
            "[servicios web ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services). "
            "Datos atribuidos a EMBL-EBI ChEMBL, CC BY-SA 3.0.\n\n"
            "[XGBoost: API Python](https://xgboost.readthedocs.io/en/stable/python/python_api.html), "
            "[PaDELPy](https://github.com/ecrl/padelpy), "
            "[RDKit](https://www.rdkit.org/docs/Cookbook.html), "
            "[TF-IDF](https://scikit-learn.org/stable/modules/generated/sklearn.feature_extraction.text.TfidfVectorizer.html), "
            "[matriz de confusión](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html).")
    md = "# Informe de rendimiento: clasificación de moléculas frente a MMP-9\n\n"
    md += "\n\n".join(f"## {title}\n\n{text}" for title, text in sections) + "\n"
    (out / "informe_rendimiento.md").write_text(md, encoding="utf-8")
    import markdown
    body = markdown.markdown(md, extensions=["tables", "fenced_code"])
    document = '''<!doctype html><html lang="es"><meta charset="utf-8"><title>MMP-9 · Rendimiento</title>
<style>body{font:16px/1.65 system-ui,sans-serif;color:#203245;max-width:1100px;margin:48px auto;padding:0 30px}
h1{font-size:30px;color:#075963}h2{font-size:21px;border-bottom:1px solid #ccd9dd;padding-bottom:8px;margin-top:35px}
table{border-collapse:collapse;width:100%;font-size:13px}td,th{border:1px solid #ccd9dd;padding:9px;text-align:left}
th{background:#eaf4f4}tr:nth-child(even){background:#f6f8fa}img{max-width:100%}a{color:#076a90}
code{font-size:12px;overflow-wrap:anywhere}@media print{body{margin:0;font-size:11px}h2{break-after:avoid}img,table{break-inside:avoid}}</style><body>'''
    (out / "informe_rendimiento.html").write_text(document + body + "</body></html>", encoding="utf-8")
    print(table.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results/main")
    generate(parser.parse_args().results)
