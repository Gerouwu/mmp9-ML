"""Exporta matrices y gráficas PNG de las seis combinaciones evaluadas."""
import json
from zipfile import ZipFile, ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import ROOT
from .experiment import METRICS

NAMES = {"LogisticRegression": "Regresión logística", "RandomForest": "Random Forest", "XGBoost": "XGBoost"}
LABELS = ["Accuracy", "F1-score", "Sensibilidad", "Especificidad"]
COLORS = ["#426b91", "#087e8b", "#d18b32", "#995879"]


def draw_matrix(ax, group, normalized=False):
    matrices = group[["tn", "fp", "fn", "tp"]].to_numpy().reshape(-1, 2, 2)
    values = (100 * matrices / matrices.sum(axis=2, keepdims=True)).mean(axis=0) if normalized else matrices.mean(axis=0)
    vmax = 100 if normalized else 144
    ax.imshow(values, cmap="Blues", vmin=0, vmax=vmax)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{values[i,j]:.2f}" + ("%" if normalized else ""),
                    ha="center", va="center", fontsize=15,
                    color="white" if values[i,j] > vmax * .55 else "#17324a")
    first = group.iloc[0]
    ax.set(xticks=[0, 1], yticks=[0, 1], xticklabels=["Inactiva", "Activa"],
           yticklabels=["Inactiva", "Activa"], xlabel="Clase predicha", ylabel="Clase real",
           title=f"{NAMES[first.model]} · {first.representation}")


def generate():
    scores = pd.read_csv(ROOT / "results/main/metrics.csv")
    assert len(scores) == 600
    groups = [(rep, model, scores[(scores.representation == rep) & (scores.model == model)])
              for rep in ["SMILES", "PaDEL"] for model in NAMES]
    assert all(len(group) == 100 for _, _, group in groups)
    out = ROOT / "reports/imagenes"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    files = []

    def save(fig, name):
        fig.savefig(out / name, dpi=220, facecolor="white")
        plt.close(fig)
        files.append(name)

    for normalized in [False, True]:
        suffix = "porcentaje" if normalized else "conteos"
        description = "Porcentaje medio por clase real" if normalized else "Conteo medio por partición de prueba"
        fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
        for ax, (rep, model, group) in zip(axes.flat, groups):
            draw_matrix(ax, group, normalized)
            single, single_ax = plt.subplots(figsize=(6, 5.6), constrained_layout=True)
            draw_matrix(single_ax, group, normalized)
            single.suptitle(f"MMP-9 · {description}\n100 repeticiones · positiva = activa", fontsize=11)
            save(single, f"matriz_{rep}_{model}_{suffix}.png")
        fig.suptitle(f"MMP-9 · Matrices de confusión\n{description} · 100 repeticiones · positiva = activa", fontsize=16)
        save(fig, f"matrices_todas_{suffix}.png")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for ax, (rep, model, group) in zip(axes.flat, groups):
        means = group[METRICS].mean().to_numpy() * 100
        stds = group[METRICS].std(ddof=1).to_numpy() * 100

        def draw_bars(target):
            target.bar(range(4), means, yerr=stds, capsize=5, color=COLORS)
            for i, (mean, std) in enumerate(zip(means, stds)):
                target.text(i, mean + std + 2, f"{mean:.2f}%", ha="center", fontsize=10)
            target.set(xticks=range(4), xticklabels=LABELS, ylim=(0, 112),
                       ylabel="Porcentaje", title=f"{NAMES[model]} · {rep}")
            target.tick_params(axis="x", labelrotation=18)
            target.grid(axis="y", alpha=.15)
            target.set_axisbelow(True)

        draw_bars(ax)
        single, single_ax = plt.subplots(figsize=(8, 5.6), constrained_layout=True)
        draw_bars(single_ax)
        single.suptitle("MMP-9 · Media ± desviación estándar · 100 repeticiones", fontsize=12)
        save(single, f"metricas_{rep}_{model}.png")
    fig.suptitle("MMP-9 · Rendimiento de las seis combinaciones\nMedia ± desviación estándar · 100 repeticiones", fontsize=16)
    save(fig, "metricas_todas.png")

    for metric, label in zip(METRICS, LABELS):
        fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
        means = [group[metric].mean() * 100 for _, _, group in groups]
        stds = [group[metric].std(ddof=1) * 100 for _, _, group in groups]
        ax.barh(range(6), means, xerr=stds, capsize=4,
                color=["#426b91"] * 3 + ["#087e8b"] * 3)
        ax.set(yticks=range(6), yticklabels=[f"{rep} · {NAMES[model]}" for rep, model, _ in groups],
               xlim=(0, 115), xlabel=f"{label} (%)", title=f"MMP-9 · Comparación de {label}\nMedia ± DE · 100 repeticiones")
        for i, (mean, std) in enumerate(zip(means, stds)):
            ax.text(mean + std + 1, i, f"{mean:.2f} ± {std:.2f}", va="center", fontsize=10)
        ax.invert_yaxis()
        save(fig, f"comparacion_{metric}.png")

    note = ("# Imágenes del experimento MMP-9\n\n"
            "Fuente: results/main/metrics.csv. Generador: mmp9/plots.py.\n\n"
            "Las matrices de conteos son medias de 100 particiones (183 moléculas por prueba), "
            "por eso tienen decimales. No representan 18.300 moléculas distintas. "
            "Las matrices porcentuales se normalizan por clase real dentro de cada repetición y luego se promedian. "
            "Filas: clase real; columnas: predicción; clase positiva: activa.\n\n"
            "Las barras muestran media ± desviación estándar muestral; no son intervalos de confianza. "
            "Las repeticiones comparten moléculas y no son independientes.\n\n" +
            "\n".join(f"- [{name}]({name})" for name in files) + "\n")
    (out / "README.md").write_text(note, encoding="utf-8")
    with ZipFile(ROOT / "reports/imagenes_mmp9.zip", "w", ZIP_DEFLATED) as archive:
        for name in files + ["README.md"]:
            archive.write(out / name, arcname=name)
    print(json.dumps({"png_count": len(files), "directory": str(out)}, indent=2))


if __name__ == "__main__":
    generate()
