"""Descriptores PaDEL 1D/2D reales y alineación por identificador molecular."""
import os
import shutil
import argparse
import hashlib

import numpy as np
import pandas as pd
from padelpy import padeldescriptor

from .data import ROOT, write_json


def configure_java():
    if shutil.which("java"):
        return
    candidates = list((ROOT / ".tools/java").glob("**/bin/java.exe"))
    if candidates:
        os.environ["PATH"] = str(candidates[0].parent) + os.pathsep + os.environ["PATH"]
    else:
        raise RuntimeError("PaDEL requiere Java 8+ en PATH (ver README.md).")


def calculate(batch_size=100):
    configure_java()
    source = ROOT / "data/processed/mmp9_binary.csv"
    df = pd.read_csv(source)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    interim = ROOT / "data/interim" / digest[:12]
    interim.mkdir(parents=True, exist_ok=True)
    parts = []
    for start in range(0, len(df), batch_size):
        batch = df.iloc[start:start + batch_size]
        inp = interim / f"batch_{start:05d}.smi"
        out = inp.with_suffix(".csv")
        if not out.exists():
            pending = inp.with_suffix(".partial.csv")
            batch[["smiles", "molecule_id"]].to_csv(inp, sep="\t", header=False, index=False)
            padeldescriptor(mol_dir=str(inp), d_file=str(pending), d_2d=True,
                            fingerprints=False, threads=4, maxruntime=30000,
                            retainorder=True, sp_timeout=1800)
            validated = pd.read_csv(pending)
            assert validated.Name.is_unique and set(validated.Name) == set(batch.molecule_id)
            pending.replace(out)
        part = pd.read_csv(out)
        assert part.Name.is_unique and set(part.Name) == set(batch.molecule_id)
        parts.append(part)
        print(f"PaDEL: {min(start + batch_size, len(df))}/{len(df)}", flush=True)
    features = pd.concat(parts).set_index("Name").reindex(df.molecule_id)
    features = features.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    # No eliminar columnas con información del test: se filtran dentro del pipeline.
    bad = features.isna().all(axis=1)
    if bad.any():
        raise RuntimeError(f"PaDEL falló completamente para: {features.index[bad].tolist()}")
    features.index.name = "molecule_id"
    features.to_csv(ROOT / "data/processed/padel_descriptors.csv.gz", compression="gzip")
    write_json(ROOT / "data/processed/features.json", {
        "representation": "PaDEL 1D/2D", "n_molecules": len(features),
        "n_descriptors": features.shape[1], "input_sha256": digest,
        "missing_fraction": float(features.isna().mean().mean()),
        "molecules_with_missing": int(features.isna().any(axis=1).sum()),
        "settings": {"d_2d": True, "d_3d": False, "fingerprints": False,
                     "maxruntime_ms": 30000, "threads": 4},
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=100)
    calculate(parser.parse_args().batch_size)
