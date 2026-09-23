"""Descarga, trazabilidad y curación de IC50 de MMP-9 humana en ChEMBL."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
API = "https://www.ebi.ac.uk/chembl/api/data/"
TARGET = "CHEMBL321"


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def session():
    client = requests.Session()
    client.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])))
    return client


def download():
    folder = ROOT / "data/raw"
    folder.mkdir(parents=True, exist_ok=True)
    client = session()

    def get(url, params=None):
        response = client.get(url, params=params, timeout=120)
        response.raise_for_status()
        return response.json()

    target = get(API + f"target/{TARGET}.json")
    assert target["organism"] == "Homo sapiens"
    assert target["pref_name"] == "Matrix metalloproteinase-9"
    assert target["target_type"] == "SINGLE PROTEIN"
    assert any(c["accession"] == "P14780" for c in target["target_components"])
    write_json(folder / "target.json", target)
    write_json(folder / "chembl_status.json", get(API + "status.json"))
    queries = {
        "activity": {"target_chembl_id": TARGET, "standard_type": "IC50", "limit": 1000},
        "assay": {"target_chembl_id": TARGET, "limit": 1000},
    }
    for resource, params in queries.items():
        rows = []
        url = API + resource + ".json"
        while url:
            payload = get(url, params)
            rows.extend(payload["activities" if resource == "activity" else "assays"])
            print(resource, len(rows), "/", payload["page_meta"]["total_count"], flush=True)
            nxt = payload["page_meta"]["next"]
            url = urljoin(API, nxt) if nxt else None
            params = None
        write_json(folder / f"{resource}.json", rows)
    write_json(folder / "provenance.json", {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "api": API, "target": TARGET, "queries": queries,
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in folder.glob("*.json") if p.name != "provenance.json"},
        "source_license": "ChEMBL data: CC BY-SA 3.0; retain attribution to EMBL-EBI ChEMBL.",
    })


def canonical_parent(smiles):
    """Normaliza y elige el fragmento padre; conserva estereoquímica."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        if not any(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except (ValueError, RuntimeError):
        return None


def activity_class(value):
    return "active" if value <= 1000 else "inactive" if value >= 10000 else "intermediate"


def curate():
    raw = ROOT / "data/raw"
    out = ROOT / "data/processed"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(json.loads((raw / "activity.json").read_text(encoding="utf-8")))
    assays = pd.DataFrame(json.loads((raw / "assay.json").read_text(encoding="utf-8")))
    stats = {"raw_activities": len(df)}
    reasons = pd.Series("", index=df.index)
    numeric = pd.to_numeric(df.standard_value, errors="coerce")
    good_assays = assays.loc[(assays.confidence_score == 9) & (assays.assay_type == "B"), "assay_chembl_id"]
    rules = [
        (df.standard_relation.eq("="), "non_exact_relation"),
        (df.standard_units.eq("nM"), "units_not_nM"),
        (numeric.gt(0) & np.isfinite(numeric), "invalid_IC50"),
        (df.data_validity_comment.fillna("").eq(""), "validity_flag"),
        (pd.to_numeric(df.potential_duplicate, errors="coerce").fillna(0).eq(0), "potential_duplicate"),
        (df.assay_chembl_id.isin(good_assays), "assay_not_direct_binding_confidence9"),
    ]
    for mask, reason in rules:
        reasons.loc[reasons.eq("") & ~mask] = reason
    df["exclusion_reason"] = reasons
    df.loc[reasons.ne("")].to_csv(out / "excluded_activities.csv", index=False)
    stats["excluded_by_reason"] = reasons[reasons.ne("")].value_counts().to_dict()
    clean = df.loc[reasons.eq("")].copy()
    clean["ic50_nM"] = numeric.loc[clean.index]
    clean["smiles"] = clean.canonical_smiles.map(canonical_parent)
    stats["invalid_smiles"] = int(clean.smiles.isna().sum())
    clean = clean.dropna(subset=["smiles"])
    clean["measurement_class"] = clean.ic50_nM.map(activity_class)
    clean.to_csv(out / "eligible_measurements.csv", index=False)
    molecules, conflicts = [], []
    for smiles, group in clean.groupby("smiles", sort=True):
        classes = set(group.measurement_class)
        # Excluir conflictos que cruzan ambos umbrales extremos.
        conflict = {"active", "inactive"}.issubset(classes)
        value = float(group.ic50_nM.median())
        row = {
            "molecule_id": min(group.molecule_chembl_id),
            "all_chembl_ids": "|".join(sorted(set(group.molecule_chembl_id))),
            "smiles": smiles, "ic50_nM": value, "pic50": 9 - np.log10(value),
            "class": activity_class(value), "n_measurements": len(group),
            "ic50_min": float(group.ic50_nM.min()), "ic50_max": float(group.ic50_nM.max()),
            "activity_ids": "|".join(group.activity_id.astype(str)),
            "assay_ids": "|".join(sorted(set(group.assay_chembl_id))),
        }
        (conflicts if conflict else molecules).append(row)
    pd.DataFrame(conflicts).to_csv(out / "conflicting_molecules.csv", index=False)
    all_molecules = pd.DataFrame(molecules).sort_values("molecule_id").reset_index(drop=True)
    all_molecules.to_csv(out / "molecules_all_classes.csv", index=False)
    binary = all_molecules.loc[all_molecules["class"].ne("intermediate")].copy()
    binary["label"] = binary["class"].eq("active").astype(int)
    assert len(binary) >= 100 and binary.label.value_counts().min() >= 10
    assert binary.smiles.is_unique and binary.molecule_id.is_unique
    binary.to_csv(out / "mmp9_binary.csv", index=False)
    stats.update(eligible_measurements=len(clean), conflicting_molecules=len(conflicts),
                 curated_molecules=len(all_molecules), binary_molecules=len(binary),
                 class_counts=all_molecules["class"].value_counts().to_dict())
    write_json(out / "curation.json", stats)
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        download()
    curate()
