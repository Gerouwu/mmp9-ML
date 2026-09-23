"""Ejecuta el proyecto completo desde la raíz con una sola instrucción."""
import argparse
import hashlib
import json

from mmp9.data import ROOT, download, curate
from mmp9.features import calculate
from mmp9.experiment import run
from mmp9.verify import verify
from mmp9.report import generate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Actualizar la copia de ChEMBL")
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--output", default="results/main")
    args = parser.parse_args()
    if args.download or not (ROOT / "data/raw/activity.json").exists():
        download()
    curate()
    source = ROOT / "data/processed/mmp9_binary.csv"
    info = ROOT / "data/processed/features.json"
    descriptors = ROOT / "data/processed/padel_descriptors.csv.gz"
    cached = json.loads(info.read_text(encoding="utf-8")) if info.exists() else {}
    if not descriptors.exists() or cached.get("input_sha256") != hashlib.sha256(source.read_bytes()).hexdigest():
        calculate()
    run(args.repeats, args.workers, args.output)
    if args.repeats > 1:
        verify(args.output)
    generate(args.output)


if __name__ == "__main__":
    main()
