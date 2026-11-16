"""Phase A of EMNGly's go/no-go check 2: it downloads an
AlphaFold model for each unique UniProt accession of the N-GlyDE independent set
(``dukkakc/DeepNGlyPred``, ``EMNgly/data/N-GlyDE/NGLYDE_independent.txt``) and runs
real Phase 1.5 (``src.utils.structure_parser.parse_structure``) over each one, in order to
produce the derived sequence + position mapping table + single-chain
PDB that ``scripts/verify_emngly_nglyde_mcc.py`` subsequently consumes.

Why AlphaFold and not real crystallographic PDBs: EMNGly's own training
dataset uses AlphaFold2 models (see the docstring of
``src/engines/_emngly_runner.py``), so reproducing the published benchmark
(MCC~=0.736 on N-GlyDE) requires the same kind of structure that the paper uses, not real
crystallographic structures (which in addition do not exist for most of the
86 proteins of this set). Go/no-go check 1 (alignment against real GlyGen
sites) already covers the case of PDBs with gaps/non-continuous numbering separately.

It runs in the pipeline's main venv (``cnb_pipeline``, which has gemmi/pandas),
NEVER in ``.venv-emngly`` (fair-esm/torch/sklearn, without gemmi).
"""

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils.structure_parser import parse_structure
from src.utils.exceptions import StructureParsingError

REPO_ROOT = Path(__file__).resolve().parent.parent
DATASET_TSV = REPO_ROOT / "EMNgly" / "data" / "N-GlyDE" / "NGLYDE_independent.txt"
STRUCTURES_DIR = REPO_ROOT / "EMNgly" / "data" / "N-GlyDE" / "structures"
# NUNCA hardcodear 'AF-{acc}-F1-model_v4.pdb' -- confirmado real 2026-08-07
# (todas las 86 proteinas de este dataset dieron 404): AlphaFold DB ya sirve
# v6, no v4. Mismo hallazgo y mismo fix ya aplicado en
# scripts/generate_deepptmpred_calibration.py (ver su docstring) -- resolver
# la URL real via el API en vez de asumir un numero de version.
ALPHAFOLD_API_URL = "https://www.alphafold.ebi.ac.uk/api/prediction/{acc}"
REQUEST_TIMEOUT_SECONDS = 30


def _download_alphafold_model(accession: str, dest: Path) -> bool:
    if dest.is_file():
        return True
    try:
        req = urllib.request.Request(ALPHAFOLD_API_URL.format(acc=accession))
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
            entries = json.loads(resp.read())
        if not entries:
            print(f"[prepare] {accession}: sin modelo AlphaFold (API vacia), se omite.")
            return False
        entry = next((e for e in entries if e.get("uniprotAccession") == accession), entries[0])
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            print(f"[prepare] {accession}: API sin 'pdbUrl', se omite.")
            return False

        with urllib.request.urlopen(pdb_url, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
            dest.write_bytes(resp.read())
        return True
    except urllib.error.HTTPError as exc:
        print(f"[prepare] {accession}: sin modelo AlphaFold ({exc.code}), se omite.")
        return False
    except urllib.error.URLError as exc:
        print(f"[prepare] {accession}: error de red ({exc}), se omite.")
        return False


def main() -> int:
    df = pd.read_csv(DATASET_TSV, sep="\t")
    accessions = sorted(df["UniProt_ID"].unique())
    print(f"[prepare] {len(accessions)} accession(es) unico(s) en el dataset.")

    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    ok, skipped = 0, []
    for i, accession in enumerate(accessions, start=1):
        out_dir = STRUCTURES_DIR / accession
        derived_fasta = out_dir / f"{accession}_derived.fasta"
        if derived_fasta.is_file():
            ok += 1
            continue

        pdb_path = out_dir / f"{accession}_af.pdb"
        out_dir.mkdir(parents=True, exist_ok=True)
        if not _download_alphafold_model(accession, pdb_path):
            skipped.append(accession)
            continue

        try:
            parse_structure(pdb_path, out_dir)
            ok += 1
        except (StructureParsingError, FileNotFoundError) as exc:
            print(f"[prepare] {accession}: fallo Fase 1.5 ({exc}), se omite.")
            skipped.append(accession)

        if i % 10 == 0:
            print(f"[prepare] {i}/{len(accessions)} procesado(s)...")
        time.sleep(0.2)  # cortesia con el servidor de AlphaFold DB

    print(f"[prepare] listo: {ok} OK, {len(skipped)} omitido(s): {skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
