"""Runner of the biological validation panel (point 8 of the post-demo-prep robustness plan).

It runs the real pipeline (PDB path: Phase 1.5 -> Phase 2 engines -> Phase 3 consensus/annotation)
over each protein of ``src/validation/biological_panel.py`` and compares the sites that
``apply_workflow_filter`` accepts (``pasa_umbral=True``) against the real ground truth
documented in the literature, reporting recall separately for tier A and tier B (see
the docstring of ``biological_panel.py`` for the reason behind this separation).

This validation measures the QUALITY OF THE ANNOTATION/CONSENSUS (Phases 2+3). Phase A/3c
(real structural modelling via PyRosetta) was removed from the project's scope
2026-08-10 (project supervision feedback), so there is no longer a later phase
to exclude here. DeepPTMPred (ESM-2 + PyRosetta features) already dominates the
run time on its own (~10+ min even for the smallest protein of the
panel, 103 residues, confirmed real 2026-08-04; it is not proportional to the
protein's size).

Usage:
    conda run -n cnb_pipeline python scripts/validate_biological_panel.py [--only NOMBRE ...]

It requires DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN pointing at real conda envs with the
engines installed (see README.md, "Installation"). It mocks nothing, it is a real run.
"""

import argparse
import sys
from pathlib import Path
from typing import List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import run_fase1_5_structure, run_fase2_pdb_motors, run_fase3_pdb_annotation
from src.config.settings import Settings
from src.validation.biological_panel import PANEL, PanelEntry


def _accepted_positions(entry: PanelEntry, output_dir: Path) -> set:
    """Corre Fase 1.5 -> 2 -> 3 (real, sin mocks) y devuelve el set (posicion, tipo_ptm) de
    todo lo que ``apply_workflow_filter`` acepto (pasa_umbral=True)."""
    record = run_fase1_5_structure(entry.pdb_path, output_dir)
    deepmvp_results, deepptmpred_results = run_fase2_pdb_motors(record, output_dir)
    filtered, _, _ = run_fase3_pdb_annotation(record, deepmvp_results, deepptmpred_results, output_dir)
    return {(int(row["posicion"]), row["tipo_ptm"]) for _, row in filtered.iterrows()}


def _recall_by_tier(entry: PanelEntry, accepted: set) -> dict:
    result = {}
    for tier in ("A", "B"):
        sites = [s for s in entry.positives if s.tier == tier]
        if not sites:
            continue
        hits = sum(1 for s in sites if (s.position, s.ptm_type) in accepted)
        result[tier] = (hits, len(sites))
    return result


def _negative_control_report(entry: PanelEntry, accepted: set) -> List[str]:
    lines = []
    for site in entry.negatives:
        flagged = (site.position, site.ptm_type) in accepted
        estado = "FALSO POSITIVO (el pipeline SI lo acepto)" if flagged else "correctamente NO aceptado"
        lines.append(f"    control negativo {site.position} ({site.ptm_type}): {estado}")
    return lines


def main(argv: List[str] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only", nargs="*", default=None,
        help="Nombres de entradas del panel a correr (por defecto, todas). Ver "
             "src/validation/biological_panel.py::PANEL para los nombres validos.",
    )
    parser.add_argument("--output-dir", default="outputs/validation_panel")
    args = parser.parse_args(argv)

    entries = [e for e in PANEL if args.only is None or e.name in args.only]
    if not entries:
        print(f"ERROR: ningun nombre de --only coincide con el panel (validos: "
              f"{[e.name for e in PANEL]})", file=sys.stderr)
        return 1

    Settings.ensure_dirs()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    total_a_hits = total_a_n = total_b_hits = total_b_n = 0
    for entry in entries:
        print(f"\n=== {entry.name} ({entry.uniprot_accession}, {entry.length} aa) ===")
        accepted = _accepted_positions(entry, output_dir / entry.name)
        recall = _recall_by_tier(entry, accepted)
        for tier, (hits, n) in recall.items():
            print(f"  Tier {tier}: {hits}/{n} sitios reales recuperados ({100 * hits / n:.0f}%)")
            if tier == "A":
                total_a_hits += hits
                total_a_n += n
            else:
                total_b_hits += hits
                total_b_n += n
        for line in _negative_control_report(entry, accepted):
            print(line)

    print("\n=== Resumen global ===")
    if total_a_n:
        print(f"Tier A: {total_a_hits}/{total_a_n} ({100 * total_a_hits / total_a_n:.0f}%)")
    if total_b_n:
        print(f"Tier B: {total_b_hits}/{total_b_n} ({100 * total_b_hits / total_b_n:.0f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
