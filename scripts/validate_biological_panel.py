"""Runner del panel de validacion biologica (punto 8 del plan de robustez post-demo-prep).

Corre el pipeline real (Camino PDB: Fase 1.5 -> Fase 2 motores -> Fase 3 consenso/anotacion)
sobre cada proteina de ``src/validation/biological_panel.py`` y compara los sitios que
``apply_workflow_filter`` acepta (``pasa_umbral=True``) contra el ground truth real
documentado en literatura, reportando recall por separado para tier A y tier B (ver
docstring de ``biological_panel.py`` para el porque de esta separacion).

Deliberadamente NO corre Fase A (modelado estructural, ``run_fase_a_pdb_modeling``) --
esta validacion mide la CALIDAD DE LA ANOTACION/CONSENSO (Fase 2+3), no el modelado
estructural de un top-N de sitios ya aceptados (ese es un problema distinto, ya validado
por separado con corridas reales sobre Tau, ver STATUS.md). Evitar Fase A tambien mantiene
el tiempo de corrida razonable -- DeepPTMPred (ESM-2 + features PyRosetta) ya domina el
tiempo por si solo (~10+ min incluso para la proteina mas pequeña del panel, 103
residuos -- confirmado real 2026-08-04, no es proporcional al tamaño de la proteina).

Uso:
    conda run -n cnb_pipeline python scripts/validate_biological_panel.py [--only NOMBRE ...]

Requiere DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN apuntando a conda envs reales con los
motores instalados (ver README.md "Instalacion") -- no mockea nada, es una corrida real.
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
    filtered, _ = run_fase3_pdb_annotation(record, deepmvp_results, deepptmpred_results, output_dir)
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
    parser.add_argument("--output-dir", default="fasta_outputs/validation_panel")
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
