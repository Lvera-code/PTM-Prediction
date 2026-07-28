"""Orquestador del pipeline de prediccion de zonas PTM.

Numeracion de fases alineada con `BCell-Epitope-Prediction` (proyecto 1):
Fase 1 (saneamiento) -> Fase 1.5 (extraccion de estructura, Camino PDB
unicamente) -> Fase 2 (motores: DeepMVP / DeepMVP+DeepPTMPred) -> Fase 3
(nucleo: anotacion + filtro). Renombrado 2026-07-28 -- antes los motores se
llamaban "Fase 3a" en vez de "Fase 2", inconsistente con el numerado del
proyecto 1 sin que hubiera una decision documentada detras; ningun cambio
de comportamiento, solo de nombres (funciones, docstrings, logs).

Camino FASTA: Fase 1 (saneamiento) -> Fase 2 (DeepMVP, unico motor) -> Fase 3
(nucleo: anotacion + filtro).
Camino PDB: Fase 1.5 (extraccion ATMSEQ + pdb de una cadena) -> Fase 2
(DeepMVP + DeepPTMPred en consenso) -> Fase 3 (nucleo: anotacion + filtro,
fusiona consenso donde ambos motores coinciden en tipo+posicion).

Diseno del nucleo de Fase 3 (B: anotacion/filtrado, D: logica de flujo) en
``01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md``
del vault.
"""

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd

from src.config.settings import Settings
from src.engines.deepmvp_engine import DeepMVPEngine
from src.engines.deepptmpred_engine import DeepPTMPredEngine
from src.engines.ptm_annotation import annotate_fasta_path, annotate_pdb_path, apply_workflow_filter
from src.utils.exceptions import PipelineError
from src.utils.fasta_parser import load_and_sanitize, write_fasta
from src.utils.input_router import route_input
from src.utils.logger_config import setup_logger
from src.utils.structure_parser import parse_structure

logger = setup_logger(__name__)


def parse_args(argv: List[str] = None) -> argparse.Namespace:
    """Define y parsea los argumentos de linea de comandos del pipeline."""
    parser = argparse.ArgumentParser(
        prog="pipeline.py",
        description="Pipeline de prediccion de zonas PTM (DeepMVP + DeepPTMPred).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True,
        help="Ruta al archivo de entrada (dentro de fasta_inputs/): FASTA (Camino FASTA, "
        "DeepMVP solo) o PDB/mmCIF (Camino PDB, consenso DeepMVP+DeepPTMPred). El tipo se "
        "detecta automaticamente (ver src.utils.input_router).",
    )
    parser.add_argument(
        "--output-dir", default=str(Settings.FASTA_OUTPUT_DIR),
        help="Carpeta donde se guardan todos los resultados del pipeline.",
    )
    return parser.parse_args(argv)


def run_fase1_fasta(input_path: Path, output_dir: Path) -> Path:
    """Camino FASTA: satura Fase 1 y escribe el FASTA saneado listo para DeepMVP."""
    records = load_and_sanitize(input_path)
    accession = records[0].accession if len(records) == 1 else input_path.stem
    clean_path = output_dir / f"{accession}_clean.fasta"
    write_fasta(records, clean_path)
    logger.info(
        "Fase 1 completa (Camino FASTA): %d registro(s) saneado(s) -> '%s'.",
        len(records), clean_path,
    )
    return clean_path


def run_fase1_5_structure(input_path: Path, output_dir: Path):
    """Camino PDB: corre Fase 1.5 y devuelve el StructureRecord (FASTA ATMSEQ + pdb_path)."""
    record = parse_structure(input_path, output_dir)
    logger.info(
        "Fase 1.5 completa (Camino PDB): accession='%s', %d residuo(s). "
        "FASTA derivado ('%s') alimenta DeepMVP; pdb_path original ('%s') alimenta DeepPTMPred.",
        record.accession, len(record.sequence), record.fasta_path, record.pdb_path,
    )
    return record


def run_fase2_fasta_motor(clean_fasta: Path, output_dir: Path) -> pd.DataFrame:
    """Camino FASTA: Fase 2, unico motor (DeepMVP). Devuelve las predicciones crudas."""
    deepmvp_results = DeepMVPEngine().run([str(clean_fasta)], output_dir=output_dir)[0]
    logger.info(
        "Fase 2 completa (Camino FASTA): %d prediccion(es) crudas de DeepMVP.",
        len(deepmvp_results),
    )
    return deepmvp_results


def run_fase3_fasta_annotation(
    records, deepmvp_results: pd.DataFrame, output_dir: Path, out_stem: str
) -> Path:
    """Camino FASTA: Fase 3, nucleo (B: anotacion + D: filtro). Escribe el reporte final.

    DeepMVP procesa el FASTA completo (posiblemente multi-accession) en una
    unica invocacion (columna ``protein`` distingue cada accession en su
    salida): se anota cada accession por separado (secuencia correcta para
    el calculo de ventanas de cada una) y se concatena, en vez de tratar el
    FASTA como una unica secuencia fusionada.
    """
    per_accession = []
    for record in records:
        subset = deepmvp_results[deepmvp_results["protein"] == record.accession]
        per_accession.append(annotate_fasta_path(record.accession, record.sequence, subset))
    annotated = pd.concat(per_accession, ignore_index=True) if per_accession else deepmvp_results
    filtered = apply_workflow_filter(annotated)

    report_path = output_dir / f"{out_stem}_ptm_sites.csv"
    filtered.to_csv(report_path, index=False)
    logger.info(
        "Fase 3 completa (Camino FASTA): %d/%d sitio(s) PTM pasan el umbral -> '%s'.",
        len(filtered), len(annotated), report_path,
    )
    return report_path


def run_fase2_pdb_motors(record, output_dir: Path):
    """Camino PDB: Fase 2, DeepMVP + DeepPTMPred. Devuelve ambos DataFrames crudos."""
    deepmvp_results = DeepMVPEngine().run([str(record.fasta_path)], output_dir=output_dir)[0]
    deepptmpred_results = DeepPTMPredEngine().run([record], output_dir=output_dir)[0]
    logger.info(
        "Fase 2 completa (Camino PDB): %d prediccion(es) crudas de DeepMVP, %d de DeepPTMPred.",
        len(deepmvp_results), len(deepptmpred_results),
    )
    return deepmvp_results, deepptmpred_results


def run_fase3_pdb_annotation(
    record, deepmvp_results: pd.DataFrame, deepptmpred_results: pd.DataFrame, output_dir: Path
) -> Path:
    """Camino PDB: Fase 3, nucleo con consenso (B: anotacion + D: filtro)."""
    annotated = annotate_pdb_path(record.accession, record.sequence, deepmvp_results, deepptmpred_results)
    filtered = apply_workflow_filter(annotated)

    report_path = output_dir / f"{record.accession}_ptm_sites.csv"
    filtered.to_csv(report_path, index=False)
    n_consenso = int(annotated["consenso"].sum())
    logger.info(
        "Fase 3 completa (Camino PDB): %d/%d sitio(s) PTM pasan el umbral (%d con consenso "
        "DeepMVP+DeepPTMPred) -> '%s'.",
        len(filtered), len(annotated), n_consenso, report_path,
    )
    return report_path


def main(argv: List[str] = None) -> int:
    args = parse_args(argv)
    Settings.ensure_dirs()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        routed = route_input(input_path)

        if routed.input_type == "fasta":
            records = load_and_sanitize(input_path)
            out_stem = records[0].accession if len(records) == 1 else input_path.stem
            clean_path = run_fase1_fasta(input_path, output_dir)
            deepmvp_results = run_fase2_fasta_motor(clean_path, output_dir)
            report_path = run_fase3_fasta_annotation(records, deepmvp_results, output_dir, out_stem)
            print(f"Camino FASTA completo. Motor: DeepMVP. Reporte: {report_path}")
        else:
            record = run_fase1_5_structure(input_path, output_dir)
            deepmvp_results, deepptmpred_results = run_fase2_pdb_motors(record, output_dir)
            report_path = run_fase3_pdb_annotation(record, deepmvp_results, deepptmpred_results, output_dir)
            print(f"Camino PDB completo. Consenso DeepMVP+DeepPTMPred. Reporte: {report_path}")

    except PipelineError as exc:
        logger.error("Pipeline detenido: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
