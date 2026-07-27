"""Orquestador del pipeline de prediccion de zonas PTM.

Estado actual (2026-07-27): Fase 1 (saneamiento FASTA) y Fase 1.5
(extraccion de estructura PDB/mmCIF) implementadas y funcionales. Fase 3
(anotacion/filtrado + logica de flujo, motores DeepMVP/DeepPTMPred) todavia
NO esta implementada -- ver
PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md
en el vault para el diseno ya decidido, pendiente de construir.

Camino FASTA: Fase 1 (este modulo) -> DeepMVP (unico motor, sin corroborar).
Camino PDB: Fase 1.5 (extrae FASTA canonico ATMSEQ de la estructura) ->
DeepMVP + DeepPTMPred en consenso.
"""

import argparse
import sys
from pathlib import Path
from typing import List

from src.config.settings import Settings
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


def main(argv: List[str] = None) -> int:
    args = parse_args(argv)
    Settings.ensure_dirs()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        routed = route_input(input_path)

        if routed.input_type == "fasta":
            run_fase1_fasta(input_path, output_dir)
            print(f"Camino FASTA: Fase 1 completa. Motor unico: DeepMVP (Fase 3 aun no implementada).")
        else:
            run_fase1_5_structure(input_path, output_dir)
            print(f"Camino PDB: Fase 1.5 completa. Consenso: DeepMVP + DeepPTMPred (Fase 3 aun no implementada).")

    except PipelineError as exc:
        logger.error("Pipeline detenido: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
