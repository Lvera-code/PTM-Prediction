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
import json
import sys
from pathlib import Path
from typing import List

import pandas as pd

from src.config.settings import Settings
from src.engines.deepmvp_engine import DeepMVPEngine
from src.engines.deepptmpred_engine import DeepPTMPredEngine
from src.engines.fase_a_engine import FaseAEngine, FaseASiteRequest
from src.engines.ptm_annotation import (
    annotate_fasta_path,
    annotate_pdb_path,
    apply_workflow_filter,
    select_fase_a_candidates,
)
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
        per_accession.append(
            annotate_fasta_path(
                record.accession, record.sequence, subset,
                enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED,
            )
        )
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
) -> "tuple[pd.DataFrame, Path]":
    """Camino PDB: Fase 3, nucleo con consenso (B: anotacion + D: filtro).

    ``record.chain_pdb_path`` (no ``record.pdb_path``, que puede tener mas
    de una cadena) se pasa a ``annotate_pdb_path`` para habilitar la
    corroboracion opcional de tipo via MeToken -- sus posiciones 1-based
    deben coincidir exactamente con ``record.sequence``, mismo criterio que
    usa ``record.chain_id`` (ver Fase 1.5). Se activa/desactiva solo con
    ``Settings.METOKEN_ENABLED`` (ver ``src/engines/ptm_annotation.py``), sin
    tocar este orquestador. La corroboracion de N-glicosilacion via
    StackGlyEmbed (``Settings.STACKGLYEMBED_ENABLED``) solo necesita
    ``record.sequence`` -- se activa igual, sin depender del PDB.

    Devuelve ``(filtered, report_path)`` -- a diferencia de antes de
    2026-08-03, expone tambien el DataFrame en memoria (no solo la ruta del
    CSV ya escrito) porque ``run_fase_a_pdb_modeling`` (paso siguiente en
    ``main()``) necesita seleccionar candidatos de ``filtered`` sin tener que
    releerlo de disco.
    """
    annotated = annotate_pdb_path(
        record.accession, record.sequence, deepmvp_results, deepptmpred_results,
        pdb_path=record.chain_pdb_path, chain_id=record.chain_id,
        enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED,
    )
    filtered = apply_workflow_filter(annotated)

    report_path = output_dir / f"{record.accession}_ptm_sites.csv"
    filtered.to_csv(report_path, index=False)
    n_consenso = int(annotated["consenso"].sum())
    logger.info(
        "Fase 3 completa (Camino PDB): %d/%d sitio(s) PTM pasan el umbral (%d con consenso "
        "DeepMVP+DeepPTMPred) -> '%s'.",
        len(filtered), len(annotated), n_consenso, report_path,
    )
    return filtered, report_path


_FASE_A_RESULT_COLUMNS = {
    "estado": "fase_a_estado",
    "clase": "fase_a_clase",
    "ddg": "fase_a_ddg",
    "glycan_tree": "fase_a_glycan_tree",
    "glygen_evidencia": "fase_a_glygen_evidencia",
    "conjugation_metrics": "fase_a_conjugation_metrics",
    "output_pdb": "fase_a_output_pdb",
}


def run_fase_a_pdb_modeling(
    record, filtered: pd.DataFrame, output_dir: Path, report_path: Path
) -> pd.DataFrame:
    """Camino PDB: Fase A, modelado estructural real de un top-N de sitios por tipo.

    Conectado al pipeline principal 2026-08-03 (ver
    ``src/engines/fase_a_engine.py`` y ``src/structural/fase_a_dispatch.py``
    para el detalle completo): revierte la decision 2026-07-27 de que D no
    rutea a Extension 3/Fase A porque esas fases no existian todavia.

    Selecciona candidatos con ``select_fase_a_candidates`` (top-N por tipo
    entre los 9/17 tipos con modulo de Fase A real, nunca todos los sitios
    aceptados -- costo computacional real, ver docstring de la funcion),
    modela cada uno via ``FaseAEngine`` (subprocess con PyRosetta, un sitio
    por proceso) y reescribe ``report_path`` con las columnas
    ``fase_a_estado``/``fase_a_clase``/``fase_a_ddg``/``fase_a_glycan_tree``/
    ``fase_a_glygen_evidencia``/``fase_a_conjugation_metrics``/``fase_a_output_pdb``
    anadidas para TODAS las filas de ``filtered`` (no solo las seleccionadas):
    las no seleccionadas quedan con ``fase_a_estado="no_seleccionado"``, para
    que el reporte final documente explicitamente por que un sitio aceptado
    no tiene modelado estructural, en vez de dejar una columna vacia
    ambigua.

    Si ``Settings.FASE_A_ENABLED`` es ``False`` o no hay candidatos, escribe
    igualmente las columnas (todas ``no_seleccionado``/``no_disponible``) para
    que el esquema del reporte final sea estable independientemente de la
    configuracion.
    """
    enriched = filtered.copy()
    for column in _FASE_A_RESULT_COLUMNS.values():
        enriched[column] = None
    enriched["fase_a_estado"] = "no_seleccionado"

    candidates = select_fase_a_candidates(filtered, Settings.FASE_A_TOP_N_PER_TYPE)
    if candidates.empty:
        logger.info("Fase A: ningun candidato seleccionable (0 sitios de los 9 tipos soportados).")
        enriched.to_csv(report_path, index=False)
        return enriched

    requests = [
        FaseASiteRequest(
            accession=record.accession,
            pdb_path=record.chain_pdb_path,
            position=int(row["posicion"]),
            ptm_type=row["tipo_ptm"],
        )
        for _, row in candidates.iterrows()
    ]
    results = FaseAEngine().run(requests, output_dir=output_dir)

    for request, result in zip(requests, results):
        mask = (enriched["posicion"] == request.position) & (enriched["tipo_ptm"] == request.ptm_type)
        for result_key, column in _FASE_A_RESULT_COLUMNS.items():
            value = result.get(result_key)
            if value is not None and not isinstance(value, (str, int, float, bool)):
                value = json.dumps(value)
            enriched.loc[mask, column] = value

    enriched.to_csv(report_path, index=False)
    n_modelado = int((enriched["fase_a_estado"] == "modelado").sum())
    logger.info(
        "Fase A completa (Camino PDB): %d/%d sitio(s) candidato(s) modelados con exito -> '%s'.",
        n_modelado, len(candidates), report_path,
    )
    return enriched


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
            filtered, report_path = run_fase3_pdb_annotation(
                record, deepmvp_results, deepptmpred_results, output_dir
            )
            run_fase_a_pdb_modeling(record, filtered, output_dir, report_path)
            print(f"Camino PDB completo. Consenso DeepMVP+DeepPTMPred + Fase A. Reporte: {report_path}")

    except PipelineError as exc:
        logger.error("Pipeline detenido: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
