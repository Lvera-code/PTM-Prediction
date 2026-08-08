"""Orquestador del pipeline de prediccion de zonas PTM.

Numeracion de fases alineada con `BCell-Epitope-Prediction` (proyecto 1):
Fase 1 (saneamiento) -> Fase 1.5 (extraccion de estructura, Camino PDB
unicamente) -> Fase 2 (motores: DeepMVP / DeepMVP+DeepPTMPred) -> Fase 3
(nucleo: anotacion + filtro) -> Fase 3b (cruces informativos: via secretora,
Kinase Library, MeToken, competencia entre PTMs) -> Fase 3c (modelado
estructural real, Camino PDB unicamente). Renombrado 2026-08-08 -- Fase 3b/3c
antes se llamaban "cruces informativos"/"Fase A" en la documentacion; sin
cambio de comportamiento, ambas siguen calculandose exactamente igual
(3b dentro de ``annotate_pdb_path``, 3c via ``run_fase_a_pdb_modeling``/
``FaseAEngine``) -- el nombre nuevo es solo para el resumen en pantalla y la
documentacion, los nombres de funciones/columnas/variables de entorno
(``fase_a_*``, ``FASE_A_ENABLED``, etc.) no cambian.

Camino FASTA: Fase 1 (saneamiento) -> Fase 2 (DeepMVP, unico motor) -> Fase 3
(nucleo: anotacion + filtro).
Camino PDB: Fase 1.5 (extraccion ATMSEQ + pdb de una cadena) -> Fase 2
(DeepMVP + DeepPTMPred en consenso) -> Fase 3 (nucleo: anotacion + filtro,
fusiona consenso donde ambos motores coinciden en tipo+posicion) -> Fase 3b
(cruces informativos) -> Fase 3c (modelado estructural real).

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
from src.utils.input_router import FASTA_EXTENSIONS, STRUCTURE_EXTENSIONS, route_input

# Analisis de coherencia biologica 2026-08-07 (cambio 1 de 3): ningun motor
# de este pipeline (DeepMVP/DeepPTMPred/EMNGly/StackGlyEmbed) modela la
# via biosintetica real del sustrato (co-expresion/co-localizacion de la
# enzima) -- todos predicen CAPACIDAD de un sitio de modificarse a partir de
# secuencia/estructura, nunca si esa PTM ocurre realmente en una celula/
# tejido/condicion dada (eso depende de contexto biologico que ningun
# predictor de secuencia puede capturar). Impreso una vez al final de cada
# corrida (no se escribe dentro del CSV del reporte -- cambiaria su esquema,
# rompiendo lectores existentes que no esperan una fila de comentario).
INTERPRETATION_DISCLAIMER = (
    "AVISO DE INTERPRETACION: este reporte predice sitios POTENCIALMENTE "
    "modificables (propiedad de secuencia/estructura), no si esa "
    "modificacion ocurre realmente en una celula/tejido/condicion "
    "especifica -- eso depende de la co-expresion y co-localizacion de la "
    "enzima real, que este pipeline no modela. Ver README.md - seccion "
    "'Alcance e interpretacion'."
)
from src.utils.logger_config import setup_logger
from src.utils.structure_parser import parse_structure

logger = setup_logger(__name__)

_SEPARATOR = "=" * 70


def _print_table(headers: List[str], rows: List[list]) -> None:
    """Imprime una tabla alineada por ancho de columna con ``print()`` plano.

    Sin dependencias nuevas (no usa ``DataFrame.to_string``, que arrastra
    indice/dtypes que no aportan nada a un resumen en pantalla) -- mismo
    espiritu que los resumenes por fase de BCell-Epitope-Prediction, que
    tambien son ``print()`` simples, no volcados de pandas.
    """
    widths = [len(str(h)) for h in headers]
    str_rows = [[str(c) for c in row] for row in rows]
    for row in str_rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def _fmt(cells: List[str]) -> str:
        return "  ".join(cell.ljust(w) for cell, w in zip(cells, widths))

    print(_fmt([str(h) for h in headers]))
    print(_fmt(["-" * w for w in widths]))
    for row in str_rows:
        print(_fmt(row))


def _ground_truth_lookup(accession: str) -> dict:
    """Si ``accession`` coincide con una entrada de
    ``src/validation/biological_panel.py`` (panel de 7 proteinas con sitios
    PTM reales documentados en literatura), devuelve ``{(posicion, tipo_ptm):
    tier}`` para marcar en el resumen en pantalla que sitios de consenso ya
    estan confirmados por PMID real -- nunca inventa literatura donde no la
    hay: para cualquier proteina fuera del panel devuelve ``{}`` y el resumen
    simplemente no muestra la columna.
    """
    from src.validation.biological_panel import PANEL

    for entry in PANEL:
        if Path(entry.pdb_filename).stem == accession:
            return {
                (site.position, site.ptm_type): site.tier
                for site in entry.sites
                if not site.is_negative
            }
    return {}


def _fmt_score(value) -> str:
    return f"{float(value):.3f}" if pd.notna(value) else "-"


def parse_args(argv: List[str] = None) -> argparse.Namespace:
    """Define y parsea los argumentos de linea de comandos del pipeline."""
    parser = argparse.ArgumentParser(
        prog="pipeline.py",
        description="Pipeline de prediccion de zonas PTM (DeepMVP + DeepPTMPred).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True,
        help="Ruta a un archivo de entrada (dentro de inputs/): FASTA (Camino FASTA, "
        "DeepMVP solo) o PDB/mmCIF (Camino PDB, consenso DeepMVP+DeepPTMPred); o a un "
        "directorio con varios de estos archivos, para correr el pipeline sobre todos en "
        "una sola invocacion (modo batch). El tipo de cada archivo se detecta "
        "automaticamente (ver src.utils.input_router).",
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
) -> "tuple[pd.DataFrame, int, Path]":
    """Camino FASTA: Fase 3, nucleo (B: anotacion + D: filtro). Escribe el reporte final.

    DeepMVP procesa el FASTA completo (posiblemente multi-accession) en una
    unica invocacion (columna ``protein`` distingue cada accession en su
    salida): se anota cada accession por separado (secuencia correcta para
    el calculo de ventanas de cada una) y se concatena, en vez de tratar el
    FASTA como una unica secuencia fusionada.

    Devuelve ``(filtered, n_evaluados, report_path)`` -- ``n_evaluados`` es
    ``len(annotated)`` (antes del filtro de umbral), para que el resumen en
    pantalla pueda mostrar "X/Y pasan el umbral" sin releer el CSV.
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
    return filtered, len(annotated), report_path


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
) -> "tuple[pd.DataFrame, int, Path]":
    """Camino PDB: Fase 3, nucleo con consenso (B: anotacion + D: filtro).

    ``record.chain_pdb_path`` (no ``record.pdb_path``, que puede tener mas
    de una cadena) se pasa a ``annotate_pdb_path`` para habilitar la
    corroboracion opcional de tipo via MeToken -- sus posiciones 1-based
    deben coincidir exactamente con ``record.sequence``, mismo criterio que
    usa ``record.chain_id`` (ver Fase 1.5). Se activa/desactiva solo con
    ``Settings.METOKEN_ENABLED`` (ver ``src/engines/ptm_annotation.py``), sin
    tocar este orquestador. ``record.position_mapping`` (misma Fase 1.5) se
    pasa tambien para habilitar el consenso REAL de N-glicosilacion via
    EMNGly+StackGlyEmbed (``Settings.EMNGLY_ENABLED``, decision 2026-08-06,
    reemplaza a CoNglyPred -- ver STATUS.md) -- necesita traducir posiciones
    ATMSEQ a numeracion real de PDB para alinear ``structure_emb``
    correctamente (ver docstring de ``_emngly_runner.py``).

    Devuelve ``(filtered, n_evaluados, report_path)`` -- a diferencia de
    antes de 2026-08-03, expone tambien el DataFrame en memoria (no solo la
    ruta del CSV ya escrito) porque ``run_fase_a_pdb_modeling`` (paso
    siguiente en ``main()``) necesita seleccionar candidatos de ``filtered``
    sin tener que releerlo de disco; ``n_evaluados`` (``len(annotated)``,
    antes del filtro de umbral) es lo que el resumen en pantalla necesita
    para mostrar "X/Y pasan el umbral".
    """
    annotated = annotate_pdb_path(
        record.accession, record.sequence, deepmvp_results, deepptmpred_results,
        pdb_path=record.chain_pdb_path, chain_id=record.chain_id,
        enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED,
        position_mapping=record.position_mapping,
    )
    filtered = apply_workflow_filter(annotated)

    report_path = output_dir / f"{record.accession}_ptm_sites.csv"
    filtered.to_csv(report_path, index=False)
    n_consenso = int(annotated["consenso"].sum())
    logger.info(
        "Fase 3 completa (Camino PDB): %d/%d sitio(s) PTM pasan el umbral (%d con consenso -- "
        "DeepMVP+DeepPTMPred para la mayoria de tipos, DeepMVP+EMNGly+StackGlyEmbed para "
        "N-glicosilacion) -> '%s'.",
        len(filtered), len(annotated), n_consenso, report_path,
    )
    return filtered, len(annotated), report_path


_FASE_A_RESULT_COLUMNS = {
    "estado": "fase_a_estado",
    "clase": "fase_a_clase",
    "ddg": "fase_a_ddg",
    "ddg_std": "fase_a_ddg_std",
    "glycan_tree": "fase_a_glycan_tree",
    "glygen_evidencia": "fase_a_glygen_evidencia",
    "conjugation_metrics": "fase_a_conjugation_metrics",
    "cadena_tipo_aviso": "fase_a_cadena_tipo_aviso",
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
    ``fase_a_estado``/``fase_a_clase``/``fase_a_ddg``/``fase_a_ddg_std``/
    ``fase_a_glycan_tree``/``fase_a_glygen_evidencia``/``fase_a_conjugation_metrics``/
    ``fase_a_cadena_tipo_aviso``/``fase_a_output_pdb``
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
        logger.info("Fase 3c: ningun candidato seleccionable (0 sitios de los 9 tipos soportados).")
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
        "Fase 3c completa (Camino PDB): %d/%d sitio(s) candidato(s) modelados con exito -> '%s'.",
        n_modelado, len(candidates), report_path,
    )
    return enriched


def _discover_batch_inputs(directory: Path) -> List[Path]:
    """Lista, en orden alfabetico, los archivos de ``directory`` con una extension reconocida
    por ``src.utils.input_router`` (FASTA o estructura). No recursivo -- mismo alcance plano
    que ``inputs/`` en el uso de un solo archivo. El contenido de cada archivo se
    revalida igual que en modo de un solo archivo (``route_input`` dentro de ``run_single_input``),
    esto es solo un filtro rapido por extension para no intentar rutear cada archivo random
    que pueda haber en la carpeta (p.ej. ``.csv``/``.log`` de una corrida anterior).
    """
    recognized = FASTA_EXTENSIONS | STRUCTURE_EXTENSIONS
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in recognized
    )


def _print_fasta_summary(filtered: pd.DataFrame, n_evaluados: int, report_path: Path) -> None:
    """Resumen limpio del Camino FASTA -- sin motor de estructura, no hay consenso
    ni Fase 3b/3c: la unica pregunta real es cuantos sitios acepta DeepMVP y de que tipo.
    """
    print(f"\n{_SEPARATOR}\nFASE 3 | Anotacion + filtro (DeepMVP)\n{_SEPARATOR}")
    print(f"{len(filtered)}/{n_evaluados} sitio(s) pasan el umbral.")
    if len(filtered):
        counts = filtered["tipo_ptm"].value_counts()
        print("Por tipo: " + " . ".join(f"{tipo} {n}" for tipo, n in counts.items()))
    print(f"\nReporte completo: {report_path}")


def _print_pdb_summary(
    record, deepmvp_results: pd.DataFrame, deepptmpred_results: pd.DataFrame,
    filtered: pd.DataFrame, n_evaluados: int, enriched: pd.DataFrame, report_path: Path,
) -> None:
    """Resumen limpio del Camino PDB -- lo unico que se ve en pantalla tras una corrida.

    Todo el detalle de subprocess/progreso por-tipo que antes saturaba la
    consola (una linea INFO por cada uno de los 17 tipos de DeepPTMPred, por
    ejemplo) ya no llega aqui -- vive solo en ``logs/ptm_pipeline.log`` (ver
    ``src/utils/logger_config.py``). Lo que se muestra es exactamente lo que
    un usuario necesita para juzgar el resultado: cuanto paso, que tiene
    consenso real, y que se modelo estructuralmente -- nunca la mecanica
    interna de como se llego ahi.
    """
    print(f"\n{_SEPARATOR}\nFASE 1.5 | Extraccion de estructura\n{_SEPARATOR}")
    print(f"Cadena '{record.chain_id}' ({len(record.sequence)} residuo(s)) -> {record.chain_pdb_path}")

    print(f"\n{_SEPARATOR}\nFASE 2 | Motores (DeepMVP + DeepPTMPred)\n{_SEPARATOR}")
    print(
        f"DeepMVP: {len(deepmvp_results)} prediccion(es) crudas | "
        f"DeepPTMPred: {len(deepptmpred_results)} prediccion(es) crudas"
    )

    print(f"\n{_SEPARATOR}\nFASE 3 | Consenso + anotacion\n{_SEPARATOR}")
    n_consenso = int(filtered["consenso"].sum()) if len(filtered) else 0
    print(f"{len(filtered)}/{n_evaluados} sitio(s) pasan el umbral | {n_consenso} con consenso real (2+ motores de acuerdo)")
    if len(filtered):
        counts = filtered["tipo_ptm"].value_counts()
        print("Por tipo: " + " . ".join(f"{tipo} {n}" for tipo, n in counts.items()))

    ground_truth = _ground_truth_lookup(record.accession)
    consenso_df = filtered[filtered["consenso"] == True].sort_values("posicion") if len(filtered) else filtered  # noqa: E712
    if len(consenso_df):
        print(f"\n-- Consenso real ({len(consenso_df)} sitio(s)) --")
        headers = ["Pos", "Res", "Tipo", "Motor", "DeepMVP", "DeepPTMPred", "Avisos"]
        if ground_truth:
            headers.append("Literatura")
        rows = []
        for _, row in consenso_df.iterrows():
            avisos = []
            if row.get("metoken_type_coincide") is False:
                avisos.append("MeToken<>")
            if pd.notna(row.get("ptm_crosstalk_aviso")):
                avisos.append("crosstalk")
            cells = [
                int(row["posicion"]), row["residuo_wt"], row["tipo_ptm"], row["motor"],
                _fmt_score(row.get("score_deepmvp")), _fmt_score(row.get("score_deepptmpred")),
                ", ".join(avisos) if avisos else "-",
            ]
            if ground_truth:
                tier = ground_truth.get((int(row["posicion"]), row["tipo_ptm"]))
                cells.append(f"tier {tier}" if tier else "-")
            rows.append(cells)
        _print_table(headers, rows)
    else:
        print("\nNingun sitio con consenso real en esta corrida.")

    print(f"\n{_SEPARATOR}\nFASE 3b | Cruces informativos\n{_SEPARATOR}")
    if len(filtered):
        nglyco = filtered[filtered["tipo_ptm"].isin(["n_linked_glycosylation", "glycosylation_n"])]
        if len(nglyco):
            via_sec = nglyco["via_secretora_evidencia"]
            print(
                f"Via secretora (N-glico, UniProt): {int((via_sec == True).sum())} con evidencia, "  # noqa: E712
                f"{int((via_sec == False).sum())} sin evidencia, {int(via_sec.isna().sum())} sin dato."  # noqa: E712
            )
        fosfo = filtered[filtered["tipo_ptm"] == "phosphorylation"]
        if len(fosfo):
            con_kinasa = int(fosfo["kinase_library_top_kinase"].notna().sum())
            print(f"Kinase Library: {con_kinasa}/{len(fosfo)} fosforilacion(es) con familia de quinasa asignada.")
        con_metoken = int(filtered["metoken_type"].notna().sum())
        if con_metoken:
            desacuerdo = int((filtered["metoken_type_coincide"] == False).sum())  # noqa: E712
            print(f"MeToken: {con_metoken} sitio(s) corroborado(s), {desacuerdo} en desacuerdo con el tipo de consenso.")
        con_crosstalk = int(filtered["ptm_crosstalk_aviso"].notna().sum())
        if con_crosstalk:
            print(f"Competencia entre PTMs: {con_crosstalk} sitio(s) con aviso de crosstalk real (ver columna 'Avisos' arriba).")
    else:
        print("Sin sitios que evaluar.")

    print(f"\n{_SEPARATOR}\nFASE 3c | Modelado estructural real\n{_SEPARATOR}")
    candidatos = enriched[enriched["fase_a_estado"] != "no_seleccionado"]
    if candidatos.empty:
        print("Ningun candidato seleccionable (0 sitios de los 9 tipos con modulo real).")
    elif (candidatos["fase_a_estado"] == "no_disponible").all():
        print(f"Desactivada para esta corrida (FASE_A_ENABLED=false) -- {len(candidatos)} candidato(s) habrian sido seleccionados.")
    else:
        n_modelado = int((candidatos["fase_a_estado"] == "modelado").sum())
        print(f"{n_modelado}/{len(candidatos)} candidato(s) modelados con exito.")
        headers = ["Pos", "Tipo", "Estado", "Detalle"]
        rows = []
        for _, row in candidatos.sort_values("posicion").iterrows():
            detalle = "-"
            if row["fase_a_estado"] == "modelado":
                if pd.notna(row.get("fase_a_ddg")):
                    detalle = f"ddG={float(row['fase_a_ddg']):.2f}"
                elif pd.notna(row.get("fase_a_glycan_tree")):
                    detalle = "glicano adjuntado"
                elif pd.notna(row.get("fase_a_conjugation_metrics")):
                    detalle = "conjugacion modelada"
            elif row["fase_a_estado"] == "error":
                detalle = "fallo (ver log)"
            rows.append([int(row["posicion"]), row["tipo_ptm"], row["fase_a_estado"], detalle])
        _print_table(headers, rows)

    print(f"\nReporte completo: {report_path}")


def run_single_input(input_path: Path, output_dir: Path) -> Path:
    """Corre el pipeline completo (Camino FASTA o PDB, segun ``route_input``) sobre un unico
    archivo de entrada. Extraido de ``main()`` para que el modo batch (``main()`` con
    ``--input`` apuntando a un directorio) pueda reusar exactamente la misma logica por
    archivo, sin duplicarla.

    Devuelve la ruta del reporte final. Cualquier ``PipelineError`` de las fases internas se
    propaga sin capturar -- el llamador decide si eso detiene todo el proceso (modo de un
    solo archivo) o solo ese archivo del batch (modo batch, ver ``main()``).
    """
    routed = route_input(input_path)

    if routed.input_type == "fasta":
        records = load_and_sanitize(input_path)
        out_stem = records[0].accession if len(records) == 1 else input_path.stem
        clean_path = run_fase1_fasta(input_path, output_dir)
        deepmvp_results = run_fase2_fasta_motor(clean_path, output_dir)
        filtered, n_evaluados, report_path = run_fase3_fasta_annotation(
            records, deepmvp_results, output_dir, out_stem
        )
        _print_fasta_summary(filtered, n_evaluados, report_path)
    else:
        record = run_fase1_5_structure(input_path, output_dir)
        deepmvp_results, deepptmpred_results = run_fase2_pdb_motors(record, output_dir)
        filtered, n_evaluados, report_path = run_fase3_pdb_annotation(
            record, deepmvp_results, deepptmpred_results, output_dir
        )
        enriched = run_fase_a_pdb_modeling(record, filtered, output_dir, report_path)
        _print_pdb_summary(
            record, deepmvp_results, deepptmpred_results, filtered, n_evaluados, enriched, report_path
        )

    print(f"\n{INTERPRETATION_DISCLAIMER}")
    return report_path


def _run_batch(input_dir: Path, output_dir: Path) -> int:
    """Modo batch: corre ``run_single_input`` sobre cada archivo reconocido de ``input_dir``.

    Un archivo que falla se registra como error y NO detiene el resto del batch -- mismo
    criterio de degradacion no fatal que ``FaseAEngine``/StackGlyEmbed/MeToken aplican
    por-sitio (un fallo individual real no debe tumbar todo el barrido). Escribe
    ``batch_summary.csv`` (columnas: archivo, estado, reporte/error) en ``output_dir``.
    Codigo de salida: 0 solo si TODOS los archivos completaron sin error, 1 si al menos uno
    fallo (incluye el caso de un directorio sin ningun archivo reconocido -- ver docstring de
    ``_discover_batch_inputs``, "vacio" tambien es un resultado que no debe pasar en silencio).

    Captura ``Exception`` en general, no solo ``PipelineError`` (bug real encontrado en
    auditoria 2026-08-07: antes, cualquier excepcion inesperada -- p.ej. un ``KeyError``/
    ``AttributeError`` real de un motor/engine con una salida malformada, no necesariamente
    un ``PipelineError`` -- escapaba del bucle entero, tumbando TODO el batch antes de
    escribir ``batch_summary.csv`` y perdiendo el registro de los archivos previos ya
    procesados con exito, pese a que sus reportes individuales SI quedaban en disco). Mismo
    patron ``# noqa: BLE001`` ya usado en ``fase_a_dispatch.run_fase_a_for_site`` para el
    mismo motivo (un item individual de un barrido no debe tumbar el resto).
    """
    inputs = _discover_batch_inputs(input_dir)
    if not inputs:
        logger.error("Modo batch: '%s' no contiene ningun archivo FASTA/PDB/mmCIF reconocido.", input_dir)
        print(f"ERROR: '{input_dir}' no contiene ningun archivo de entrada reconocido.", file=sys.stderr)
        return 1

    logger.info("Modo batch: %d archivo(s) encontrado(s) en '%s'.", len(inputs), input_dir)

    rows = []
    for input_path in inputs:
        try:
            report_path = run_single_input(input_path, output_dir)
            rows.append({"archivo": input_path.name, "estado": "ok", "detalle": str(report_path)})
        except Exception as exc:  # noqa: BLE001 -- un archivo individual no debe tumbar el resto del batch
            logger.error("Modo batch: '%s' fallo, continua con el resto -- %s", input_path.name, exc)
            rows.append({
                "archivo": input_path.name, "estado": "error",
                "detalle": f"{type(exc).__name__}: {exc}",
            })

    summary_path = output_dir / "batch_summary.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)

    n_ok = sum(1 for row in rows if row["estado"] == "ok")
    print(f"Modo batch completo: {n_ok}/{len(rows)} archivo(s) OK. Resumen: {summary_path}")
    return 0 if n_ok == len(rows) else 1


def main(argv: List[str] = None) -> int:
    args = parse_args(argv)
    Settings.ensure_dirs()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        return _run_batch(input_path, output_dir)

    try:
        run_single_input(input_path, output_dir)
    except PipelineError as exc:
        logger.error("Pipeline detenido: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
