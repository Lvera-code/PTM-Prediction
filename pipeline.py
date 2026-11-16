"""Orchestrator of the PTM zone prediction pipeline.

Phase numbering aligned with `BCell-Epitope-Prediction` (project 1):
Phase 1 (sanitation) -> Phase 1.5 (structure extraction, PDB path
only) -> Phase 2 (engines: DeepMVP / DeepMVP+DeepPTMPred) -> Phase 3
(core: annotation + filter) -> Phase 3b (informative cross-references: secretory
pathway, Kinase Library, MeToken, competition between PTMs).

FASTA path: Phase 1 (sanitation) -> Phase 2 (DeepMVP, the only engine) -> Phase 3
(core: annotation + filter).
PDB path: Phase 1.5 (ATMSEQ extraction + single-chain pdb) -> Phase 2
(DeepMVP + DeepPTMPred in consensus) -> Phase 3 (core: annotation + filter,
merging consensus where both engines agree on type+position) -> Phase 3b
(informative cross-references).

The design of the Phase 3 core (B: annotation/filtering, D: flow logic) is in
``01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md``
in the vault.

Phase 3c (real structural modelling via PyRosetta, formerly "Phase A") was
removed from the project's scope by a supervision decision (feedback
2026-08-10, see ``01-Proyectos/PTM-Prediction/Decisiones/`` in the vault, the note
dated 2026-08-10). The core (Phases 1-3b) is confirmed feature-complete and the project
moves on to integration into Scipion.
"""

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd

from src.config.settings import Settings
from src.engines.deepmvp_engine import DeepMVPEngine
from src.engines.deepptmpred_engine import DeepPTMPredEngine
from src.engines.ptm_annotation import (
    annotate_fasta_path,
    annotate_pdb_path,
    apply_workflow_filter,
)
from src.utils.exceptions import PipelineError
from src.utils.fasta_parser import load_and_sanitize, write_fasta
from src.utils.input_router import FASTA_EXTENSIONS, STRUCTURE_EXTENSIONS, route_input

# Biological coherence analysis 2026-08-07 (change 1 of 3): no engine
# of this pipeline (DeepMVP/DeepPTMPred/EMNGly) models the
# substrate's real biosynthetic pathway (co-expression/co-localization of the
# enzyme). They all predict a site's CAPACITY to be modified from
# sequence/structure, never whether that PTM actually occurs in a given cell/
# tissue/condition (which depends on biological context that no
# sequence predictor can capture). Printed once at the end of each
# run (it is not written inside the report CSV, as that would change its schema and
# break existing readers that do not expect a comment row).
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
    """Prints a column-width-aligned table with plain ``print()``.

    Without new dependencies (it does not use ``DataFrame.to_string``, which carries
    an index/dtypes that add nothing to an on-screen summary), in the same
    spirit as the per-phase summaries of BCell-Epitope-Prediction, which
    are also simple ``print()`` calls rather than pandas dumps.
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
    """``{(posicion, tipo_ptm): (tier, fuente, pmids)}``, used to mark in the
    on-screen summary which consensus sites are already confirmed by a real PMID.
    It never invents literature where there is none.

    Two sources, in order of precedence (decision 2026-08-09, vault):
    1. ``src/validation/biological_panel.py::PANEL``, the small subset
       (~20-30 proteins) curated manually: every PMID verified individually
       against NCBI eutils, and the only ones that additionally run the complete
       pipeline to measure recall (``scripts/validate_biological_panel.py``). Marked
       ``fuente="curado"``.
    2. ``src/validation/dbptm_lookup.py``, a large database derived from dbPTM
       (thousands of reviewed human proteins), lookup only, never running
       engines. Marked ``fuente="dbPTM"``, which is a WEAKER signal than
       "curado": the tier there is only a count of PMIDs without individual
       verification, and the real verification of the curated subset
       found non-existent or misattributed PMIDs in the dbPTM draft
       (see the decision of 2026-08-09), so a "tier A" from this source does not carry
       the same weight as a "curado" one.
    For any protein outside both it returns ``{}`` and the summary
    simply does not show the columns.
    """
    from src.validation.biological_panel import PANEL
    from src.validation.dbptm_lookup import lookup_ground_truth

    for entry in PANEL:
        if Path(entry.pdb_filename).stem == accession:
            return {
                (site.position, site.ptm_type): (site.tier, "curado", site.pmids)
                for site in entry.sites
                if not site.is_negative
            }
    return {
        key: (tier, "dbPTM", pmids)
        for key, (tier, pmids) in lookup_ground_truth(accession).items()
    }


def _fmt_pmids(pmids: tuple, max_shown: int = 1) -> str:
    """``"20123963+7"``: the first ``max_shown`` PMIDs plus a count of the rest.

    Default ``max_shown=1`` (decision 2026-08-10): it was verified that in the
    curated panel the PMIDs are almost always sorted ascending (~92% of the 484
    sites, oldest ``pmids`` first; see ``import_dbptm_panel.py`` for
    the dbPTM fallback, which ALWAYS sorts this way), so the first one is a
    reasonable citation ("the study that reported it first", not necessarily "the most
    cited"). The remaining ~8% (5 proteins verified by a subagent that
    deliberately placed the dedicated/landmark paper first) is also
    well represented by the first one. It is never "the most important" in a
    strict sense, but it is a defensible citation to show on its own.
    """
    shown = ",".join(str(p) for p in pmids[:max_shown])
    resto = len(pmids) - max_shown
    return f"{shown}+{resto}" if resto > 0 else shown


def _crosstalk_count(aviso: str) -> int:
    """The number of competing types in ``ptm_crosstalk_aviso`` (the fixed format of
    ``_add_ptm_crosstalk_flag``: ``"Compite con: tipo1, tipo2 (...)"``).

    It only counts, and names none as the "main" one: it was verified in the code
    (``_add_ptm_crosstalk_flag``, ``sorted(set(...))``) that the list is
    sorted ALPHABETICALLY, not by biological relevance, so showing the
    first one would suggest a bias that does not exist (for example "acetylation" wins
    almost always by alphabetical order, not because it is the dominant competitor).
    Decision 2026-08-10.
    """
    marker = "Compite con: "
    if marker not in aviso:
        return 0
    return len(aviso.split(marker, 1)[1].split(" (mismo residuo", 1)[0].split(", "))


def _abbrev_motor(motor: str) -> str:
    """For on-screen display only. NEVER touch the DataFrame's real ``motor``
    column, which ``ptm_annotation.py`` uses as logic (for example
    ``result["motor"] == "DeepMVP"`` for nglyco rows, line 836)."""
    return motor.replace("DeepPTMPred", "DPTMP").replace("DeepMVP", "DMVP")


def _fmt_score(value) -> str:
    return f"{float(value):.3f}" if pd.notna(value) else "-"


def parse_args(argv: List[str] = None) -> argparse.Namespace:
    """Defines and parses the pipeline's command line arguments."""
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
    """FASTA path: completes Phase 1 and writes the sanitized FASTA ready for DeepMVP."""
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
    """FASTA path: Phase 3, the core (B: annotation + D: filter). It writes the final report.

    DeepMVP processes the complete FASTA (possibly multi-accession) in a
    single invocation (the ``protein`` column distinguishes each accession in its
    output), so each accession is annotated separately (with the correct sequence for
    computing each one's windows) and then concatenated, instead of treating the
    FASTA as a single fused sequence.

    It returns ``(filtered, n_evaluados, report_path)``, where ``n_evaluados`` is
    ``len(annotated)`` (before the threshold filter), so that the on-screen
    summary can show "X/Y pass the threshold" without re-reading the CSV.
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
    """PDB path: Phase 3, the core with consensus (B: annotation + D: filter).

    ``record.chain_pdb_path`` (not ``record.pdb_path``, which may have more
    than one chain) is passed to ``annotate_pdb_path`` in order to enable the
    optional type corroboration via MeToken: its 1-based positions
    must match ``record.sequence`` exactly, the same criterion that
    ``record.chain_id`` uses (see Phase 1.5). It is enabled/disabled solely with
    ``Settings.METOKEN_ENABLED`` (see ``src/engines/ptm_annotation.py``), without
    touching this orchestrator. ``record.position_mapping`` (the same Phase 1.5) is
    also passed in order to enable the REAL N-glycosylation consensus via
    EMNGly (``Settings.EMNGLY_ENABLED``, decision 2026-08-06, replacing
    CoNglyPred): it needs to translate ATMSEQ positions to the
    real PDB numbering in order to align ``structure_emb`` correctly
    (see the docstring of ``_emngly_runner.py``).

    It returns ``(filtered, n_evaluados, report_path)``, and also exposes the
    in-memory DataFrame (not only the path of the CSV already written) because
    ``_print_pdb_summary`` reads it directly for the on-screen summary,
    without re-reading it from disk; ``n_evaluados`` (``len(annotated)``, before the
    threshold filter) is what that summary needs in order to show "X/Y pass
    the threshold".
    """
    annotated = annotate_pdb_path(
        record.accession, record.sequence, deepmvp_results, deepptmpred_results,
        pdb_path=record.chain_pdb_path, chain_id=record.chain_id,
        position_mapping=record.position_mapping,
    )
    filtered = apply_workflow_filter(annotated)

    report_path = output_dir / f"{record.accession}_ptm_sites.csv"
    filtered.to_csv(report_path, index=False)
    n_consenso = int(annotated["consenso"].sum())
    logger.info(
        "Fase 3 completa (Camino PDB): %d/%d sitio(s) PTM pasan el umbral (%d con consenso -- "
        "DeepMVP+DeepPTMPred para la mayoria de tipos, DeepMVP+EMNGly para "
        "N-glicosilacion) -> '%s'.",
        len(filtered), len(annotated), n_consenso, report_path,
    )
    return filtered, len(annotated), report_path


def _discover_batch_inputs(directory: Path) -> List[Path]:
    """Lists, in alphabetical order, the files of ``directory`` with an extension recognized
    by ``src.utils.input_router`` (FASTA or structure). Not recursive, the same flat scope
    as ``inputs/`` in single-file use. The content of each file is
    revalidated exactly as in single-file mode (``route_input`` inside ``run_single_input``);
    this is only a quick filter by extension so as not to attempt routing every arbitrary file
    that may be in the folder (for example a ``.csv``/``.log`` from a previous run).
    """
    recognized = FASTA_EXTENSIONS | STRUCTURE_EXTENSIONS
    return sorted(
        p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in recognized
    )


def _print_fasta_summary(filtered: pd.DataFrame, n_evaluados: int, report_path: Path) -> None:
    """A clean summary of the FASTA path. Without a structural engine there is no consensus
    and no Phase 3b/3c: the only real question is how many sites DeepMVP accepts and of what type.
    """
    print(f"\n{_SEPARATOR}\nFASE 3 | Anotacion + filtro (DeepMVP)\n{_SEPARATOR}")
    print(f"{len(filtered)}/{n_evaluados} sitio(s) pasan el umbral.")
    if len(filtered):
        counts = filtered["tipo_ptm"].value_counts()
        print("Por tipo: " + " . ".join(f"{tipo} {n}" for tipo, n in counts.items()))
    print(f"\nReporte completo: {report_path}")


def _print_pdb_summary(
    record, deepmvp_results: pd.DataFrame, deepptmpred_results: pd.DataFrame,
    filtered: pd.DataFrame, n_evaluados: int, report_path: Path,
) -> None:
    """A clean summary of the PDB path, the only thing shown on screen after a run.

    All the per-type subprocess/progress detail that used to saturate the
    console (one INFO line for each of DeepPTMPred's 17 types, for
    example) no longer reaches here: it lives only in ``logs/ptm_pipeline.log`` (see
    ``src/utils/logger_config.py``). What is shown is exactly what
    a user needs in order to judge the result: how much passed, what has
    a real consensus, and what was modelled structurally, never the internal
    mechanics of how that was reached.
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
        if ground_truth:
            # The source is constant per accession, as _ground_truth_lookup never
            # mixes PANEL (curated) and the dbPTM fallback in the same run, so
            # repeating it on every row only wastes width. It is printed once.
            fuente = next(iter(ground_truth.values()))[1]
            if fuente == "curado":
                print("Literatura: fuente = panel curado (cada PMID verificado a mano contra NCBI eutils).")
            else:
                print(
                    "Literatura: fuente = base dbPTM (tier por conteo de PMIDs, NUNCA verificado "
                    "individualmente -- tratar como señal mas debil que el panel curado)."
                )
        headers = ["Pos", "Res", "Tipo", "Motor", "DMVP", "DPTMP", "Avisos"]
        mostrar_quinasa = (
            (consenso_df["tipo_ptm"] == "phosphorylation").any()
            and consenso_df["kinase_library_top_kinase"].notna().any()
        )
        if mostrar_quinasa:
            headers.append("Quinasa")
        if ground_truth:
            headers.extend(["Literatura", "PMIDs"])
        rows = []
        for _, row in consenso_df.iterrows():
            avisos = []
            crosstalk_aviso = row.get("ptm_crosstalk_aviso")
            if pd.notna(crosstalk_aviso):
                avisos.append(f"crosstalk:+{_crosstalk_count(crosstalk_aviso)}")
            if row.get("metoken_type_coincide") is False:
                metoken_tipo = row.get("metoken_type")
                avisos.append(f"MeToken<>{metoken_tipo}" if pd.notna(metoken_tipo) else "MeToken<>")
            cells = [
                int(row["posicion"]), row["residuo_wt"], row["tipo_ptm"], _abbrev_motor(row["motor"]),
                _fmt_score(row.get("score_deepmvp")), _fmt_score(row.get("score_deepptmpred")),
                ", ".join(avisos) if avisos else "-",
            ]
            if mostrar_quinasa:
                kinasa = row.get("kinase_library_top_kinase")
                cells.append(kinasa if pd.notna(kinasa) else "-")
            if ground_truth:
                entry = ground_truth.get((int(row["posicion"]), row["tipo_ptm"]))
                if entry:
                    tier, _fuente, pmids = entry
                    cells.extend([f"tier {tier}", _fmt_pmids(pmids)])
                else:
                    cells.extend(["-", "-"])
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

    print(f"\nReporte completo: {report_path}")


def run_single_input(input_path: Path, output_dir: Path) -> Path:
    """Runs the complete pipeline (FASTA or PDB path, according to ``route_input``) over a single
    input file. Extracted from ``main()`` so that batch mode (``main()`` with
    ``--input`` pointing at a directory) can reuse exactly the same per-file
    logic without duplicating it.

    It returns the path of the final report. Any ``PipelineError`` from the internal phases is
    propagated uncaught: the caller decides whether that stops the whole process (single-file
    mode) or only that file of the batch (batch mode, see ``main()``).
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
        _print_pdb_summary(
            record, deepmvp_results, deepptmpred_results, filtered, n_evaluados, report_path
        )

    print(f"\n{INTERPRETATION_DISCLAIMER}")
    return report_path


def _run_batch(input_dir: Path, output_dir: Path) -> int:
    """Batch mode: runs ``run_single_input`` over each recognized file of ``input_dir``.

    A file that fails is recorded as an error and does NOT stop the rest of the batch, the same
    non-fatal degradation criterion that MeToken applies per site (one
    real individual failure must not bring down the whole sweep). It writes
    ``batch_summary.csv`` (columns: file, status, report/error) in ``output_dir``.
    Exit code: 0 only if ALL the files completed without error, 1 if at least one
    failed (this includes the case of a directory without any recognized file; see the docstring of
    ``_discover_batch_inputs``, as "empty" is also a result that must not pass silently).

    It catches ``Exception`` in general, not only ``PipelineError`` (a real bug found in the
    2026-08-07 audit: previously, any unexpected exception, for example a real ``KeyError``/
    ``AttributeError`` from an engine with malformed output, not necessarily
    a ``PipelineError``, escaped the whole loop, bringing down the ENTIRE batch before
    writing ``batch_summary.csv`` and losing the record of the previously
    processed files, even though their individual reports WERE on disk).
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
