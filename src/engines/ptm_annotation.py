"""Phase 3 (core): annotation/filtering (B) + flow decision logic (D).

The design was closed in the architecture decision of 2026-07-27 (see
``01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md``
in the vault): the residue position as the primary unit, a window only for
types with a defined biological motif, a per-tool threshold with the raw
score always preserved, and types without corroboration included in the
core marked as such. D (``apply_workflow_filter``) is a single-responsibility
filter over ``pasa_umbral``: it did not route to Phase A/Extension 3
while that phase existed, and that decision was closed completely when Phase
A/3c was removed from the project's scope (project supervision feedback,
2026-08-10, see
the project decision record).

DeepMVP <-> DeepPTMPred type correspondence (verified by reading both
repos on 2026-07-27, not assumed): of the 6 biological types that
DeepMVP covers (8 residue-specific models), 7 have a direct equivalent in
DeepPTMPred; ``phosphorylation_y`` (phosphorylation on tyrosine) has NO
equivalent model in DeepPTMPred (its only phosphorylation model
covers S/T, not Y), which is why, even on the PDB path, DeepMVP's
``phosphorylation_y`` predictions never have a possible consensus. This is not a
scope exclusion but a real coverage gap between the two engines.
The 10 types exclusive to DeepPTMPred (hydroxylation,
gamma_carboxyglutamic_acid, malonylation, crotonylation, succinylation,
glutathionylation, s_nitrosylation, glutarylation, citrullination,
o_linked_glycosylation) likewise have no possible consensus for the
inverse reason (DeepMVP does not cover them).

## Optional corroboration of TYPE (MeToken, decision 2026-08-01)

``annotate_pdb_path`` optionally accepts ``pdb_path``/``chain_id`` (both
``None``/default if not passed, with behaviour identical to that before
this improvement): if they are provided and ``Settings.METOKEN_ENABLED`` is ``True``, then for
each row with ``pasa_umbral=True`` (sites that the consensus has ALREADY accepted)
MeToken is invoked (``src/engines/metoken_engine.py``, see its docstring and that of
``src/engines/_metoken_runner.py`` for the complete detail of what it is and
why) and 3 purely informative columns are added: ``metoken_type`` (the
type with the highest probability according to MeToken at that position, out of the 24
real classes, excluding the null class and the grouped "rare" class, see the runner),
``metoken_probability`` (its raw probability) and ``metoken_type_coincide``
(``True``/``False`` according to whether ``tipo_ptm`` has a known equivalent in MeToken and
matches or not, and ``None`` if there is no mapped equivalent, see
``CANONICAL_TO_METOKEN_TYPE``, or if MeToken could not evaluate that position).
MeToken NEVER changes ``pasa_umbral``/``consenso``, the same
non-deciding pattern as the secretory pathway
(``src/structural/uniprot_localization_client.py``): a failure (the repo not
installed, the subprocess failing, a timeout) leaves the 3 columns as ``None``
for the whole table, without affecting the rest of the report.

## Real N-GLYCOSYLATION consensus, PDB path (EMNGly, decision 2026-08-06)

``annotate_pdb_path`` promotes EMNGly (``src/engines/emngly_engine.py``,
the replacement for CoNglyPred, confirmed as having no published weights anywhere,
see the project decision record, "Decision 2") to a REAL CONSENSUS ENGINE alongside DeepMVP,
specifically for the rows with ``tipo_ptm`` in ``_NGLYCO_TYPES`` and
``motor == 'DeepMVP'`` (candidates proposed by DeepMVP that DeepPTMPred
never merges for this type, see ``CONSENSUS_EXCLUDED_TYPES`` below;
DeepPTMPred is still reported separately, unchanged). It is activated
automatically when ``pdb_path`` is not ``None`` and
``Settings.EMNGLY_ENABLED`` (the same activation pattern as MeToken).
StackGlyEmbed ran here as a third consensus engine from
2026-08-06 and was removed from the project on 2026-08-10 by project decision
(it depended on a sibling project's venv, real friction for the Scipion
integration); the consensus degrades naturally to 2 engines,
without code changes beyond removing its branch (see
``_apply_nglyco_consensus``).

The rule (provisional, see ``Settings.NGLYCO_CONSENSUS_MIN_ENGINES``): of the
engines that managed to evaluate the position (DeepMVP always available;
EMNGly degrades to absent without raising if it is not installed, see
``_validate_installation``), ``pasa_umbral`` = at least 1 passes its own
threshold (generalizing the 2-engine OR rule already used for the other
types), and ``consenso`` = at least ``Settings.NGLYCO_CONSENSUS_MIN_ENGINES``
(default 2) pass, which in practice requires BOTH engines to pass.
``motor`` is rewritten to reflect exactly which engines managed to
evaluate the position (``'DeepMVP+EMNGly'``, or ``'DeepMVP'`` if EMNGly is
degraded, in which case ``consenso`` is always ``False``, as it is impossible
to reach the minimum with a single engine).

## Optional corroboration of the SECRETORY PATHWAY, N-glycosylation, both paths
   (biological coherence analysis 2026-08-07)

Neither DeepMVP nor DeepPTMPred nor EMNGly models the substrate's real
biosynthetic pathway. They all predict from sequence/structure whether an
N-X-[S/T] sequon IS glycosylatable, never whether the protein actually transits the
ER/Golgi (a real chemical requirement of N-glycosylation). For
each row with ``tipo_ptm`` in ``_NGLYCO_TYPES`` and ``pasa_umbral=True`` (a
site that the consensus has ALREADY accepted), UniProt is queried
(``src/structural/uniprot_localization_client.py``) ONCE per accession and
the PURELY INFORMATIVE column ``via_secretora_evidencia`` is added:
``True`` (UniProt reports a localization consistent with the secretory pathway),
``False`` (UniProt DOES have real localization data but NONE of it is
consistent, so real evidence against, not merely absence), or ``None`` (no
localization data in UniProt, an unrecognized accession, which is the most
common case here since ``accession`` normally comes from the name of the input
file rather than a real UniProt ID, or a network failure). It NEVER changes
``pasa_umbral``/``consenso``, the same non-deciding pattern as MeToken.
It deliberately does NOT apply to
``o_linked_glycosylation``: there are two biologically distinct
O-glycosylation pathways (cytoplasmic/nuclear O-GlcNAc vs mucin-type
O-glycosylation in the secretory pathway) and this client does not distinguish them, so the
scope is limited to N-glycosylation so as not to risk an incorrect biological
claim.

## COMPETITION/CROSSTALK warning between PTMs of the same residue (biological
   coherence analysis 2026-08-07)

Several groups of this project's PTM types chemically modify the
SAME functional group of a residue and are mutually exclusive in the
same molecule at the same instant; they literally cannot coexist
on the same atom. The core scores each type/position
independently, so if the consensus accepts 2+ competing types
at exactly the same position, today they are presented as equally
valid findings without any signal that they are exclusive. The groups are
real and conservative (see ``_PTM_COMPETITION_GROUPS``, based on
established PTM crosstalk literature, deliberately not exhaustive, and each
group is additionally verified against the row's real ``residuo_wt``, never
assumed merely from the type):
- Lysine (lysine acyl): acetylation, ubiquitination, sumoylation,
  Lys methylation, malonylation, glutarylation, succinylation,
  crotonylation, all occupying the same epsilon-amino group.
- Cysteine (thiol): S-nitrosylation, glutathionylation.
- Arginine (guanidino): Arg methylation, citrullination (citrullination
  can additionally block prior methylation biologically, a crosstalk
  documented in the PAD/PRMT literature).
- Serine/Threonine (hydroxyl): phosphorylation, O-glycosylation (the
  phospho/O-GlcNAc "Yin-Yang" hypothesis, widely documented).

``_add_ptm_crosstalk_flag`` (see below) adds the PURELY
INFORMATIVE column ``ptm_crosstalk_aviso`` to ANY row with ``pasa_umbral=True``
that shares group+position+accession with at least one other row also
``pasa_umbral=True``. It NEVER changes ``pasa_umbral``/``consenso``. It
runs on both paths (FASTA and PDB), after all the other
corroborations, so as to see each row's final ``pasa_umbral`` state
(for example after the N-glycosylation consensus).

## Optional corroboration of KINASE SPECIFICITY, phosphorylation, both
   paths (biological coherence analysis 2026-08-07, point 5)

Neither DeepMVP nor DeepPTMPred distinguishes WHICH kinase family phosphorylates a
site; both predict "phosphorylatable in general". Unlike the
polyubiquitin chain type (a later cellular event, not derivable from the
structure, a line of reasoning that lived in Phase A/3c, removed from scope
2026-08-10), kinase specificity IS a local sequence property
around the site, and there is
a real, published and downloadable source (Johnson et al. 2023 Nature +
Yaron-Barir et al. 2024 Nature, packaged in ``kinase-library``, see the
docstring of ``src/engines/_kinase_library_runner.py``). For each row with
``tipo_ptm`` in ``_PHOSPHO_TYPES`` and ``pasa_umbral=True``, 4
PURELY INFORMATIVE columns are added: ``kinase_library_top_kinase``,
``kinase_library_top_family``, ``kinase_library_percentile`` and
``kinase_library_top3_kinases``. They NEVER change ``pasa_umbral``/``consenso``.
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from src.config.settings import Settings
from src.engines.emngly_engine import get_emngly_predictions
from src.engines.kinase_library_engine import get_kinase_corroboration
from src.engines.metoken_engine import get_type_corroboration
from src.structural.uniprot_localization_client import (
    UniProtLookupError,
    lookup_secretory_pathway_evidence,
)
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

# Tipo canonico del pipeline -> etiqueta real de MeToken
# (src/constant.py::PTMtype_list del repo clonado, verificado leyendo el
# archivo real, no adivinado). MeToken no distingue metilacion de Lys vs Arg
# (una sola clase "Methylation") -- ambos tipos del pipeline mapean ahi.
# crotonylation/glutarylation/citrullination NO tienen equivalente en
# MeToken (ausentes de sus 24 clases reales): deliberadamente fuera de este
# dict, 'metoken_type_coincide' queda 'None' (no evaluable) para esos 3.
CANONICAL_TO_METOKEN_TYPE = {
    "phosphorylation": "Phosphorylation",
    # 'phosphorylation_y' (DeepMVP crudo, sin equivalente en DeepPTMPred --
    # ver DEEPMVP_TO_CANONICAL_TYPE) SI tiene equivalente aqui: MeToken no
    # distingue residuo para fosforilacion (una sola clase "Phosphorylation"
    # para S/T/Y, confirmado en src/constant.py -- no hay "Phosphorylation_Y"),
    # mismo criterio que lys/arg_methylation -> "Methylation" abajo.
    "phosphorylation_y": "Phosphorylation",
    "acetylation": "Acetylation",
    "ubiquitination": "Ubiquitination",
    "hydroxylation": "Hydroxylation",
    "gamma_carboxyglutamic_acid": "Gamma-carboxyglutamic acid",
    "lys_methylation": "Methylation",
    "arg_methylation": "Methylation",
    "malonylation": "Malonylation",
    "succinylation": "Succinylation",
    "glutathionylation": "Glutathionylation",
    "sumoylation": "Sumoylation",
    "s_nitrosylation": "S-nitrosylation",
    "o_linked_glycosylation": "O-linked Glycosylation",
    "n_linked_glycosylation": "N-linked Glycosylation",
}

# Nombre canonico = el que usa DeepPTMPred (superset de tipos). Mapea el
# nombre especifico de modelo de DeepMVP (residuo incluido en el nombre) a
# su equivalente. 'phosphorylation_y' se deja fuera deliberadamente: no
# tiene equivalente, se reporta con su propio nombre (ver docstring modulo).
DEEPMVP_TO_CANONICAL_TYPE = {
    "acetylation_k": "acetylation",
    "glycosylation_n": "n_linked_glycosylation",
    "methylation_k": "lys_methylation",
    "methylation_r": "arg_methylation",
    "phosphorylation_st": "phosphorylation",
    "sumoylation_k": "sumoylation",
    "ubiquitination_k": "ubiquitination",
}

# Tipos con motivo de secuencia biologico definido -> funcion que calcula la
# ventana (o None si la posicion no cumple el motivo). Unico caso hoy:
# sequon de N-glicosilacion N-X-[S/T] (X != P por convencion, aplicada
# aqui). Cualquier tipo no listado aqui no tiene ventana (motivo puntual:
# fosforilacion/acetilacion/metilacion/ubiquitinacion/etc.).
_NGLYCO_TYPES = {"n_linked_glycosylation", "glycosylation_n"}

# tipo_ptm (crudo de DeepMVP o canonico de DeepPTMPred) de fosforilacion --
# usado por ``_add_kinase_library_corroboration`` (analisis de coherencia
# biologica 2026-08-07, punto 5). Incluye 'phosphorylation_y' (Tirosina, sin
# equivalente en DeepPTMPred, ver DEEPMVP_TO_CANONICAL_TYPE) porque Kinase
# Library cubre el kinoma Tyr completo ademas del Ser/Thr (Yaron-Barir et al.
# 2024 Nature) -- a diferencia de ``_PTM_COMPETITION_GROUPS`` mas abajo, que
# deliberadamente la excluye por un motivo distinto (grupo quimico sin
# competencia real conocida).
_PHOSPHO_TYPES = {"phosphorylation", "phosphorylation_st", "phosphorylation_y"}

# tipo_ptm (nombre CRUDO de DeepMVP o CANONICO de DeepPTMPred, ambos
# posibles en esta columna, ver DEEPMVP_TO_CANONICAL_TYPE) -> (grupo de
# competencia quimica real, residuos que ese grupo puede ocupar). Usado por
# ``_add_ptm_crosstalk_flag`` -- ver docstring del modulo para la
# justificacion biologica de cada grupo. Deliberadamente conservador: solo
# incluye tipos de este proyecto con competencia real bien establecida en la
# literatura por el MISMO grupo funcional de un residuo. Cada grupo se
# valida ademas contra ``residuo_wt`` real de la fila (nunca se asume el
# residuo objetivo solo por el nombre del tipo).
_PTM_COMPETITION_GROUPS = {
    # Lisina (acilo-lisina): grupo epsilon-amino, unico entre todos estos.
    "acetylation": ("K_acilo", {"K"}),
    "acetylation_k": ("K_acilo", {"K"}),
    "ubiquitination": ("K_acilo", {"K"}),
    "ubiquitination_k": ("K_acilo", {"K"}),
    "sumoylation": ("K_acilo", {"K"}),
    "sumoylation_k": ("K_acilo", {"K"}),
    "lys_methylation": ("K_acilo", {"K"}),
    "methylation_k": ("K_acilo", {"K"}),
    "malonylation": ("K_acilo", {"K"}),
    "glutarylation": ("K_acilo", {"K"}),
    "succinylation": ("K_acilo", {"K"}),
    "crotonylation": ("K_acilo", {"K"}),
    # Cisteina (tiol).
    "s_nitrosylation": ("C_tiol", {"C"}),
    "glutathionylation": ("C_tiol", {"C"}),
    # Arginina (guanidino).
    "arg_methylation": ("R_guanidino", {"R"}),
    "methylation_r": ("R_guanidino", {"R"}),
    "citrullination": ("R_guanidino", {"R"}),
    # Serina/Treonina (hidroxilo) -- NUNCA incluye 'phosphorylation_y'
    # (fosforilacion en Tirosina, sin equivalente en DeepPTMPred, ver
    # DEEPMVP_TO_CANONICAL_TYPE): ocupa un grupo quimico distinto (hidroxilo
    # fenolico, no alcoholico), sin competencia real conocida con los tipos
    # de este proyecto.
    "phosphorylation": ("ST_hidroxilo", {"S", "T"}),
    "phosphorylation_st": ("ST_hidroxilo", {"S", "T"}),
    "o_linked_glycosylation": ("ST_hidroxilo", {"S", "T"}),
}

# Tipos excluidos deliberadamente de la fusion DeepMVP+DeepPTMPred (decision
# 2026-08-01, ver el registro de decisiones del proyecto, "investigacion de n_linked_glycosylation y los 4
# tipos mediocres"): DeepPTMPred no tiene poder discriminativo real para
# estos tipos (verificado contra las metricas de entrenamiento del propio
# repo, AUC 0.495 -- no es un problema de umbral). Ambos motores se siguen
# reportando, nunca fusionados en una fila de consenso entre ELLOS DOS. Esto
# NO significa que el tipo se quede sin consenso posible: desde 2026-08-06,
# 'n_linked_glycosylation'/'glycosylation_n' en Camino PDB tiene un consenso
# real DISTINTO (DeepMVP+EMNGly, ver ``_apply_nglyco_consensus`` mas abajo)
# que reemplaza el rol que hubiera tenido DeepPTMPred si no estuviera muerto
# para este tipo.
CONSENSUS_EXCLUDED_TYPES = {"n_linked_glycosylation"}

OUTPUT_COLUMNS = [
    "accession", "posicion", "residuo_wt", "tipo_ptm", "motor",
    "score_deepmvp", "score_deepptmpred", "consenso", "ventana", "camino",
    "pasa_umbral",
]


def _nglyco_window(sequence: str, position: int) -> Optional[str]:
    """Sequon N-X-[S/T] (X != P) centrado en ``position`` (1-based). ``None`` si no aplica."""
    idx = position - 1
    if idx < 0 or idx + 2 >= len(sequence):
        return None
    motif = sequence[idx : idx + 3]
    if len(motif) == 3 and motif[0] == "N" and motif[1] != "P" and motif[2] in ("S", "T"):
        return motif
    return None


def _window_for(tipo_ptm: str, sequence: str, position: int) -> Optional[str]:
    if tipo_ptm in _NGLYCO_TYPES:
        return _nglyco_window(sequence, position)
    return None


def annotate_fasta_path(accession: str, sequence: str, deepmvp_df: pd.DataFrame) -> pd.DataFrame:
    """Fase 3 nucleo (B), Camino FASTA: solo DeepMVP, sin consenso posible.

    Args:
        accession: Accession de la proteina (misma que ``deepmvp_df['protein']``).
        sequence: Secuencia saneada de Fase 1 (para calcular ventanas).
        deepmvp_df: Salida cruda de ``DeepMVPEngine.run()`` para este accession
            (columnas ``protein|aa|pos|x|y_pred|fpr|ptm``).

    Returns:
        DataFrame con ``OUTPUT_COLUMNS``, una fila por sitio candidato
        reportado por DeepMVP, sin ningun filtrado (ver
        :func:`apply_workflow_filter` para la capa D).
    """
    rows = []
    for _, r in deepmvp_df.iterrows():
        rows.append({
            "accession": accession,
            "posicion": int(r["pos"]),
            "residuo_wt": r["aa"],
            "tipo_ptm": r["ptm"],
            "motor": "DeepMVP",
            "score_deepmvp": float(r["y_pred"]),
            "score_deepptmpred": None,
            "consenso": False,
            "ventana": _window_for(r["ptm"], sequence, int(r["pos"])),
            "camino": "FASTA",
            "pasa_umbral": bool(r["fpr"] <= Settings.DEEPMVP_MAX_FPR),
        })
    result = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    if Settings.SECRETORY_PATHWAY_CHECK_ENABLED:
        try:
            result = _add_secretory_pathway_evidence(result, accession)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, el cliente ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo evidencia de via secretora para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
            )

    if Settings.KINASE_LIBRARY_ENABLED:
        try:
            result = _add_kinase_library_corroboration(result, sequence)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, kinase_library_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion Kinase Library para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
            )

    if Settings.PTM_CROSSTALK_CHECK_ENABLED:
        result = _add_ptm_crosstalk_flag(result)

    return result


def annotate_pdb_path(
    accession: str, sequence: str, deepmvp_df: pd.DataFrame, deepptmpred_df: pd.DataFrame,
    pdb_path: Optional[Path] = None, chain_id: str = "A",
    position_mapping: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Fase 3 nucleo (B), Camino PDB: consenso DeepMVP + DeepPTMPred donde exista.

    Args:
        accession: Accession de la proteina.
        sequence: Secuencia ATMSEQ de Fase 1.5 (para calcular ventanas).
        deepmvp_df: Salida cruda de ``DeepMVPEngine.run()`` para este accession.
        deepptmpred_df: Salida cruda de ``DeepPTMPredEngine.run()`` para este
            accession (columnas ``protein_id|position|residue|probability|ptm_type``).
        pdb_path: Opcional -- PDB de una sola cadena (``record.chain_pdb_path``
            de Fase 1.5, para que las posiciones 1-based coincidan con
            ``sequence``) para habilitar la corroboracion informativa de
            tipo via MeToken Y el consenso real de N-glicosilacion via
            EMNGly (ver docstring del modulo). ``None`` (default) desactiva
            ambas sin cambiar ningun otro comportamiento -- identico a antes
            de esta mejora.
        chain_id: Cadena a leer del PDB si ``pdb_path`` no es ``None``
            (default ``"A"``).
        position_mapping: Opcional -- tabla de Fase 1.5
            (``StructureRecord.position_mapping``, columnas
            ``fasta_position``/``pdb_seqid``) necesaria para que EMNGly
            alinee ``structure_emb`` con la numeracion real del PDB (ver
            docstring de ``_emngly_runner.py``). ``None`` (default) hace que
            EMNGly se omita del consenso de N-glicosilacion aunque
            ``pdb_path`` este presente (degradacion no fatal, igual que
            EMNGly no instalado).

    Returns:
        DataFrame con ``OUTPUT_COLUMNS``. Sitios reportados por ambos
        motores para el mismo tipo canonico y posicion se fusionan en una
        fila (``motor='DeepMVP+DeepPTMPred'``, ``consenso`` refleja si
        ambos pasan su propio umbral). Sitios reportados por un solo motor
        (incluye ``phosphorylation_y`` de DeepMVP y los 10 tipos exclusivos
        de DeepPTMPred, ver docstring del modulo) se incluyen igual,
        marcados ``consenso=False``. Si ``pdb_path`` no es ``None`` y
        ``Settings.METOKEN_ENABLED``, incluye ademas ``metoken_type``/
        ``metoken_probability``/``metoken_type_coincide`` (ver docstring del
        modulo, independiente de todo lo demas). Para filas nglyco
        (``tipo_ptm`` en ``_NGLYCO_TYPES``), si el consenso real esta activo
        (``pdb_path`` no ``None`` y ``Settings.EMNGLY_ENABLED``, ver
        ``_apply_nglyco_consensus``): ``motor``/``pasa_umbral``/``consenso``
        reflejan la fusion DeepMVP+EMNGly, y se agrega ``score_emngly``.
        Columna ausente (mismas ``OUTPUT_COLUMNS`` de siempre) si el
        consenso no esta activo.
    """
    ptmpred_lookup = {
        (int(r["position"]), r["ptm_type"]): r
        for _, r in deepptmpred_df.iterrows()
    }
    matched_keys = set()
    rows = []

    for _, r in deepmvp_df.iterrows():
        pos = int(r["pos"])
        tipo_deepmvp = r["ptm"]
        tipo_canonico = DEEPMVP_TO_CANONICAL_TYPE.get(tipo_deepmvp, tipo_deepmvp)
        key = (pos, tipo_canonico)
        ptmpred_row = (
            ptmpred_lookup.get(key) if tipo_canonico not in CONSENSUS_EXCLUDED_TYPES else None
        )

        deepmvp_pasa = bool(r["fpr"] <= Settings.DEEPMVP_MAX_FPR)

        if ptmpred_row is not None:
            matched_keys.add(key)
            deepptmpred_pasa = bool(
                ptmpred_row["probability"] >= Settings.deepptmpred_threshold_for(tipo_canonico)
            )
            rows.append({
                "accession": accession,
                "posicion": pos,
                "residuo_wt": r["aa"],
                "tipo_ptm": tipo_canonico,
                "motor": "DeepMVP+DeepPTMPred",
                "score_deepmvp": float(r["y_pred"]),
                "score_deepptmpred": float(ptmpred_row["probability"]),
                "consenso": bool(deepmvp_pasa and deepptmpred_pasa),
                "ventana": _window_for(tipo_canonico, sequence, pos),
                "camino": "PDB",
                "pasa_umbral": bool(deepmvp_pasa or deepptmpred_pasa),
            })
        else:
            rows.append({
                "accession": accession,
                "posicion": pos,
                "residuo_wt": r["aa"],
                "tipo_ptm": tipo_deepmvp,
                "motor": "DeepMVP",
                "score_deepmvp": float(r["y_pred"]),
                "score_deepptmpred": None,
                "consenso": False,
                "ventana": _window_for(tipo_deepmvp, sequence, pos),
                "camino": "PDB",
                "pasa_umbral": deepmvp_pasa,
            })

    for (pos, tipo_canonico), r in ptmpred_lookup.items():
        if (pos, tipo_canonico) in matched_keys:
            continue
        rows.append({
            "accession": accession,
            "posicion": pos,
            "residuo_wt": r["residue"],
            "tipo_ptm": tipo_canonico,
            "motor": "DeepPTMPred",
            "score_deepmvp": None,
            "score_deepptmpred": float(r["probability"]),
            "consenso": False,
            "ventana": _window_for(tipo_canonico, sequence, pos),
            "camino": "PDB",
            "pasa_umbral": bool(r["probability"] >= Settings.deepptmpred_threshold_for(tipo_canonico)),
        })

    result = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    # Orden real corregido en auditoria 2026-08-07: el consenso de N-glicosilacion
    # (_apply_nglyco_consensus) DEBE correr antes que MeToken -- MeToken calcula sus
    # filas elegibles a partir de 'pasa_umbral' EN EL MOMENTO en que corre (ver
    # docstring de _add_metoken_corroboration), y el consenso de nglyco puede PROMOVER
    # una fila de pasa_umbral=False a True (DeepMVP solo no paso su propio umbral, pero
    # EMNGly si). Con el orden viejo (MeToken primero), esas filas
    # promovidas nunca recibian corroboracion MeToken -- quedaban en None
    # permanentemente pese a que la fila ya paso el consenso, contradiciendo la propia
    # documentacion de MeToken de ser "independiente de todo lo demas". Mismo criterio
    # ya aplicado a '_add_ptm_crosstalk_flag' (corre ultimo, tras ver el pasa_umbral
    # final de cada fila).
    if pdb_path is not None and Settings.EMNGLY_ENABLED:
        try:
            result = _apply_nglyco_consensus(result, accession, sequence, pdb_path, position_mapping)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, emngly_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado aplicando el consenso de N-glicosilacion (EMNGly) "
                "para '%s' (no fatal, filas nglyco quedan como motor='DeepMVP' sin cambios): %s",
                accession, exc,
            )

    if pdb_path is not None and Settings.METOKEN_ENABLED:
        try:
            result = _add_metoken_corroboration(result, pdb_path, chain_id)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, metoken_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion MeToken para '%s' (no fatal, no afecta "
                "consenso/pasa_umbral): %s", accession, exc,
            )

    if Settings.SECRETORY_PATHWAY_CHECK_ENABLED:
        try:
            result = _add_secretory_pathway_evidence(result, accession)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, el cliente ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo evidencia de via secretora para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
            )

    if Settings.KINASE_LIBRARY_ENABLED:
        try:
            result = _add_kinase_library_corroboration(result, sequence)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, kinase_library_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion Kinase Library para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
            )

    if Settings.PTM_CROSSTALK_CHECK_ENABLED:
        result = _add_ptm_crosstalk_flag(result)

    return result


def _add_metoken_corroboration(result: pd.DataFrame, pdb_path: Path, chain_id: str) -> pd.DataFrame:
    """Anade ``metoken_type``/``metoken_probability``/``metoken_type_coincide`` a las filas elegibles.

    Puramente informativo (ver docstring del modulo): NUNCA modifica
    ``pasa_umbral``/``consenso`` de ninguna fila, solo rellena estas 3
    columnas nuevas para las filas con ``pasa_umbral=True`` (sitios que el
    consenso YA acepto, decision 2026-08-01). Si MeToken no devuelve nada
    (no instalado, subproceso fallido, timeout -- ver
    ``metoken_engine.get_type_corroboration``, degrada sin lanzar), las 3
    columnas quedan en ``None`` para toda la tabla.

    ``tipo_ptm`` puede venir como nombre CANONICO (fila con consenso, o
    exclusiva de DeepPTMPred) o como nombre CRUDO especifico de DeepMVP
    (``acetylation_k``, ``methylation_r``, etc. -- filas DeepMVP-solo sin
    match de DeepPTMPred, ver ``annotate_pdb_path``) -- se normaliza via
    ``DEEPMVP_TO_CANONICAL_TYPE`` antes de buscarlo en
    ``CANONICAL_TO_METOKEN_TYPE``, si no fuera asi ``metoken_type_coincide``
    quedaria ``None`` incorrectamente para toda fila DeepMVP-solo (hallazgo
    real durante el testing de este wiring, no solo teorico).
    """
    result = result.copy()
    result["metoken_type"] = None
    result["metoken_probability"] = None
    result["metoken_type_coincide"] = None

    eligible_positions = sorted(set(result.loc[result["pasa_umbral"], "posicion"]))
    if not eligible_positions:
        return result

    corroboration = get_type_corroboration(pdb_path, eligible_positions, chain_id=chain_id)
    if not corroboration:
        return result

    for idx, row in result.iterrows():
        if not row["pasa_umbral"]:
            continue
        site = corroboration.get(int(row["posicion"]))
        if site is None:
            continue
        result.at[idx, "metoken_type"] = site["metoken_type"]
        result.at[idx, "metoken_probability"] = site["metoken_probability"]
        tipo_canonico = DEEPMVP_TO_CANONICAL_TYPE.get(row["tipo_ptm"], row["tipo_ptm"])
        expected = CANONICAL_TO_METOKEN_TYPE.get(tipo_canonico)
        if expected is not None:
            result.at[idx, "metoken_type_coincide"] = bool(expected == site["metoken_type"])

    return result


def _add_secretory_pathway_evidence(result: pd.DataFrame, accession: str) -> pd.DataFrame:
    """Anade ``via_secretora_evidencia`` a las filas de N-glicosilacion elegibles.

    Puramente informativo (ver docstring del modulo): NUNCA modifica
    ``pasa_umbral``/``consenso``. Una unica consulta a UniProt por accession
    (nunca por fila/posicion -- la localizacion subcelular es una propiedad
    de la proteina completa, no de un sitio individual), reutilizada para
    todas las filas elegibles (``tipo_ptm`` en ``_NGLYCO_TYPES`` -- excluye
    deliberadamente ``o_linked_glycosylation``, ver docstring del modulo --
    con ``pasa_umbral=True``).

    Un accession no reconocido por UniProt (el caso mas comun aqui, ver
    docstring del modulo) o un fallo de red dejan la columna en ``None`` para
    toda la tabla -- nunca ``False`` por defecto, para no afirmar "sin via
    secretora" cuando en realidad es "no se pudo verificar".
    """
    result = result.copy()
    result["via_secretora_evidencia"] = None

    eligible_mask = result["tipo_ptm"].isin(_NGLYCO_TYPES) & result["pasa_umbral"]
    if not eligible_mask.any():
        return result

    try:
        evidencia = lookup_secretory_pathway_evidence(accession)
    except UniProtLookupError as exc:
        logger.warning(
            "No se pudo consultar UniProt para '%s' -- 'via_secretora_evidencia' queda sin "
            "determinar para todas las filas de N-glicosilacion (no fatal): %s", accession, exc,
        )
        return result

    result.loc[eligible_mask, "via_secretora_evidencia"] = evidencia
    return result


def _add_kinase_library_corroboration(result: pd.DataFrame, sequence: str) -> pd.DataFrame:
    """Anade ``kinase_library_top_kinase``/``kinase_library_top_family``/
    ``kinase_library_percentile``/``kinase_library_top3_kinases`` a las filas de fosforilacion elegibles.

    Analisis de coherencia biologica 2026-08-07, punto 5 (ver docstring de
    ``src/engines/_kinase_library_runner.py`` para la fuente real usada --
    Johnson et al. 2023 Nature + Yaron-Barir et al. 2024 Nature). Puramente
    informativo: NUNCA modifica ``pasa_umbral``/``consenso``, solo rellena
    estas 4 columnas nuevas para las filas de fosforilacion (``tipo_ptm`` en
    ``_PHOSPHO_TYPES``) con ``pasa_umbral=True``. Si Kinase Library no
    devuelve nada (entorno no instalado, subproceso fallido, timeout -- ver
    ``kinase_library_engine.get_kinase_corroboration``, degrada sin lanzar),
    las 4 columnas quedan en ``None`` para toda la tabla.

    A diferencia de ``metoken_type_coincide``, no hay ninguna columna
    "coincide": Kinase Library no predice SI el sitio se fosforila (eso ya
    lo decidio el consenso), solo POR QUIEN -- no hay una prediccion previa
    con la que comparar.
    """
    result = result.copy()
    result["kinase_library_top_kinase"] = None
    result["kinase_library_top_family"] = None
    result["kinase_library_percentile"] = None
    result["kinase_library_top3_kinases"] = None

    eligible_mask = result["tipo_ptm"].isin(_PHOSPHO_TYPES) & result["pasa_umbral"]
    eligible_positions = sorted(set(result.loc[eligible_mask, "posicion"]))
    if not eligible_positions:
        return result

    corroboration = get_kinase_corroboration(sequence, eligible_positions)
    if not corroboration:
        return result

    for idx, row in result.loc[eligible_mask].iterrows():
        site = corroboration.get(int(row["posicion"]))
        if site is None:
            continue
        result.at[idx, "kinase_library_top_kinase"] = site["kinase_library_top_kinase"]
        result.at[idx, "kinase_library_top_family"] = site["kinase_library_top_family"]
        result.at[idx, "kinase_library_percentile"] = site["kinase_library_percentile"]
        result.at[idx, "kinase_library_top3_kinases"] = site["kinase_library_top3_kinases"]

    return result


def _add_ptm_crosstalk_flag(result: pd.DataFrame) -> pd.DataFrame:
    """Anade ``ptm_crosstalk_aviso`` cuando 2+ tipos en competencia real coinciden en el mismo residuo.

    Puramente informativo (ver docstring del modulo): NUNCA modifica
    ``pasa_umbral``/``consenso``. Solo considera filas ``pasa_umbral=True``
    (sitios que el consenso YA acepto) -- agrupadas por
    ``(accession, posicion, grupo de competencia)`` via
    ``_PTM_COMPETITION_GROUPS``, validando ademas que ``residuo_wt`` de cada
    fila coincida con el residuo real que ese grupo puede ocupar (nunca se
    asume el residuo objetivo solo por el nombre del tipo). Si 2+ tipos
    DISTINTOS del mismo grupo pasan en la misma posicion, cada fila
    involucrada recibe un aviso listando los otros tipos en competencia.
    """
    result = result.copy()
    result["ptm_crosstalk_aviso"] = None

    eligible = result[result["pasa_umbral"]]
    if eligible.empty:
        return result

    groups: dict = {}
    for idx, row in eligible.iterrows():
        group_info = _PTM_COMPETITION_GROUPS.get(row["tipo_ptm"])
        if group_info is None:
            continue
        group_name, expected_residues = group_info
        if row["residuo_wt"] not in expected_residues:
            continue
        key = (row["accession"], int(row["posicion"]), group_name)
        groups.setdefault(key, []).append((idx, row["tipo_ptm"]))

    for members in groups.values():
        distinct_types = sorted({DEEPMVP_TO_CANONICAL_TYPE.get(t, t) for _, t in members})
        if len(distinct_types) < 2:
            continue
        for idx, tipo in members:
            tipo_normalizado = DEEPMVP_TO_CANONICAL_TYPE.get(tipo, tipo)
            others = [t for t in distinct_types if t != tipo_normalizado]
            result.at[idx, "ptm_crosstalk_aviso"] = (
                f"Compite con: {', '.join(others)} (mismo residuo, mutuamente excluyentes en una "
                "misma molecula/instante -- ver docstring del modulo)"
            )

    return result


def _apply_nglyco_consensus(
    result: pd.DataFrame, accession: str, sequence: str, pdb_path: Path,
    position_mapping: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Consenso real de 'n_linked_glycosylation'/'glycosylation_n' (Camino PDB, decision 2026-08-06).

    A diferencia de ``_add_metoken_corroboration`` (puramente informativa,
    evalua solo filas YA con ``pasa_umbral=True``), esta funcion DECIDE
    ``pasa_umbral``/``consenso``/``motor`` -- evalua TODAS las filas nglyco
    propuestas por DeepMVP (``motor == 'DeepMVP'``, el unico motor de origen
    posible para este tipo, ver ``CONSENSUS_EXCLUDED_TYPES``), no solo las
    que ya pasaban el umbral de DeepMVP en solitario. StackGlyEmbed corria
    aqui como tercer motor de consenso hasta 2026-08-10 (eliminado del
    proyecto, ver el registro de decisiones del proyecto) -- el consenso queda en DeepMVP+EMNGly.

    No-op (devuelve ``result`` sin tocar, sin agregar columnas) si no hay
    ninguna fila nglyco elegible -- evita ensuciar el esquema de reportes que
    no tienen ningun candidato de N-glicosilacion.
    """
    nglyco_mask = result["tipo_ptm"].isin(_NGLYCO_TYPES) & (result["motor"] == "DeepMVP")
    if not nglyco_mask.any():
        return result

    result = result.copy()
    result["score_emngly"] = None

    positions = sorted(set(result.loc[nglyco_mask, "posicion"].astype(int)))

    emngly_results = {}
    if position_mapping is not None and not position_mapping.empty:
        mapping_dir = Settings.EMNGLY_CACHE_DIR
        mapping_dir.mkdir(parents=True, exist_ok=True)
        mapping_csv_path = mapping_dir / f"{accession}_position_mapping.csv"
        position_mapping.to_csv(mapping_csv_path, index=False)
        emngly_results = get_emngly_predictions(
            accession, sequence, positions, pdb_path, mapping_csv_path,
        )
    else:
        logger.warning(
            "'position_mapping' no disponible para '%s' -- EMNGly se omite del consenso de "
            "N-glicosilacion (no fatal, degrada a solo DeepMVP).", accession,
        )

    for idx, row in result.loc[nglyco_mask].iterrows():
        pos = int(row["posicion"])
        engines_ran = ["DeepMVP"]
        n_pass = int(bool(row["pasa_umbral"]))  # umbral de DeepMVP ya calculado (fpr)

        emngly_site = emngly_results.get(pos)
        if emngly_site is not None:
            engines_ran.append("EMNGly")
            score = emngly_site["emngly_probability"]
            result.at[idx, "score_emngly"] = score
            if score >= Settings.EMNGLY_MIN_PROBABILITY:
                n_pass += 1

        result.at[idx, "motor"] = "+".join(engines_ran)
        result.at[idx, "pasa_umbral"] = bool(n_pass >= 1)
        result.at[idx, "consenso"] = bool(n_pass >= Settings.NGLYCO_CONSENSUS_MIN_ENGINES)
        # Bug real encontrado en auditoria 2026-08-07: esta funcion nunca canonizaba
        # tipo_ptm de 'glycosylation_n' (nombre crudo de DeepMVP, unico origen posible
        # aqui, ver nglyco_mask arriba) a 'n_linked_glycosylation'. _NGLYCO_TYPES/
        # _PTM_COMPETITION_GROUPS ya trataban ambos nombres como equivalentes, asi que
        # renombrar aqui no rompe ninguna otra verificacion de elegibilidad (el
        # consumidor original del nombre canonico, la seleccion de candidatos de Fase
        # A/3c, ya no existe -- eliminada del alcance 2026-08-10).
        result.at[idx, "tipo_ptm"] = "n_linked_glycosylation"

    return result


def apply_workflow_filter(annotated_df: pd.DataFrame) -> pd.DataFrame:
    """Fase 3 nucleo (D): filtro de responsabilidad unica sobre la tabla de B.

    Regla generica (decision 2026-07-27): mantiene solo las filas donde
    ``pasa_umbral`` es ``True`` (ya calculado en B usando el umbral propio de
    cada motor -- fpr calibrado para DeepMVP, probabilidad para DeepPTMPred).

    Args:
        annotated_df: Salida de :func:`annotate_fasta_path` o
            :func:`annotate_pdb_path`.

    Returns:
        Subconjunto de ``annotated_df`` que pasa el umbral, mismas columnas,
        orden preservado.
    """
    return annotated_df[annotated_df["pasa_umbral"]].reset_index(drop=True)
