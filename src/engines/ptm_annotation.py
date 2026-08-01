"""Fase 3 (nucleo): anotacion/filtrado (B) + logica de decision de flujo (D).

Diseno cerrado en la decision de arquitectura 2026-07-27 (ver
``01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md``
en el vault): posicion de residuo como unidad primaria, ventana solo para
tipos con motivo biologico definido, umbral propio por herramienta con
score crudo siempre conservado, tipos sin corroboracion incluidos en el
nucleo marcados como tal, D como filtro de responsabilidad unica sin rutas
a fases todavia inexistentes (Extension 3 / Fase A).

Correspondencia de tipos DeepMVP <-> DeepPTMPred (verificada leyendo ambos
repos el 2026-07-27, no asumida): de los 6 tipos biologicos que cubre
DeepMVP (8 modelos especificos de residuo), 7 tienen equivalente directo en
DeepPTMPred; ``phosphorylation_y`` (fosforilacion en tirosina) NO tiene
ningun modelo equivalente en DeepPTMPred (su unico modelo de fosforilacion
cubre S/T, no Y) -- por eso, incluso en Camino PDB, las predicciones de
``phosphorylation_y`` de DeepMVP nunca tienen consenso posible, no es un
descarte de alcance sino un hueco real de cobertura entre ambos motores.
Los 10 tipos exclusivos de DeepPTMPred (hydroxylation,
gamma_carboxyglutamic_acid, malonylation, crotonylation, succinylation,
glutathionylation, s_nitrosylation, glutarylation, citrullination,
o_linked_glycosylation) tampoco tienen consenso posible por el motivo
inverso (DeepMVP no los cubre).

## Corroboracion opcional de TIPO (MeToken, decision 2026-08-01)

``annotate_pdb_path`` acepta opcionalmente ``pdb_path``/``chain_id`` (ambos
``None``/default si no se pasan -- comportamiento identico al de antes de
esta mejora): si se proveen y ``Settings.METOKEN_ENABLED`` es ``True``, para
cada fila con ``pasa_umbral=True`` (sitios que el consenso YA acepto) se
invoca MeToken (``src/engines/metoken_engine.py``, ver su docstring y el de
``src/engines/_metoken_runner.py`` para el detalle completo de que es y por
que) y se anaden 3 columnas puramente informativas: ``metoken_type`` (el
tipo con mayor probabilidad segun MeToken en esa posicion, de las 24 clases
reales -- excluye la clase null y la clase "rare" agrupada, ver el runner),
``metoken_probability`` (su probabilidad cruda) y ``metoken_type_coincide``
(``True``/``False`` si ``tipo_ptm`` tiene equivalente conocido en MeToken y
coincide o no, ``None`` si no hay equivalente mapeado -- ver
``CANONICAL_TO_METOKEN_TYPE`` -- o si MeToken no pudo evaluar esa posicion).
MeToken NUNCA cambia ``pasa_umbral``/``consenso`` -- mismo patron
no-decisorio que GlyGen (``src/structural/glygen_client.py``): un fallo
(repo no instalado, subproceso revienta, timeout) deja las 3 columnas en
``None`` para toda la tabla, sin afectar el resto del reporte.
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from src.config.settings import Settings
from src.engines.metoken_engine import get_type_corroboration
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

# Tipos excluidos deliberadamente del consenso (decision 2026-08-01, ver
# STATUS.md "investigacion de n_linked_glycosylation y los 4 tipos
# mediocres"): DeepPTMPred no tiene poder discriminativo real para estos
# tipos (verificado contra las metricas de entrenamiento del propio repo,
# AUC 0.495 -- no es un problema de umbral). Ambos motores se siguen
# reportando, nunca fusionados en una fila de consenso.
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
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def annotate_pdb_path(
    accession: str, sequence: str, deepmvp_df: pd.DataFrame, deepptmpred_df: pd.DataFrame,
    pdb_path: Optional[Path] = None, chain_id: str = "A",
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
            tipo via MeToken (ver docstring del modulo). ``None`` (default)
            desactiva la corroboracion sin cambiar ningun otro
            comportamiento -- identico a antes de esta mejora.
        chain_id: Cadena a leer del PDB si ``pdb_path`` no es ``None``
            (default ``"A"``).

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
        modulo) -- ausentes (mismas ``OUTPUT_COLUMNS`` de siempre) en caso
        contrario.
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

    if pdb_path is not None and Settings.METOKEN_ENABLED:
        try:
            result = _add_metoken_corroboration(result, pdb_path, chain_id)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, metoken_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion MeToken para '%s' (no fatal, no afecta "
                "consenso/pasa_umbral): %s", accession, exc,
            )

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


def apply_workflow_filter(annotated_df: pd.DataFrame) -> pd.DataFrame:
    """Fase 3 nucleo (D): filtro de responsabilidad unica sobre la tabla de B.

    Regla generica (decision 2026-07-27): mantiene solo las filas donde
    ``pasa_umbral`` es ``True`` (ya calculado en B usando el umbral propio de
    cada motor -- fpr calibrado para DeepMVP, probabilidad para DeepPTMPred).
    Deliberadamente NO rutea a Extension 3 (DeltaDeltaG) ni a Fase A
    (modelado estructural): esas fases no existen todavia (ver decision).

    Args:
        annotated_df: Salida de :func:`annotate_fasta_path` o
            :func:`annotate_pdb_path`.

    Returns:
        Subconjunto de ``annotated_df`` que pasa el umbral, mismas columnas,
        orden preservado.
    """
    return annotated_df[annotated_df["pasa_umbral"]].reset_index(drop=True)
