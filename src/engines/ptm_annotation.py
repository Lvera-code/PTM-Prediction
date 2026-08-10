"""Fase 3 (nucleo): anotacion/filtrado (B) + logica de decision de flujo (D).

Diseno cerrado en la decision de arquitectura 2026-07-27 (ver
``01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md``
en el vault): posicion de residuo como unidad primaria, ventana solo para
tipos con motivo biologico definido, umbral propio por herramienta con
score crudo siempre conservado, tipos sin corroboracion incluidos en el
nucleo marcados como tal. D (``apply_workflow_filter``) es un filtro de
responsabilidad unica sobre ``pasa_umbral`` -- no ruteaba a Fase A/Extension 3
mientras existio esa fase, y esa decision quedo cerrada del todo cuando Fase
A/3c se elimino del alcance del proyecto (feedback de Carlos, 2026-08-10, ver
STATUS.md).

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
no-decisorio que la via secretora
(``src/structural/uniprot_localization_client.py``): un fallo (repo no
instalado, subproceso revienta, timeout) deja las 3 columnas en ``None``
para toda la tabla, sin afectar el resto del reporte.

## Corroboracion opcional de N-GLICOSILACION -- Camino FASTA (StackGlyEmbed, decision 2026-08-01)

``annotate_fasta_path`` acepta opcionalmente ``enable_stackglyembed``
(``False`` por defecto -- comportamiento identico al de antes de esta
mejora): si es ``True`` y ``Settings.STACKGLYEMBED_ENABLED``, para cada fila
con ``tipo_ptm`` en {``n_linked_glycosylation``, ``glycosylation_n``} (ver
``_NGLYCO_TYPES``) y ``pasa_umbral=True`` se invoca StackGlyEmbed
(``src/engines/stackglyembed_engine.py``) y se anaden 3 columnas PURAMENTE
INFORMATIVAS: ``stackglyembed_veredicto`` (``'Glicosilado'``/``'No
glicosilado'``), ``stackglyembed_score`` (su probabilidad cruda) y
``stackglyembed_coincide`` (``True`` si el veredicto es ``'Glicosilado'``).
En Camino FASTA esto NUNCA cambia ``pasa_umbral``/``consenso`` -- mismo
patron no-decisorio que MeToken (EMNGly no puede correr aqui, exige un PDB
real via MIF): un fallo deja las 3 columnas en ``None``.

## Consenso real de N-GLICOSILACION -- Camino PDB (EMNGly + StackGlyEmbed, decision 2026-08-06)

``annotate_pdb_path`` promueve StackGlyEmbed del mismo rol informativo de
arriba a MOTOR DE CONSENSO REAL, junto con EMNGly (``src/engines/emngly_engine.py``,
reemplazo de CoNglyPred -- confirmado sin pesos publicados en ningun sitio,
ver STATUS.md "Decision 2"), especificamente para las filas con ``tipo_ptm``
en ``_NGLYCO_TYPES`` y ``motor == 'DeepMVP'`` (candidatos propuestos por
DeepMVP que DeepPTMPred nunca fusiona para este tipo, ver
``CONSENSUS_EXCLUDED_TYPES`` mas abajo -- DeepPTMPred sigue reportandose
por separado, sin cambios). Se activa automaticamente cuando ``pdb_path`` no
es ``None`` y ``Settings.EMNGLY_ENABLED`` (mismo patron de activacion que
MeToken, sin parametro nuevo en la firma), independientemente del flag
``enable_stackglyembed`` -- en ese caso el pathway informativo de arriba se
omite para Camino PDB (no se invoca dos veces el mismo subproceso de
StackGlyEmbed, ver ``_apply_nglyco_consensus``).

Regla (provisional, ver ``Settings.NGLYCO_CONSENSUS_MIN_ENGINES``): de los
motores que lograron evaluar la posicion (DeepMVP siempre disponible;
EMNGly/StackGlyEmbed degradan a ausentes sin lanzar si no estan instalados,
ver sus respectivos ``_validate_installation``), ``pasa_umbral`` = al menos
1 pasa su propio umbral (generaliza la regla OR de 2 motores ya usada para
el resto de tipos), ``consenso`` = al menos
``Settings.NGLYCO_CONSENSUS_MIN_ENGINES`` (default 2) pasan. ``motor`` se
reescribe para reflejar exactamente que motores lograron evaluar la
posicion (p. ej. ``'DeepMVP+EMNGly+StackGlyEmbed'``, o
``'DeepMVP+StackGlyEmbed'`` si EMNGly esta degradado). Si ninguno de los 2
motores nuevos logra evaluar la posicion, el comportamiento es identico al
de antes de esta mejora (``motor='DeepMVP'``, ``consenso`` siempre
``False`` -- imposible alcanzar el minimo con un unico motor).

## Corroboracion opcional de VIA SECRETORA -- N-glicosilacion, ambos caminos
   (analisis de coherencia biologica 2026-08-07)

Ni DeepMVP ni DeepPTMPred ni EMNGly ni StackGlyEmbed modelan la via
biosintetica real del sustrato -- todos predicen desde secuencia/estructura
si un sequon N-X-[S/T] ES glicosilable, nunca si la proteina realmente
transita el RE/Golgi (requisito quimico real de la N-glicosilacion). Para
cada fila con ``tipo_ptm`` en ``_NGLYCO_TYPES`` y ``pasa_umbral=True`` (un
sitio que el consenso YA acepto), se consulta UniProt
(``src/structural/uniprot_localization_client.py``) UNA VEZ por accession y
se anade la columna PURAMENTE INFORMATIVA ``via_secretora_evidencia``:
``True`` (UniProt reporta localizacion consistente con la via secretora),
``False`` (UniProt SI tiene datos reales de localizacion pero NINGUNO es
consistente -- evidencia real en contra, no solo ausencia), o ``None`` (sin
datos de localizacion en UniProt, accession no reconocido -- el caso mas
comun aqui, ya que ``accession`` normalmente viene del nombre del archivo de
entrada, no de un ID UniProt real -- o fallo de red). NUNCA cambia
``pasa_umbral``/``consenso`` -- mismo patron no-decisorio que MeToken/
StackGlyEmbed-informativo. Deliberadamente NO aplica a
``o_linked_glycosylation``: existen dos vias de O-glicosilacion
biologicamente distintas (O-GlcNAc citoplasmatica/nuclear vs O-glicosilacion
de tipo mucina en la via secretora) y este cliente no las distingue -- alcance
limitado a N-glicosilacion para no arriesgar una afirmacion biologica
incorrecta.

## Aviso de COMPETENCIA/CROSSTALK entre PTMs del mismo residuo (analisis de
   coherencia biologica 2026-08-07)

Varios grupos de tipos de PTM de este proyecto modifican quimicamente el
MISMO grupo funcional de un residuo y son mutuamente excluyentes en una
misma molecula en un mismo instante -- no pueden coexistir literalmente
sobre el mismo atomo. El nucleo puntua cada tipo/posicion de forma
independiente, asi que si el consenso acepta 2+ tipos en competencia
exactamente en la misma posicion, hoy se presentan como hallazgos
igualmente validos sin ninguna senal de que son excluyentes. Grupos
reales, conservadores (ver ``_PTM_COMPETITION_GROUPS`` -- basados en
literatura establecida de PTM crosstalk, no exhaustivos a proposito, cada
grupo ademas se verifica contra el ``residuo_wt`` real de la fila, nunca
se asume solo por el tipo):
- Lisina (acilo-lisina): acetilacion, ubiquitinacion, sumoilacion,
  metilacion de Lys, malonilacion, glutarilacion, succinilacion,
  crotonilacion -- todas ocupan el mismo grupo epsilon-amino.
- Cisteina (tiol): S-nitrosilacion, glutationilacion.
- Arginina (guanidino): metilacion de Arg, citrulinacion (la citrulinacion
  ademas puede bloquear biologicamente la metilacion previa -- crosstalk
  documentado en la literatura de PAD/PRMT).
- Serina/Treonina (hidroxilo): fosforilacion, O-glicosilacion (la hipotesis
  "Yin-Yang" fosfo/O-GlcNAc, ampliamente documentada).

``_add_ptm_crosstalk_flag`` (ver mas abajo) anade la columna PURAMENTE
INFORMATIVA ``ptm_crosstalk_aviso`` a CUALQUIER fila con ``pasa_umbral=True``
que comparta grupo+posicion+accession con al menos otra fila tambien
``pasa_umbral=True`` -- NUNCA cambia ``pasa_umbral``/``consenso``. Se
ejecuta en ambos caminos (FASTA y PDB), tras todas las demas
corroboraciones, para ver el estado final de ``pasa_umbral`` de cada fila
(p. ej. tras el consenso de N-glicosilacion).

## Corroboracion opcional de ESPECIFICIDAD DE QUINASA -- fosforilacion, ambos
   caminos (analisis de coherencia biologica 2026-08-07, punto 5)

Ni DeepMVP ni DeepPTMPred distinguen QUE familia de quinasa fosforila un
sitio -- ambos predicen "fosforilable en general". A diferencia del tipo de
cadena de poliubiquitina (un evento celular posterior, no derivable de la
estructura -- razonamiento que vivia en Fase A/3c, eliminada del alcance
2026-08-10), la especificidad de quinasa SI es una propiedad local de
secuencia alrededor del sitio -- existe
una fuente real, publicada y descargable (Johnson et al. 2023 Nature +
Yaron-Barir et al. 2024 Nature, empaquetadas en ``kinase-library``, ver
docstring de ``src/engines/_kinase_library_runner.py``). Para cada fila con
``tipo_ptm`` en ``_PHOSPHO_TYPES`` y ``pasa_umbral=True``, se anaden 4
columnas PURAMENTE INFORMATIVAS: ``kinase_library_top_kinase``,
``kinase_library_top_family``, ``kinase_library_percentile`` y
``kinase_library_top3_kinases`` -- NUNCA cambian ``pasa_umbral``/``consenso``.
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from src.config.settings import Settings
from src.engines.emngly_engine import get_emngly_predictions
from src.engines.kinase_library_engine import get_kinase_corroboration
from src.engines.metoken_engine import get_type_corroboration
from src.engines.stackglyembed_engine import get_nglyco_corroboration
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
# 2026-08-01, ver STATUS.md "investigacion de n_linked_glycosylation y los 4
# tipos mediocres"): DeepPTMPred no tiene poder discriminativo real para
# estos tipos (verificado contra las metricas de entrenamiento del propio
# repo, AUC 0.495 -- no es un problema de umbral). Ambos motores se siguen
# reportando, nunca fusionados en una fila de consenso entre ELLOS DOS. Esto
# NO significa que el tipo se quede sin consenso posible: desde 2026-08-06,
# 'n_linked_glycosylation'/'glycosylation_n' en Camino PDB tiene un consenso
# real DISTINTO (DeepMVP+EMNGly+StackGlyEmbed, ver
# ``_apply_nglyco_consensus`` mas abajo) que reemplaza el rol que hubiera
# tenido DeepPTMPred si no estuviera muerto para este tipo.
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


def annotate_fasta_path(
    accession: str, sequence: str, deepmvp_df: pd.DataFrame, enable_stackglyembed: bool = False,
) -> pd.DataFrame:
    """Fase 3 nucleo (B), Camino FASTA: solo DeepMVP, sin consenso posible.

    Args:
        accession: Accession de la proteina (misma que ``deepmvp_df['protein']``).
        sequence: Secuencia saneada de Fase 1 (para calcular ventanas).
        deepmvp_df: Salida cruda de ``DeepMVPEngine.run()`` para este accession
            (columnas ``protein|aa|pos|x|y_pred|fpr|ptm``).
        enable_stackglyembed: ``False`` (default) desactiva la corroboracion
            de N-glicosilacion via StackGlyEmbed sin cambiar ningun otro
            comportamiento -- identico a antes de esta mejora. ``True``
            (y ``Settings.STACKGLYEMBED_ENABLED``) la habilita, ver
            docstring del modulo.

    Returns:
        DataFrame con ``OUTPUT_COLUMNS``, una fila por sitio candidato
        reportado por DeepMVP, sin ningun filtrado (ver
        :func:`apply_workflow_filter` para la capa D). Si
        ``enable_stackglyembed`` y ``Settings.STACKGLYEMBED_ENABLED``,
        incluye ademas ``stackglyembed_veredicto``/``stackglyembed_score``/
        ``stackglyembed_coincide`` -- ausentes (mismas ``OUTPUT_COLUMNS`` de
        siempre) en caso contrario.
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

    if enable_stackglyembed and Settings.STACKGLYEMBED_ENABLED:
        try:
            result = _add_stackglyembed_corroboration(result, sequence)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, stackglyembed_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion StackGlyEmbed para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
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


def annotate_pdb_path(
    accession: str, sequence: str, deepmvp_df: pd.DataFrame, deepptmpred_df: pd.DataFrame,
    pdb_path: Optional[Path] = None, chain_id: str = "A", enable_stackglyembed: bool = False,
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
        enable_stackglyembed: ``False`` (default) desactiva la corroboracion
            INFORMATIVA de N-glicosilacion via StackGlyEmbed sin cambiar
            ningun otro comportamiento -- identico a antes de esta mejora.
            Se ignora (no se invoca dos veces) cuando el consenso real de
            N-glicosilacion esta activo (``pdb_path`` no ``None`` y
            ``Settings.EMNGLY_ENABLED``, ver docstring del modulo) -- ese
            camino ya promueve StackGlyEmbed a motor de consenso.
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
        reflejan la fusion de hasta 3 motores, y se agregan ``score_emngly``/
        ``stackglyembed_veredicto``/``stackglyembed_score``/
        ``stackglyembed_coincide``. Si ese consenso NO esta activo pero
        ``enable_stackglyembed`` y ``Settings.STACKGLYEMBED_ENABLED``, se
        agregan las mismas 3 columnas de StackGlyEmbed pero puramente
        informativas (comportamiento identico a antes de esta mejora, ver
        docstring del modulo). Columnas ausentes (mismas ``OUTPUT_COLUMNS``
        de siempre) si ninguna condicion se cumple.
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
    # EMNGly/StackGlyEmbed si). Con el orden viejo (MeToken primero), esas filas
    # promovidas nunca recibian corroboracion MeToken -- quedaban en None
    # permanentemente pese a que la fila ya paso el consenso, contradiciendo la propia
    # documentacion de MeToken de ser "independiente de todo lo demas". Mismo criterio
    # ya aplicado a '_add_ptm_crosstalk_flag' (corre ultimo, tras ver el pasa_umbral
    # final de cada fila).
    nglyco_consensus_active = pdb_path is not None and Settings.EMNGLY_ENABLED
    if nglyco_consensus_active:
        try:
            result = _apply_nglyco_consensus(result, accession, sequence, pdb_path, position_mapping)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, emngly_engine/stackglyembed_engine ya degradan internamente
            logger.warning(
                "Fallo inesperado aplicando el consenso de N-glicosilacion (EMNGly+StackGlyEmbed) "
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

    if not nglyco_consensus_active and enable_stackglyembed and Settings.STACKGLYEMBED_ENABLED:
        try:
            result = _add_stackglyembed_corroboration(result, sequence)
        except Exception as exc:  # noqa: BLE001 -- doble seguro, stackglyembed_engine ya degrada internamente
            logger.warning(
                "Fallo inesperado anadiendo corroboracion StackGlyEmbed para '%s' (no fatal, no "
                "afecta consenso/pasa_umbral): %s", accession, exc,
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


def _add_stackglyembed_corroboration(result: pd.DataFrame, sequence: str) -> pd.DataFrame:
    """Anade ``stackglyembed_veredicto``/``stackglyembed_score``/``stackglyembed_coincide`` a las filas elegibles.

    Puramente informativo (ver docstring del modulo): NUNCA modifica
    ``pasa_umbral``/``consenso`` de ninguna fila, solo rellena estas 3
    columnas nuevas para las filas de N-glicosilacion (``tipo_ptm`` en
    ``_NGLYCO_TYPES`` -- cubre tanto el nombre crudo de DeepMVP
    ``glycosylation_n`` como el canonico ``n_linked_glycosylation`` de
    DeepPTMPred, ver ``annotate_pdb_path``) con ``pasa_umbral=True``. Si
    StackGlyEmbed no devuelve nada (venv/pickles no instalados, subproceso
    fallido, timeout -- ver ``stackglyembed_engine.get_nglyco_corroboration``,
    degrada sin lanzar), las 3 columnas quedan en ``None`` para toda la
    tabla.

    ``stackglyembed_coincide`` es ``True`` si el veredicto es
    ``'Glicosilado'`` (StackGlyEmbed corrobora el candidato que
    DeepMVP/DeepPTMPred ya propusieron), ``False`` si es
    ``'No glicosilado'`` (StackGlyEmbed discrepa) -- a diferencia de
    ``metoken_type_coincide``, no hay ambiguedad de "tipo sin equivalente"
    porque StackGlyEmbed solo predice N-glicosilacion, siempre comparable.
    """
    result = result.copy()
    result["stackglyembed_veredicto"] = None
    result["stackglyembed_score"] = None
    result["stackglyembed_coincide"] = None

    eligible_mask = result["tipo_ptm"].isin(_NGLYCO_TYPES) & result["pasa_umbral"]
    eligible_positions = sorted(set(result.loc[eligible_mask, "posicion"]))
    if not eligible_positions:
        return result

    corroboration = get_nglyco_corroboration(sequence, eligible_positions)
    if not corroboration:
        return result

    for idx, row in result.loc[eligible_mask].iterrows():
        site = corroboration.get(int(row["posicion"]))
        if site is None:
            continue
        result.at[idx, "stackglyembed_veredicto"] = site["stackglyembed_veredicto"]
        result.at[idx, "stackglyembed_score"] = site["stackglyembed_score"]
        result.at[idx, "stackglyembed_coincide"] = bool(site["stackglyembed_veredicto"] == "Glicosilado")

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

    A diferencia de ``stackglyembed_coincide``, no hay ninguna columna
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

    A diferencia de ``_add_metoken_corroboration``/``_add_stackglyembed_corroboration``
    (puramente informativos, evaluan solo filas YA con ``pasa_umbral=True``),
    esta funcion DECIDE ``pasa_umbral``/``consenso``/``motor`` -- evalua TODAS
    las filas nglyco propuestas por DeepMVP (``motor == 'DeepMVP'``, el unico
    motor de origen posible para este tipo, ver ``CONSENSUS_EXCLUDED_TYPES``),
    no solo las que ya pasaban el umbral de DeepMVP en solitario.

    No-op (devuelve ``result`` sin tocar, sin agregar columnas) si no hay
    ninguna fila nglyco elegible -- evita ensuciar el esquema de reportes que
    no tienen ningun candidato de N-glicosilacion.
    """
    nglyco_mask = result["tipo_ptm"].isin(_NGLYCO_TYPES) & (result["motor"] == "DeepMVP")
    if not nglyco_mask.any():
        return result

    result = result.copy()
    result["score_emngly"] = None
    result["stackglyembed_veredicto"] = None
    result["stackglyembed_score"] = None
    result["stackglyembed_coincide"] = None

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
            "N-glicosilacion (no fatal, degrada a DeepMVP+StackGlyEmbed).", accession,
        )

    stackglyembed_results = {}
    if Settings.STACKGLYEMBED_ENABLED:
        stackglyembed_results = get_nglyco_corroboration(sequence, positions)

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

        sge_site = stackglyembed_results.get(pos)
        if sge_site is not None:
            engines_ran.append("StackGlyEmbed")
            sge_pasa = bool(sge_site["stackglyembed_veredicto"] == "Glicosilado")
            result.at[idx, "stackglyembed_veredicto"] = sge_site["stackglyembed_veredicto"]
            result.at[idx, "stackglyembed_score"] = sge_site["stackglyembed_score"]
            result.at[idx, "stackglyembed_coincide"] = sge_pasa
            if sge_pasa:
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
