"""Fase 1.5: extraccion de secuencia (ATMSEQ) y mapeo de posiciones desde una estructura.

Se invoca UNICAMENTE para archivos que ``src.utils.input_router`` clasifico
como ``"structure"`` (PDB o mmCIF) -- el Camino PDB del proyecto, unico que
habilita el consenso DeepMVP + DeepPTMPred (DeepPTMPred exige ``pdb_path``
obligatorio, sin modo solo-secuencia).

Responsabilidad exclusiva: parsear la estructura con ``gemmi`` (100% local,
sin llamadas de red), elegir UNA cadena de referencia, extraer su secuencia
realmente resuelta en los registros ATOM/HETATM (ATMSEQ, NO SEQRES: es la
secuencia sobre la que DeepMVP/DeepPTMPred van a predecir de verdad y la que
debe coincidir con lo que se muestra en el reporte final de Fase 3) y
construir la tabla de mapeo de posiciones PDB <-> FASTA derivado (necesaria
para que las posiciones que reporta Fase 3 sean interpretables contra la
estructura original).

Resolucion de residuos via CCD (``gemmi.find_tabulated_residue``): esto
resuelve automaticamente residuos modificados (MSE->M, SEP->S, TPO->T,
PTR->Y, CSO->C, etc. -- ``one_letter_code`` viene en MINUSCULA para residuos
no estandar con padre canonico, de ahi el ``.upper()``). Un residuo sin
codigo de una sola letra resoluble usa 'X'.

Deliberadamente NO se usa ``ResidueSpan.make_one_letter_sequence()``: ante
un nombre de residuo no reconocido por el CCD, esa funcion puede insertar
caracteres extra (p. ej. un '-' de relleno) que desalinean el conteo de
caracteres respecto al numero real de residuos del polimero -- inviable para
construir un mapeo de posiciones 1:1. En su lugar se itera
``Chain.get_polymer()`` residuo por residuo, resolviendo cada uno
individualmente y garantizando exactamente 1 caracter por residuo.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List

import gemmi
import pandas as pd

from src.config.settings import Settings
from src.utils.exceptions import StructureParsingError
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

POSITION_MAPPING_COLUMNS = [
    "accession", "chain_id", "pdb_seqid", "insertion_code", "fasta_position", "residue_letter",
]


@dataclass(frozen=True)
class StructureRecord:
    """Resultado de Fase 1.5 para un archivo de estructura de entrada.

    Attributes:
        accession: Identificador derivado del nombre de archivo (``path.stem``),
            mismo rol que el accession de un registro FASTA en el resto del
            pipeline.
        pdb_path: Ruta a la estructura original, sin modificar (puede tener
            multiples cadenas). Es lo que se le pasa a DeepPTMPred (exige
            ``pdb_path`` obligatorio).
        chain_pdb_path: Ruta a un PDB derivado que contiene UNICAMENTE la
            cadena elegida como referencia.
        fasta_path: Ruta al FASTA canonico derivado (ATMSEQ de la cadena
            elegida), ya escrito a disco -- lo que recibe DeepMVP en Camino
            PDB (motor 1 de 2 del consenso).
        chain_id: Nombre de la cadena elegida como referencia.
        sequence: Secuencia ATMSEQ de esa cadena (1 caracter por residuo,
            mismo orden que ``position_mapping``).
        position_mapping: DataFrame con columnas
            ``POSITION_MAPPING_COLUMNS`` (una fila por residuo).
    """

    accession: str
    pdb_path: Path
    chain_pdb_path: Path
    fasta_path: Path
    chain_id: str
    sequence: str
    position_mapping: pd.DataFrame


def _resolve_residue_letter(resname: str) -> str:
    """Resuelve un nombre de residuo (codigo CCD de 3 letras) a su letra canonica.

    Usa ``gemmi.find_tabulated_residue`` (nunca ``None``: siempre devuelve un
    ``ResidueInfo``, con ``one_letter_code`` en blanco si el nombre no se
    reconoce). El codigo viene en mayuscula para residuos estandar y en
    minuscula para residuos modificados con padre canonico conocido (p. ej.
    MSE -> 'm'); se normaliza a mayuscula. Cualquier resultado que no sea una
    unica letra alfabetica cae a 'X'.
    """
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code.strip().upper() if info is not None else ""
    return code if len(code) == 1 and code.isalpha() else "X"


def _sanitize_accession(accession: str, source_name: str) -> str:
    """Sanea un accession derivado de ``path.stem`` (mismo criterio que ``fasta_parser.py``).

    ``path.stem`` nunca contiene separadores de ruta reales, pero un nombre
    de archivo literal como ``"..pdb"`` produce ``stem == ".."`` (verificado:
    ``Path("..pdb").stem == ".."``) -- si ese valor se usa despues para
    construir ``base_output_dir / record.accession`` (como hace
    ``DeepPTMPredEngine.run()``), el resultado escapa el directorio de salida
    (riesgo real de path traversal, ver STATUS.md - auditoria 2026-07-28,
    item 1). El accession de un registro FASTA ya sanea '/' y '\\'
    (``fasta_parser.py``); este helper aplica el mismo criterio aqui, mas el
    caso especifico de accession vacio/'.'/'..' que solo puede darse via
    nombre de archivo, no via cabecera FASTA.
    """
    sane = accession.replace("/", "_").replace("\\", "_")
    if sane in ("", ".", ".."):
        sane = "UNKNOWN"
    if sane != accession:
        logger.warning(
            "Accession '%s' derivado del nombre de archivo '%s' es inseguro para construir "
            "rutas de salida (vacio, '.', '..' o con separador de ruta): renombrado a '%s' "
            "por seguridad.",
            accession, source_name, sane,
        )
    return sane


def _select_chain(model: "gemmi.Model", strategy: str, explicit_chain_id: str) -> "gemmi.Chain":
    """Elige la cadena de referencia segun ``strategy``, siempre logueando el motivo.

    Solo se consideran cadenas con al menos 1 residuo en su polimero (via
    ``Chain.get_polymer()``, que usa la entidad/conectividad resuelta por
    gemmi para distinguir polipeptido de aguas/ligandos/iones): una cadena de
    solo aguas o solo heteroatomos no es una candidata valida.

    Raises:
        StructureParsingError: Si ninguna cadena tiene un polimero no vacio,
            o si ``strategy == 'explicit'`` y ``explicit_chain_id`` esta vacio
            o no existe / no es una cadena polimero valida.
    """
    candidates = [(chain, chain.get_polymer()) for chain in model]
    candidates = [(chain, poly) for chain, poly in candidates if poly.length() > 0]

    if not candidates:
        raise StructureParsingError(
            f"El modelo 1 de la estructura no tiene ninguna cadena con al menos un residuo de "
            f"aminoacido valido en su polimero (cadenas encontradas: {[c.name for c in model]}, "
            "todas vacias o compuestas solo por aguas/heteroatomos). No hay cadena de referencia "
            "posible: revisa que el archivo contenga una cadena proteica real."
        )

    if strategy == "explicit":
        if not explicit_chain_id:
            raise StructureParsingError(
                "PDB_CHAIN_SELECTION_STRATEGY='explicit' requiere PDB_EXPLICIT_CHAIN_ID "
                "configurado (no puede estar vacio)."
            )
        for chain, poly in candidates:
            if chain.name == explicit_chain_id:
                logger.info(
                    "Cadena '%s' elegida por estrategia 'explicit' (PDB_EXPLICIT_CHAIN_ID=%s), "
                    "%d residuo(s) en el polimero.",
                    chain.name, explicit_chain_id, poly.length(),
                )
                return chain
        raise StructureParsingError(
            f"PDB_EXPLICIT_CHAIN_ID='{explicit_chain_id}' no coincide con ninguna cadena "
            f"polimero valida (candidatas: {[c.name for c, _ in candidates]})."
        )

    if strategy != "longest":
        raise StructureParsingError(
            f"PDB_CHAIN_SELECTION_STRATEGY='{strategy}' no reconocida (valores validos: "
            "'longest', 'explicit')."
        )

    chosen_chain, chosen_poly = max(candidates, key=lambda item: item[1].length())
    others = ", ".join(f"{c.name}={p.length()}aa" for c, p in candidates if c.name != chosen_chain.name)
    logger.info(
        "Cadena '%s' elegida por estrategia 'longest' (%d residuo(s) en el polimero)%s.",
        chosen_chain.name, chosen_poly.length(),
        f"; otras candidatas: {others}" if others else " (unica cadena candidata)",
    )
    return chosen_chain


def parse_structure(
    path: Path,
    output_dir: Path,
    chain_selection_strategy: str = Settings.PDB_CHAIN_SELECTION_STRATEGY,
    explicit_chain_id: str = Settings.PDB_EXPLICIT_CHAIN_ID,
) -> StructureRecord:
    """Fase 1.5: parsea una estructura, extrae ATMSEQ de una cadena y persiste los artefactos.

    Args:
        path: Ruta al archivo PDB o mmCIF de entrada (ya clasificado como
            ``"structure"`` por ``src.utils.input_router``).
        output_dir: Carpeta donde se escriben el FASTA derivado y la tabla de
            mapeo de posiciones (mismo ``output_dir`` que usan las demas
            fases, p. ej. ``Outputs/``).
        chain_selection_strategy: ``'longest'`` (default) o ``'explicit'``.
        explicit_chain_id: ID de cadena a usar si la estrategia es
            ``'explicit'``.

    Returns:
        :class:`StructureRecord` con la cadena elegida, su secuencia ATMSEQ y
        el mapeo de posiciones, con los artefactos ya escritos en
        ``output_dir``.

    Raises:
        FileNotFoundError: Si ``path`` no existe.
        StructureParsingError: Si el archivo no puede parsearse, o si no hay
            ninguna cadena con un polimero de aminoacidos valido (ver
            :func:`_select_chain`).
    """
    if not path.is_file():
        raise FileNotFoundError(f"No se encontro el archivo de estructura de entrada: {path}")

    accession = _sanitize_accession(path.stem, path.name)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()
    except (RuntimeError, ValueError) as exc:
        raise StructureParsingError(
            f"No se pudo parsear '{path.name}' como PDB/mmCIF: {exc}. Verifica que el archivo "
            "tenga sintaxis valida (registros ATOM/HEADER para PDB legado, o un bloque 'data_' "
            "para mmCIF)."
        ) from exc

    if len(structure) == 0:
        raise StructureParsingError(f"'{path.name}' no contiene ningun modelo (MODEL) parseable.")

    model = structure[0]
    if len(structure) > 1:
        logger.info(
            "'%s' contiene %d modelos (probable estructura NMR/ensamble): se usa unicamente "
            "el modelo 1.",
            path.name, len(structure),
        )

    chain = _select_chain(model, chain_selection_strategy, explicit_chain_id)

    residues = list(chain.get_polymer())
    letters: List[str] = []
    mapping_rows: List[dict] = []
    for fasta_position, residue in enumerate(residues, start=1):
        letter = _resolve_residue_letter(residue.name)
        letters.append(letter)
        mapping_rows.append(
            {
                "accession": accession,
                "chain_id": chain.name,
                "pdb_seqid": residue.seqid.num,
                "insertion_code": residue.seqid.icode.strip() if residue.seqid.icode else "",
                "fasta_position": fasta_position,
                "residue_letter": letter,
            }
        )

    sequence = "".join(letters)
    n_unmapped = sum(1 for letter in letters if letter == "X")
    if n_unmapped:
        logger.warning(
            "'%s' (cadena '%s'): %d/%d residuo(s) no se pudieron resolver a una letra canonica "
            "via CCD, se reportan como 'X'.",
            path.name, chain.name, n_unmapped, len(letters),
        )

    fasta_path = output_dir / f"{accession}_derived.fasta"
    with fasta_path.open("w", encoding="utf-8") as fh:
        fh.write(f">{accession}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i : i + 60] + "\n")

    position_mapping = pd.DataFrame.from_records(mapping_rows, columns=POSITION_MAPPING_COLUMNS)
    mapping_path = output_dir / f"{accession}_position_mapping.csv"
    position_mapping.to_csv(mapping_path, index=False)

    chain_pdb_path = output_dir / f"{accession}_chain_{chain.name}.pdb"
    single_chain_structure = structure.clone()
    single_chain_model = single_chain_structure[0]
    for other_name in [c.name for c in single_chain_model if c.name != chain.name]:
        single_chain_model.remove_chain(other_name)
    single_chain_structure.write_pdb(str(chain_pdb_path))

    logger.info(
        "Fase 1.5 completa para '%s': accession='%s', cadena='%s', %d residuo(s) -> '%s' / '%s' / '%s'.",
        path.name, accession, chain.name, len(sequence), fasta_path, mapping_path, chain_pdb_path,
    )

    return StructureRecord(
        accession=accession,
        pdb_path=path,
        chain_pdb_path=chain_pdb_path,
        fasta_path=fasta_path,
        chain_id=chain.name,
        sequence=sequence,
        position_mapping=position_mapping,
    )
