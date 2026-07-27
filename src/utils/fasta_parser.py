"""Fase 1: Saneamiento de FASTA crudo (Camino FASTA, motor unico DeepMVP).

Responsabilidad exclusiva: leer un FASTA tal como llega del usuario,
separarlo en registros (cabecera, secuencia) y producir una version saneada
-mayusculas, sin saltos de linea internos- lista para DeepMVP (unico motor
del Camino FASTA; el Camino PDB usa ``src.utils.structure_parser`` en su
lugar, no este modulo).

Politica de residuos no canonicos: por defecto se RECHAZA (fatal) cualquier
caracter fuera de los 20 aminoacidos estandar, mismo motivo de fondo que en
``B-Cell-Epitope-Prediction``: eliminar o sustituir un residuo desplaza la
numeracion de posicion y puede fusionar residuos que en la proteina real no
son vecinos, fabricando una zona PTM en la costura que no existe en ninguna
secuencia real. Esta es una politica CONSERVADORA por defecto, no una
verificacion confirmada del comportamiento real de DeepMVP ante 'X'/'U'/'O'
-- a diferencia de BepiPred-3.0 en proyecto 1 (confirmado por lectura directa
del codigo que rechaza en bloque), la tolerancia real de DeepMVP a residuos
no canonicos todavia no se ha verificado leyendo su repo. Revisar y ajustar
esta politica cuando se construya ``src/engines/deepmvp_engine.py``.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from src.utils.exceptions import FastaFormatError, InvalidSequenceError
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

CANONICAL_AMINOACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class FastaRecord:
    """Un registro FASTA ya saneado, listo para Fase 3 (DeepMVP)."""

    header: str
    accession: str
    sequence: str


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    """Separa un archivo FASTA crudo en pares ``(cabecera, secuencia_cruda)``.

    Args:
        path: Ruta al archivo FASTA de entrada.

    Returns:
        Lista de tuplas ``(header_sin_'>', secuencia_sin_saltos_de_linea)``
        en el mismo orden en que aparecen en el archivo.

    Raises:
        FileNotFoundError: Si ``path`` no existe.
        FastaFormatError: Si el archivo no contiene ningun registro valido
            (no empieza con '>' o esta vacio). Es un error fatal.
    """
    if not path.is_file():
        raise FileNotFoundError(f"No se encontro el archivo FASTA de entrada: {path}")

    raw_text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line.strip() for line in raw_text.splitlines() if line.strip()]

    if not lines or not lines[0].startswith(">"):
        raise FastaFormatError(
            f"'{path.name}' no cumple la sintaxis FASTA minima (debe iniciar con '>')."
        )

    records: List[Tuple[str, str]] = []
    header: str = ""
    seq_chunks: List[str] = []

    for line in lines:
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq_chunks)))
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)

    if header:
        records.append((header, "".join(seq_chunks)))

    if not records:
        raise FastaFormatError(f"'{path.name}' no contiene ningun registro FASTA valido.")

    return records


def sanitize_sequence(raw_sequence: str) -> Tuple[str, List[str]]:
    """Normaliza a mayusculas y detecta (sin corregir) residuos no canonicos.

    Deliberadamente NO elimina ni sustituye los caracteres no canonicos: ver
    politica documentada en el docstring del modulo.

    Args:
        raw_sequence: Secuencia de aminoacidos sin procesar.

    Returns:
        Tupla ``(secuencia_en_mayusculas, caracteres_no_canonicos_encontrados)``,
        con la lista de caracteres invalidos en orden de aparicion (vacia si
        la secuencia es 100% canonica).
    """
    upper = raw_sequence.upper()
    invalid_chars = [c for c in upper if c not in CANONICAL_AMINOACIDS]
    return upper, invalid_chars


def load_and_sanitize(path: Path) -> List[FastaRecord]:
    """Lee y sanea un FASTA completo, descartando registros vacios y rechazando residuos no canonicos.

    Args:
        path: Ruta al archivo FASTA de entrada (dentro de ``fasta_inputs/``).

    Returns:
        Lista de :class:`FastaRecord` saneados.

    Raises:
        FastaFormatError: Si el archivo no tiene sintaxis FASTA valida
            (fatal); si algun registro contiene residuos no canonicos
            (fatal, ver politica en el docstring del modulo); o si dos o mas
            registros comparten el mismo accession (primer token de la
            cabecera). Este ultimo caso es fatal por diseno: Fase 3 agrupa
            resultados por accession, y dos proteinas distintas con el mismo
            accession se fusionarian en una unica cadena de reporte.
        InvalidSequenceError: Si NINGUN registro tiene secuencia (fatal a
            nivel de archivo). El descarte de registros individuales
            genuinamente vacios solo se loggea como warning y no detiene el
            resto del lote.
    """
    raw_records = parse_fasta(path)
    sane_records: List[FastaRecord] = []

    for header, raw_seq in raw_records:
        if not raw_seq:
            logger.warning("Registro '%s' descartado: no tiene ninguna secuencia asociada.", header)
            continue

        upper_seq, invalid_chars = sanitize_sequence(raw_seq)
        if invalid_chars:
            raise FastaFormatError(
                f"Registro '{header}' en '{path.name}' contiene {len(invalid_chars)} residuo(s) no "
                f"canonico(s) ({sorted(set(invalid_chars))}): rechazado por politica conservadora "
                "por defecto (ambiguedades IUPAC X/B/Z/J, selenocisteina U, pirrolisina O, gaps '-', "
                "stops '*', digitos, etc.). Sustituye manualmente ese residuo por su mejor "
                "aproximacion canonica (o elimina el registro) en el FASTA de entrada y vuelve a "
                "intentarlo."
            )

        accession = header.split()[0] if header else "UNKNOWN"
        if "/" in accession or "\\" in accession:
            sane_accession = accession.replace("/", "_").replace("\\", "_")
            logger.warning(
                "Accession '%s' en '%s' contiene un separador de ruta ('/' o '\\'): "
                "renombrado a '%s' por seguridad (evita construir rutas de archivo "
                "inesperadas al escribir outputs derivados de este accession).",
                accession, path.name, sane_accession,
            )
            accession = sane_accession
        sane_records.append(FastaRecord(header=header, accession=accession, sequence=upper_seq))

    if not sane_records:
        raise InvalidSequenceError(f"'{path.name}' no contiene ningun registro con secuencia.")

    accession_counts: dict = {}
    for record in sane_records:
        accession_counts[record.accession] = accession_counts.get(record.accession, 0) + 1
    duplicates = sorted(acc for acc, count in accession_counts.items() if count > 1)
    if duplicates:
        raise FastaFormatError(
            f"'{path.name}' contiene registros con accession duplicado: {duplicates}. "
            "Cada accession debe identificar una unica secuencia fisica (Fase 3 agrupa por "
            "accession); renombra las cabeceras duplicadas en el FASTA de entrada y vuelve a "
            "intentarlo."
        )

    return sane_records


def write_fasta(records: List[FastaRecord], out_path: Path, line_width: int = 60) -> None:
    """Escribe una lista de :class:`FastaRecord` saneados como FASTA valido.

    Escribe unicamente ``record.accession`` (primer token de la cabecera
    original, sin espacios) como cabecera, descartando el resto de la
    descripcion libre -- evita que un motor que solo conserve el primer
    token de la cabecera (patron comun) reporte un accession distinto al que
    usa el resto del pipeline.

    Args:
        records: Registros a escribir, en orden.
        out_path: Ruta de salida (se sobreescribe si ya existe).
        line_width: Ancho de linea para el envoltorio de la secuencia.
    """
    with out_path.open("w", encoding="utf-8") as fh:
        for record in records:
            fh.write(f">{record.accession}\n")
            seq = record.sequence
            for i in range(0, len(seq), line_width):
                fh.write(seq[i : i + line_width] + "\n")
