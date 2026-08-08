"""Fase 1: Saneamiento de FASTA crudo (Camino FASTA, motor unico DeepMVP).

Responsabilidad exclusiva: leer un FASTA tal como llega del usuario,
separarlo en registros (cabecera, secuencia) y producir una version saneada
-mayusculas, sin saltos de linea internos- lista para DeepMVP (unico motor
del Camino FASTA; el Camino PDB usa ``src.utils.structure_parser`` en su
lugar, no este modulo).

Politica de residuos no canonicos (relajada 2026-07-27 para igualar la
tolerancia real de DeepMVP, verificada leyendo directamente
``lib/PeptideEncode.py::encodePeptideOneHot`` del repo
github.com/bzhanglab/DeepMVP -- no una suposicion): DeepMVP NUNCA aborta
ante un residuo fuera de su alfabeto (el ``exit(1)`` correspondiente esta
comentado en el codigo real). 'X' se codifica como vector de ceros (sin
señal, no es un error); cualquier otro caracter no reconocido (ambiguedades
IUPAC Z/J, pirrolisina O, gaps '-', stops '*', digitos, etc.) se codifica
como vector de 0.5 con un aviso impreso, y el resto de la prediccion sigue.
El alfabeto que SI reconoce con codificacion propia (no degradada) son los
20 estandar mas U (selenocisteina) y B (ambiguedad Asx), confirmado en
``letterDict`` del repo.

En consecuencia, este modulo ya NO rechaza (fatal) ningun registro por
residuos no canonicos -- solo lo reporta como warning, dejando pasar el
registro sin modificar su secuencia (se mantiene el principio de nunca
eliminar ni sustituir un residuo: eso desplazaria la numeracion de posicion
y podria fusionar residuos que en la proteina real no son vecinos). Unico
consumidor del Camino FASTA es DeepMVP (decision 2026-07-26, asimetria
aceptada), asi que igualar exactamente su tolerancia real es correcto aqui
y no una relajacion generica del pipeline.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from src.utils.exceptions import FastaFormatError, InvalidSequenceError
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

CANONICAL_AMINOACIDS = set("ACDEFGHIKLMNPQRSTVWY")
# Alfabeto que DeepMVP codifica con señal propia (no degradada), confirmado
# leyendo lib/PeptideEncode.py::letterDict: los 20 estandar mas U
# (selenocisteina) y B (ambiguedad Asx). 'X' se trata aparte (vector de
# ceros, ver sanitize_sequence): no es señal degradada, es "sin señal".
DEEPMVP_KNOWN_AMINOACIDS = CANONICAL_AMINOACIDS | {"U", "B"}


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
    """Normaliza a mayusculas y detecta (sin corregir) residuos que DeepMVP degrada.

    Deliberadamente NO elimina ni sustituye ningun caracter: ver politica
    documentada en el docstring del modulo. 'X' se excluye del reporte de
    "degradados" porque DeepMVP lo trata como "sin señal" (vector de ceros),
    no como un error -- solo se listan los caracteres que DeepMVP codifica
    como vector de 0.5 (fuera de los 20 estandar + U + B).

    Args:
        raw_sequence: Secuencia de aminoacidos sin procesar.

    Returns:
        Tupla ``(secuencia_en_mayusculas, caracteres_degradados_encontrados)``,
        con la lista de caracteres en orden de aparicion (vacia si la
        secuencia es 100% alfabeto conocido de DeepMVP o 'X').
    """
    upper = raw_sequence.upper()
    degraded_chars = [c for c in upper if c not in DEEPMVP_KNOWN_AMINOACIDS and c != "X"]
    return upper, degraded_chars


def load_and_sanitize(path: Path) -> List[FastaRecord]:
    """Lee y sanea un FASTA completo, descartando registros vacios y rechazando residuos no canonicos.

    Args:
        path: Ruta al archivo FASTA de entrada (dentro de ``inputs/``).

    Returns:
        Lista de :class:`FastaRecord` saneados.

    Raises:
        FastaFormatError: Si el archivo no tiene sintaxis FASTA valida
            (fatal); o si dos o mas registros comparten el mismo accession
            (primer token de la cabecera). Este ultimo caso es fatal por
            diseno: Fase 3 agrupa resultados por accession, y dos proteinas
            distintas con el mismo accession se fusionarian en una unica
            cadena de reporte. Los residuos no canonicos YA NO son fatales
            (ver politica relajada en el docstring del modulo) -- solo se
            reportan como warning.
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

        upper_seq, degraded_chars = sanitize_sequence(raw_seq)
        if degraded_chars:
            logger.warning(
                "Registro '%s' en '%s' contiene %d residuo(s) fuera del alfabeto conocido de "
                "DeepMVP (%s): se acepta igual, sin modificar la secuencia. DeepMVP los codifica "
                "como vector degradado (0.5) en vez de fallar (verificado leyendo "
                "lib/PeptideEncode.py::encodePeptideOneHot).",
                header, path.name, len(degraded_chars), sorted(set(degraded_chars)),
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
