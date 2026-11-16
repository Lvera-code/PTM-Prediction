"""Phase 1: raw FASTA sanitation (FASTA path, DeepMVP as the only engine).

Its exclusive responsibility: to read a FASTA exactly as it arrives from the user,
split it into records (header, sequence) and produce a sanitized version
(upper case, without internal line breaks) ready for DeepMVP (the only engine
of the FASTA path; the PDB path uses ``src.utils.structure_parser``
instead, not this module).

Policy on non-canonical residues (relaxed 2026-07-27 to match DeepMVP's
real tolerance, verified by reading
``lib/PeptideEncode.py::encodePeptideOneHot`` of the repo
github.com/bzhanglab/DeepMVP directly, not an assumption): DeepMVP NEVER aborts
on a residue outside its alphabet (the corresponding ``exit(1)`` is
commented out in the real code). 'X' is encoded as a zero vector (no
signal, not an error); any other unrecognized character (the IUPAC
ambiguities Z/J, pyrrolysine O, gaps '-', stops '*', digits, etc.) is encoded
as a vector of 0.5 with a printed warning, and the rest of the prediction continues.
The alphabet it does recognize with its own (non-degraded) encoding is the
20 standard residues plus U (selenocysteine) and B (Asx ambiguity), confirmed in
the repo's ``letterDict``.

Consequently, this module NO longer fatally rejects any record for
non-canonical residues; it only reports it as a warning, letting the
record through without modifying its sequence (the principle of never
removing or substituting a residue is maintained: doing so would shift the position
numbering and could fuse residues that are not neighbours in the real protein).
The only consumer of the FASTA path is DeepMVP (decision 2026-07-26, an
accepted asymmetry), so matching its real tolerance exactly is correct here
and not a generic relaxation of the pipeline.
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
