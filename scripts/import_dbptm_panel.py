"""Importa el panel de validacion biologica desde dbPTM (decision 2026-08-09, vault).

## Por que existe

`src/validation/biological_panel.py` era un panel de 7 proteinas curadas a mano,
cada sitio con su PMID verificado uno por uno -- representativo del rigor del
proyecto (nunca fabricar datos cientificos), pero chico. Este script construye
una base de datos GRANDE de sitios PTM humanos reales derivada de dbPTM
(Academia Sinica/CUHK), que alimenta la columna "Literatura" del CLI
(`pipeline.py::_ground_truth_lookup`, via `src/validation/dbptm_lookup.py`) para
CUALQUIER proteina humana, no solo las 7-31 curadas a mano que ademas corren el
pipeline completo para medir recall (ver `scripts/select_recall_subset.py`).

## Hechos verificados en vivo, no asumidos de la documentacion (2026-08-09)

- Descarga real: `https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/{type}.gz`
  (el mirror `awi.cuhk.edu.cn` devolvio HTTP 403 a fetch automatizado -- descartado).
- Cada `.gz` es un tar de UN SOLO miembro (`tarfile.open(mode="r:gz")`, no
  `gzip.open` directo), sin cabecera, 6 columnas tab-delimited:
  `nombre_proteina, accession, posicion(1-based), tipo_dbptm, pmids(';'-separados),
  ventana_secuencia(21 chars, centro=indice 10 0-based, '-'=relleno de extremo)`.
  Formato confirmado identico en los 16 archivos relevantes (los 16, no una
  muestra -- ver verificacion real hecha durante la implementacion).
- dbPTM solo tiene "Methylation" generico (no separa Lys/Arg) -- se resuelve
  mirando el residuo en el centro de la ventana (K -> lys_methylation,
  R -> arg_methylation), sin necesitar el PDB.
- Filtrado a humano *reviewed*: ~29% de las filas sobreviven (verificado real
  contra Phosphorylation: 462965/1615054), via
  `src/validation/uniprot_human_accessions.py`.
- **Hallazgo real durante la primera corrida completa (2026-08-09)**: la
  columna de PMIDs no siempre son solo enteros separados por ';'. Tokens no
  numericos reales encontrados: `"UniProtKB CARBOHYD"` (23330 apariciones,
  sobre todo en N/O-linked Glycosylation -- una referencia a la propia
  anotacion curada de UniProt, no un PMID), `"-"` (sin referencia), `"N.N."`,
  `"doi:..."` (DOI en vez de PMID), `"<numero>PubMed"` (122 casos, PMID real
  con el sufijo "PubMed" pegado -- se recupera quitando el sufijo) y
  `"<numero>?"` (3 casos, el propio dbPTM marca la cita como incierta -- se
  descarta el token, no se recupera, para no tratar una cita dudosa como
  confirmada). Una fila puede mezclar PMIDs reales con tokens no numericos en
  la misma lista -- se filtra token por token, la fila solo se descarta si
  NINGUN token sobrevive (antes se descartaba la fila entera si CUALQUIER
  token no era numerico, perdiendo PMIDs reales -- bug real, corregido).

Corre en el venv principal del pipeline (`cnb_pipeline`, stdlib urllib/sqlite3,
sin dependencia nueva).
"""

import io
import json
import re
import sqlite3
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, NamedTuple, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.config.settings import Settings
from src.validation.uniprot_human_accessions import load_or_fetch_human_reviewed_accessions

DBPTM_BASE_URL = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/{type}.gz"
REQUEST_TIMEOUT_SECONDS = 120  # Phosphorylation.gz solo pesa 33MB comprimido, confirmado real
WINDOW_CENTER_INDEX = 10  # confirmado real: ventana de 21 chars, centro = indice 10 (0-based)

# 16 archivos dbPTM a descargar. 15 mapean 1:1 a un tipo canonico; "Methylation"
# cubre 2 (lys_methylation/arg_methylation), separados por el residuo en el
# centro de la ventana. Nombres de columna 4 confirmados reales byte-a-byte
# 2026-08-09 contra los 16 archivos descargados (no una muestra).
DBPTM_TYPE_TO_CANONICAL: Dict[str, Tuple[str, ...]] = {
    "Phosphorylation": ("phosphorylation",),
    "Acetylation": ("acetylation",),
    "Ubiquitination": ("ubiquitination",),
    "Hydroxylation": ("hydroxylation",),
    "Gamma-carboxyglutamic acid": ("gamma_carboxyglutamic_acid",),
    "Malonylation": ("malonylation",),
    "Crotonylation": ("crotonylation",),
    "Succinylation": ("succinylation",),
    "Glutathionylation": ("glutathionylation",),
    "Sumoylation": ("sumoylation",),
    "S-nitrosylation": ("s_nitrosylation",),
    "Glutarylation": ("glutarylation",),
    "Citrullination": ("citrullination",),
    "O-linked Glycosylation": ("o_linked_glycosylation",),
    "N-linked Glycosylation": ("n_linked_glycosylation",),
    "Methylation": ("lys_methylation", "arg_methylation"),
}


class DbptmRow(NamedTuple):
    accession: str  # canonicalizado, sin sufijo de isoforma
    position: int
    dbptm_type: str
    pmids: Tuple[int, ...]
    residue: str


class LookupSite(NamedTuple):
    residue: str
    tier: str  # "A" (2+ PMIDs independientes) o "B" (1 PMID) -- misma semantica que biological_panel.py
    pmids: Tuple[int, ...]


def _canonicalize_accession(accession: str) -> str:
    """Quita sufijo de isoforma ('-2' etc), mismo criterio que uniprot_localization_client.py."""
    return accession.split("-")[0]


def _download_type_file(dbptm_type: str, dest: Path) -> bool:
    """Descarga+extrae el archivo de ``dbptm_type`` a ``dest`` (texto plano, sin tar/gz).

    Idempotente (skip si ``dest.is_file()``). Nunca lanza -- reporta y devuelve
    False, mismo idioma que ``prepare_emngly_nglyde_structures.py::_download_alphafold_model``.
    """
    if dest.is_file():
        return True
    url = DBPTM_BASE_URL.format(type=urllib.parse.quote(dbptm_type))
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
            raw_gz = resp.read()
        with tarfile.open(fileobj=io.BytesIO(raw_gz), mode="r:gz") as tf:
            members = tf.getnames()
            if not members:
                print(f"[import_dbptm] {dbptm_type}: archivo tar vacio, se omite.")
                return False
            member = tf.extractfile(members[0])
            dest.write_bytes(member.read())
        return True
    except urllib.error.HTTPError as exc:
        print(f"[import_dbptm] {dbptm_type}: HTTP {exc.code}, se omite.")
        return False
    except urllib.error.URLError as exc:
        print(f"[import_dbptm] {dbptm_type}: error de red ({exc.reason}), se omite.")
        return False
    except tarfile.TarError as exc:
        print(f"[import_dbptm] {dbptm_type}: tar invalido ({exc}), se omite.")
        return False


_PMID_WITH_SUFFIX_RE = re.compile(r"^(\d+)PubMed$")


def _parse_pmids(pmids_str: str) -> Tuple[int, ...]:
    """Extrae los PMIDs numericos reales de ``pmids_str``, token por token.

    dbPTM mezcla PMIDs reales con tokens no numericos en la misma lista
    ';'-separada (ver hallazgo real documentado en el docstring del modulo):
    enteros puros se toman tal cual; ``"<numero>PubMed"`` se recupera quitando
    el sufijo (formato, no una referencia distinta); cualquier otra cosa
    (``"-"``, ``"N.N."``, ``"doi:..."``, ``"UniProtKB CARBOHYD"``,
    ``"<numero>?"``) se descarta -- deliberado para ``"<numero>?"`` tambien:
    el propio dbPTM marca esa cita como incierta, no se trata como confirmada.
    """
    pmids = set()
    for token in pmids_str.split(";"):
        token = token.strip()
        if not token:
            continue
        if token.isdigit():
            pmids.add(int(token))
            continue
        match = _PMID_WITH_SUFFIX_RE.match(token)
        if match:
            pmids.add(int(match.group(1)))
    return tuple(sorted(pmids))


def _iter_rows(raw_path: Path, dbptm_type: str) -> Iterator[DbptmRow]:
    """Parsea ``raw_path`` linea por linea. Filas malformadas se saltan y cuentan
    (no deberian ocurrir segun el formato confirmado, pero no se asume ciegamente)."""
    malformed = 0
    with raw_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                malformed += 1
                continue
            _name, accession, position_str, file_type, pmids_str, window = fields
            if len(window) != 21 or window[WINDOW_CENTER_INDEX] == "-":
                malformed += 1
                continue
            try:
                position = int(position_str)
            except ValueError:
                malformed += 1
                continue
            pmids = _parse_pmids(pmids_str)
            if not pmids:
                malformed += 1
                continue
            yield DbptmRow(
                accession=_canonicalize_accession(accession),
                position=position,
                dbptm_type=file_type,
                pmids=pmids,
                residue=window[WINDOW_CENTER_INDEX],
            )
    if malformed:
        print(f"[import_dbptm] {dbptm_type}: {malformed} fila(s) malformada(s) descartada(s).")


def _split_methylation(residue: str) -> Optional[str]:
    """'K' -> lys_methylation, 'R' -> arg_methylation, otro residuo -> None (descartado)."""
    if residue == "K":
        return "lys_methylation"
    if residue == "R":
        return "arg_methylation"
    return None


def build_lookup(
    human_accessions: set, raw_dir: Path,
) -> Dict[Tuple[str, int, str], LookupSite]:
    """Descarga+parsea los 16 archivos, filtra a humano, agrega por
    (accession, posicion, tipo_canonico) unificando PMIDs, deriva tier.

    Si dos filas agregadas para la misma clave discrepan en residuo (no deberia
    pasar si la posicion es real, pero no se asume), se descarta la clave entera
    -- mejor omitir un sitio dudoso que servir un dato incoherente.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    residues_by_key: Dict[Tuple[str, int, str], set] = defaultdict(set)
    pmids_by_key: Dict[Tuple[str, int, str], set] = defaultdict(set)
    n_files_ok = 0

    for i, (dbptm_type, canonical_types) in enumerate(DBPTM_TYPE_TO_CANONICAL.items(), start=1):
        dest = raw_dir / f"{dbptm_type.replace('/', '_')}.tsv"
        if not _download_type_file(dbptm_type, dest):
            continue
        n_files_ok += 1
        n_human = 0
        for row in _iter_rows(dest, dbptm_type):
            if row.accession not in human_accessions:
                continue
            if canonical_types == ("lys_methylation", "arg_methylation"):
                canonical = _split_methylation(row.residue)
                if canonical is None:
                    continue
            else:
                canonical = canonical_types[0]
            key = (row.accession, row.position, canonical)
            residues_by_key[key].add(row.residue)
            pmids_by_key[key].update(row.pmids)
            n_human += 1
        print(f"[import_dbptm] {dbptm_type}: {n_human} fila(s) tras filtro humano.")
        time.sleep(0.2)  # cortesia con el servidor

    lookup: Dict[Tuple[str, int, str], LookupSite] = {}
    n_inconsistent = 0
    for key, residues in residues_by_key.items():
        if len(residues) != 1:
            n_inconsistent += 1
            continue
        pmids = tuple(sorted(pmids_by_key[key]))
        tier = "A" if len(pmids) >= 2 else "B"
        lookup[key] = LookupSite(residue=next(iter(residues)), tier=tier, pmids=pmids)

    if n_inconsistent:
        print(f"[import_dbptm] {n_inconsistent} clave(s) con residuo inconsistente, descartada(s).")
    print(f"[import_dbptm] {n_files_ok}/{len(DBPTM_TYPE_TO_CANONICAL)} archivo(s) procesado(s) OK.")
    return lookup


def write_sqlite(lookup: Dict[Tuple[str, int, str], LookupSite], db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.is_file():
        db_path.unlink()
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute(
            """
            CREATE TABLE dbptm_sites (
                accession TEXT NOT NULL,
                position  INTEGER NOT NULL,
                ptm_type  TEXT NOT NULL,
                residue   TEXT NOT NULL,
                tier      TEXT NOT NULL CHECK(tier IN ('A','B')),
                pmids     TEXT NOT NULL,
                PRIMARY KEY (accession, position, ptm_type)
            )
            """
        )
        conn.execute("CREATE INDEX idx_dbptm_accession ON dbptm_sites(accession)")
        conn.executemany(
            "INSERT INTO dbptm_sites VALUES (?, ?, ?, ?, ?, ?)",
            (
                (acc, pos, ptm, site.residue, site.tier, json.dumps(list(site.pmids)))
                for (acc, pos, ptm), site in lookup.items()
            ),
        )
        conn.commit()
    finally:
        conn.close()


def main() -> int:
    Settings.ensure_dirs()
    accessions_cache = Settings.DBPTM_DATA_DIR / "human_reviewed_accessions.txt"
    human_accessions = load_or_fetch_human_reviewed_accessions(accessions_cache)
    print(f"[import_dbptm] {len(human_accessions)} accession(es) humana(s) reviewed.")

    raw_dir = Settings.DBPTM_DATA_DIR / "raw"
    lookup = build_lookup(human_accessions, raw_dir)

    db_path = Settings.DBPTM_DATA_DIR / "lookup.sqlite3"
    write_sqlite(lookup, db_path)

    n_accessions = len({acc for acc, _, _ in lookup})
    print(f"[import_dbptm] listo: {len(lookup)} sitio(s) unico(s) en {n_accessions} accession(es) humana(s) -> {db_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
