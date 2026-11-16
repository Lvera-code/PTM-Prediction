"""Selects the recall subset (~20-30 proteins) of the dbPTM panel.

## Why it exists

`src/validation/biological_panel.py::PANEL` is the SMALL subset that also
runs the complete pipeline (`scripts/validate_biological_panel.py`) to measure
real recall, which is why it cannot be the whole dbPTM database (~19000 proteins,
days of compute). This script systematically selects ~20-30 candidates
from `data/dbptm/lookup.sqlite3` (decision 2026-08-09, vault), replacing
the manual choice of the current 7 proteins.

## Method

1. Eligibility: at least 1 tier A site, real UniProt length <= 1500 aa
   (for PDB handling/download reasons, NOT for engine time, which does
   not scale strongly with length: DeepPTMPred dominates with 17 ESM-2 calls per
   protein, independent of size).
2. Greedy set-cover over the 17 canonical types: at each step it chooses the
   accession that covers the most types not yet covered with tier A evidence
   (tie-break: the most tier A sites in total). It keeps adding high-value
   candidates beyond the minimum coverage up to a target of ~25 (range
   20-30), giving redundancy to the most common types, just as the current
   panel already does with p53+histones.
3. `--emit-python`: for each candidate, it emits up to `MAX_SITES_PER_ACCESSION_TYPE`
   (5) sites per canonical type, prioritizing tier A and breaking ties by more
   PMIDs. Without this cap, heavily studied types (phosphorylation, ubiquitination)
   can bring hundreds of tier A sites for a single protein (a real case
   found 2026-08-09: MAP4/P27816 gave 549 candidate sites without the cap).
   The format is ready to paste into `biological_panel.py`.
4. `--verify`: it downloads the real AlphaFold PDB per candidate and checks the
   length plus the expected residue at each candidate site BEFORE accepting it.
   The project's hard rule documents 3 real numbering traps
   (histones, prothrombin, EPO) that this is designed to catch during
   curation, not later in CI.

**The output is a DRAFT, not final**: the project's hard rule (never
fabricate scientific data) requires that every site of the recall panel have its
PMID verified by title via NCBI eutils, which is stricter than "dbPTM says
so". This script resolves the systematic SELECTION; the verification of
evidence remains manual before pasting the result into
`biological_panel.py`.
"""

import argparse
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.config.settings import Settings
from src.utils.exceptions import StructureParsingError
from src.utils.structure_parser import parse_structure
from src.validation.dbptm_lookup import DBPTM_LOOKUP_PATH

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUTS_DIR = REPO_ROOT / "inputs"

ALPHAFOLD_API_URL = "https://www.alphafold.ebi.ac.uk/api/prediction/{acc}"
REQUEST_TIMEOUT_SECONDS = 30

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_QUERY = "organism_id:9606 AND reviewed:true"
UNIPROT_TIMEOUT_SECONDS = 120

# Solo manejo de PDB/descarga -- NO por tiempo de motor (ver docstring del modulo).
MAX_PROTEIN_LENGTH = 1500
TARGET_SUBSET_SIZE = 25
MAX_SUBSET_SIZE = 30
MAX_SITES_PER_ACCESSION_TYPE = 5
# Umbral para descartar una candidata entera tras --verify: si menos de la
# mitad de sus sitios candidatos coinciden con el PDB real, es mas probable
# una trampa de numeracion sistematica que datos puntuales malos.
MIN_VERIFIED_FRACTION = 0.5

CANONICAL_TYPES: Tuple[str, ...] = Settings.DEEPPTMPRED_PTM_TYPES


class DbptmLookupError(Exception):
    """La base dbPTM (data/dbptm/lookup.sqlite3) no existe todavia."""


def fetch_human_reviewed_lengths() -> Dict[str, Tuple[str, int]]:
    """``{accession: (entry_name, length)}`` para todo humano reviewed, un solo request."""
    url = (
        f"{UNIPROT_STREAM_URL}?query={urllib.parse.quote(UNIPROT_QUERY)}"
        "&fields=accession,id,length&format=tsv"
    )
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(request, timeout=UNIPROT_TIMEOUT_SECONDS) as response:
        body = response.read().decode("utf-8")
    result = {}
    for line in body.splitlines()[1:]:  # salta cabecera "Entry\tEntry Name\tLength"
        fields = line.split("\t")
        if len(fields) != 3:
            continue
        accession, entry_name, length_str = fields
        try:
            result[accession] = (entry_name, int(length_str))
        except ValueError:
            continue
    return result


def load_or_fetch_human_reviewed_lengths(cache_path: Path, force_refresh: bool = False) -> Dict[str, Tuple[str, int]]:
    if cache_path.is_file() and not force_refresh:
        result = {}
        for line in cache_path.read_text().splitlines():
            accession, entry_name, length_str = line.split("\t")
            result[accession] = (entry_name, int(length_str))
        return result
    lengths = fetch_human_reviewed_lengths()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(
        "\n".join(f"{acc}\t{name}\t{length}" for acc, (name, length) in sorted(lengths.items())) + "\n"
    )
    return lengths


def load_tier_a_types(db_path: Path) -> Dict[str, Set[str]]:
    """``{accession: {tipos canonicos con >=1 sitio tier A}}``."""
    if not db_path.is_file():
        raise DbptmLookupError(
            f"{db_path} no existe -- correr scripts/import_dbptm_panel.py primero."
        )
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = conn.execute("SELECT DISTINCT accession, ptm_type FROM dbptm_sites WHERE tier = 'A'").fetchall()
    finally:
        conn.close()
    result: Dict[str, Set[str]] = defaultdict(set)
    for accession, ptm_type in rows:
        result[accession].add(ptm_type)
    return result


def load_tier_a_counts(db_path: Path) -> Dict[str, int]:
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT accession, COUNT(*) FROM dbptm_sites WHERE tier = 'A' GROUP BY accession"
        ).fetchall()
    finally:
        conn.close()
    return dict(rows)


def select_recall_subset(
    tier_a_types: Dict[str, Set[str]],
    tier_a_counts: Dict[str, int],
    lengths: Dict[str, Tuple[str, int]],
    target_size: int = TARGET_SUBSET_SIZE,
    max_size: int = MAX_SUBSET_SIZE,
) -> Tuple[List[str], Set[str]]:
    """Greedy set-cover sobre CANONICAL_TYPES. Devuelve (seleccionadas, tipos_cubiertos)."""
    eligible = {
        acc for acc, types in tier_a_types.items()
        if types and acc in lengths and lengths[acc][1] <= MAX_PROTEIN_LENGTH
    }
    covered: Set[str] = set()
    selected: List[str] = []
    pool = set(eligible)

    while pool and len(selected) < max_size:
        best = max(pool, key=lambda acc: (len(tier_a_types[acc] - covered), tier_a_counts.get(acc, 0)))
        new_coverage = len(tier_a_types[best] - covered)
        if new_coverage == 0 and len(selected) >= target_size:
            break
        selected.append(best)
        covered |= tier_a_types[best]
        pool.discard(best)

    return selected, covered


def load_candidate_sites(accession: str, db_path: Path) -> List[Tuple[str, int, str, str, Tuple[int, ...]]]:
    """Sitios candidatos para ``accession``, hasta ``MAX_SITES_PER_ACCESSION_TYPE``
    por tipo canonico (preferidos tier A, desempate por mas PMIDs).

    Sin este tope, tipos muy estudiados (fosforilacion, ubiquitinacion) pueden
    traer cientos de sitios tier A para una sola proteina (caso real
    encontrado 2026-08-09: MAP4/P27816 daba 549 sitios candidatos) -- rompe el
    objetivo de un borrador chico y revisable. Devuelve tuplas
    ``(ptm_type, position, residue, tier, pmids)``.
    """
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT ptm_type, position, residue, tier, pmids FROM dbptm_sites WHERE accession = ? "
            "ORDER BY ptm_type, tier, position",
            (accession,),
        ).fetchall()
    finally:
        conn.close()

    by_type: Dict[str, List[Tuple[str, int, str, str, Tuple[int, ...]]]] = defaultdict(list)
    for ptm_type, position, residue, tier, pmids_json in rows:
        by_type[ptm_type].append((ptm_type, position, residue, tier, tuple(json.loads(pmids_json))))

    result: List[Tuple[str, int, str, str, Tuple[int, ...]]] = []
    for ptm_type, sites in by_type.items():
        ranked = sorted(sites, key=lambda s: (s[3] != "A", -len(s[4])))
        result.extend(ranked[:MAX_SITES_PER_ACCESSION_TYPE])
    return result


def _download_alphafold_model(accession: str, dest: Path) -> bool:
    """Idempotente, nunca lanza -- mismo idioma que
    ``prepare_emngly_nglyde_structures.py::_download_alphafold_model``."""
    if dest.is_file():
        return True
    try:
        req = urllib.request.Request(ALPHAFOLD_API_URL.format(acc=accession))
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
            entries = json.loads(resp.read())
        if not entries:
            print(f"[select_recall] {accession}: sin modelo AlphaFold (API vacia), se omite.")
            return False
        entry = next((e for e in entries if e.get("uniprotAccession") == accession), entries[0])
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            print(f"[select_recall] {accession}: API sin 'pdbUrl', se omite.")
            return False
        with urllib.request.urlopen(pdb_url, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
            dest.write_bytes(resp.read())
        return True
    except urllib.error.HTTPError as exc:
        print(f"[select_recall] {accession}: sin modelo AlphaFold ({exc.code}), se omite.")
        return False
    except urllib.error.URLError as exc:
        print(f"[select_recall] {accession}: error de red ({exc.reason}), se omite.")
        return False


def verify_candidate(
    accession: str, entry_name: str, sites: List[Tuple[str, int, str, str, Tuple[int, ...]]], output_dir: Path,
) -> Tuple[bool, List[Tuple[str, int, str, str, Tuple[int, ...]]], int, Path]:
    """Descarga el PDB real y verifica longitud + residuo esperado por sitio.

    Devuelve (aceptar_candidata, sitios_verificados, longitud_real_del_pdb,
    pdb_path). La longitud devuelta es SIEMPRE la de la secuencia ATMSEQ
    parseada del PDB real (no el largo canonico de UniProt) -- son casi
    siempre iguales pero es la primera la que
    ``tests/test_biological_panel.py::test_pdb_existe_y_tiene_la_longitud_esperada``
    exige que coincida, asi que es la unica correcta para el campo ``length``.
    Un sitio se descarta si el residuo no coincide con la secuencia real
    (posible trampa de numeracion puntual). La candidata entera se descarta si
    sobrevive menos de MIN_VERIFIED_FRACTION de sus sitios (trampa
    sistematica, ver docstring del modulo).
    """
    slug = entry_name.split("_")[0].lower()
    pdb_path = INPUTS_DIR / f"{slug}_{accession}.pdb"
    if not _download_alphafold_model(accession, pdb_path):
        return False, [], 0, pdb_path

    try:
        record = parse_structure(pdb_path, output_dir)
    except (StructureParsingError, FileNotFoundError) as exc:
        print(f"[select_recall] {accession}: fallo Fase 1.5 ({exc}), se descarta.")
        return False, [], 0, pdb_path

    verified = [
        site for site in sites
        if 1 <= site[1] <= len(record.sequence) and record.sequence[site[1] - 1] == site[2]
    ]
    fraction = len(verified) / len(sites) if sites else 0.0
    if fraction < MIN_VERIFIED_FRACTION:
        print(
            f"[select_recall] {accession}: solo {len(verified)}/{len(sites)} sitio(s) "
            f"coinciden con el PDB real ({fraction:.0%}) -- posible trampa de numeracion, se descarta."
        )
        return False, verified, len(record.sequence), pdb_path

    if len(verified) < len(sites):
        print(f"[select_recall] {accession}: {len(sites) - len(verified)} sitio(s) descartado(s) por no coincidir con el PDB real.")
    return True, verified, len(record.sequence), pdb_path


def emit_python(accession: str, entry_name: str, length: int, sites: List[Tuple[str, int, str, str, Tuple[int, ...]]]) -> str:
    slug = entry_name.split("_")[0].lower()
    lines = [f'# BORRADOR -- revisar cada PMID a mano antes de pegar (ver docstring del modulo).']
    lines.append(f'_{slug.upper()} = PanelEntry(')
    lines.append(f'    name="{slug}", uniprot_accession="{accession}", pdb_filename="{slug}_{accession}.pdb", length={length},')
    lines.append('    sites=(')
    for ptm_type, position, residue, tier, pmids in sites:
        pmids_repr = ", ".join(str(p) for p in pmids)
        lines.append(
            f'        GroundTruthSite({position}, "{residue}", "{ptm_type}", "{tier}", ({pmids_repr}{"," if len(pmids) == 1 else ""})),'
        )
    lines.append('    ),')
    lines.append(')')
    return "\n".join(lines)


def main(argv: List[str] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-size", type=int, default=TARGET_SUBSET_SIZE)
    parser.add_argument("--emit-python", action="store_true", help="Emite el borrador PanelEntry/GroundTruthSite.")
    parser.add_argument("--verify", action="store_true", help="Descarga PDB real de AlphaFold y verifica cada sitio.")
    parser.add_argument("--output", type=Path, default=None, help="Archivo de salida para --emit-python (default: stdout).")
    args = parser.parse_args(argv)

    Settings.ensure_dirs()

    try:
        tier_a_types = load_tier_a_types(DBPTM_LOOKUP_PATH)
    except DbptmLookupError as exc:
        print(f"[select_recall] {exc}", file=sys.stderr)
        return 1
    tier_a_counts = load_tier_a_counts(DBPTM_LOOKUP_PATH)
    lengths = load_or_fetch_human_reviewed_lengths(Settings.DBPTM_DATA_DIR / "human_reviewed_lengths.tsv")

    selected, covered = select_recall_subset(tier_a_types, tier_a_counts, lengths, target_size=args.target_size)
    uncovered = set(CANONICAL_TYPES) - covered
    print(f"[select_recall] {len(selected)} candidata(s) seleccionada(s), {len(covered)}/{len(CANONICAL_TYPES)} tipo(s) cubierto(s).")
    if uncovered:
        print(f"[select_recall] tipo(s) SIN cobertura tier A en ninguna candidata elegible: {sorted(uncovered)}")

    accepted_drafts = []
    output_dir = Settings.FASTA_OUTPUT_DIR / "select_recall_subset"
    for i, accession in enumerate(selected, start=1):
        entry_name, length = lengths[accession]
        sites = load_candidate_sites(accession, DBPTM_LOOKUP_PATH)
        print(f"[select_recall] ({i}/{len(selected)}) {accession} ({entry_name}, {length} aa): {len(sites)} sitio(s) candidato(s).")

        if args.verify:
            ok, sites, verified_length, _pdb_path = verify_candidate(accession, entry_name, sites, output_dir)
            time.sleep(0.2)  # cortesia con el servidor de AlphaFold DB
            if not ok:
                continue
            length = verified_length  # longitud real del PDB, no la de UniProt (ver docstring de verify_candidate)

        if args.emit_python:
            accepted_drafts.append(emit_python(accession, entry_name, length, sites))

    if args.emit_python:
        draft = "\n\n".join(accepted_drafts)
        if args.output:
            args.output.write_text(draft + "\n")
            print(f"[select_recall] borrador escrito en {args.output}")
        else:
            print("\n" + draft)

    return 0


if __name__ == "__main__":
    sys.exit(main())
