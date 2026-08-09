"""Tests de src/validation/dbptm_lookup.py contra un sqlite fixture chico armado a mano."""

import sqlite3

from src.validation.dbptm_lookup import lookup_ground_truth


def _build_fixture_db(db_path):
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        """
        CREATE TABLE dbptm_sites (
            accession TEXT NOT NULL, position INTEGER NOT NULL, ptm_type TEXT NOT NULL,
            residue TEXT NOT NULL, tier TEXT NOT NULL, pmids TEXT NOT NULL,
            PRIMARY KEY (accession, position, ptm_type)
        )
        """
    )
    conn.executemany(
        "INSERT INTO dbptm_sites VALUES (?, ?, ?, ?, ?, ?)",
        [
            ("P04637", 15, "phosphorylation", "S", "A", "[111, 222]"),
            ("P04637", 382, "acetylation", "K", "B", "[333]"),
            ("P00000", 1, "phosphorylation", "S", "A", "[444, 555]"),
        ],
    )
    conn.commit()
    conn.close()


def test_lookup_devuelve_sitios_de_la_accession_pedida(tmp_path):
    db_path = tmp_path / "lookup.sqlite3"
    _build_fixture_db(db_path)
    assert lookup_ground_truth("P04637", db_path) == {
        (15, "phosphorylation"): ("A", (111, 222)),
        (382, "acetylation"): ("B", (333,)),
    }


def test_lookup_no_mezcla_accessions_distintas(tmp_path):
    db_path = tmp_path / "lookup.sqlite3"
    _build_fixture_db(db_path)
    assert lookup_ground_truth("P00000", db_path) == {(1, "phosphorylation"): ("A", (444, 555))}


def test_lookup_devuelve_vacio_para_accession_sin_datos(tmp_path):
    db_path = tmp_path / "lookup.sqlite3"
    _build_fixture_db(db_path)
    assert lookup_ground_truth("Q99999", db_path) == {}


def test_lookup_devuelve_vacio_si_no_existe_el_archivo(tmp_path):
    assert lookup_ground_truth("P04637", tmp_path / "no_existe.sqlite3") == {}
