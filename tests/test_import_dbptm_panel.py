"""Tests de scripts/import_dbptm_panel.py contra fixtures sinteticos chicos.

Nunca descarga los archivos reales de dbPTM (cientos de MB) en CI -- la
verificacion del formato real (16/16 archivos, 6 columnas tab-delimited) se
hizo a mano contra el servidor real al implementar, documentada en el
docstring del modulo.
"""

import sqlite3
from pathlib import Path
from unittest.mock import patch

from scripts.import_dbptm_panel import (
    DbptmRow,
    LookupSite,
    _canonicalize_accession,
    _iter_rows,
    _parse_pmids,
    _split_methylation,
    build_lookup,
    write_sqlite,
)


def test_canonicalize_accession_quita_sufijo_isoforma():
    assert _canonicalize_accession("P04637-2") == "P04637"
    assert _canonicalize_accession("P04637") == "P04637"


def test_parse_pmids_solo_enteros():
    assert _parse_pmids("111;222;333") == (111, 222, 333)


def test_parse_pmids_dedup_y_ordena():
    assert _parse_pmids("222;111;111") == (111, 222)


def test_parse_pmids_recupera_sufijo_pubmed_pegado():
    # Caso real encontrado en Ubiquitination.tsv 2026-08-09: '14563840PubMed'.
    assert _parse_pmids("14563840PubMed;111") == (111, 14563840)


def test_parse_pmids_descarta_tokens_no_pmid_pero_conserva_los_reales():
    # Caso real encontrado en N-linked Glycosylation.tsv 2026-08-09: una fila
    # mezcla PMIDs reales con 'UniProtKB CARBOHYD' -- el bug original
    # descartaba la fila ENTERA, perdiendo los PMIDs reales.
    assert _parse_pmids("17989695;19159218;UniProtKB CARBOHYD;16335952") == (
        16335952, 17989695, 19159218,
    )


def test_parse_pmids_descarta_cita_marcada_incierta_por_dbptm():
    # '<numero>?' (3 casos reales en Phosphorylation.tsv): el propio dbPTM
    # marca la cita como incierta -- se descarta, no se recupera como valida.
    assert _parse_pmids("18321849?") == ()


def test_parse_pmids_descarta_dash_doi_y_nn():
    assert _parse_pmids("-") == ()
    assert _parse_pmids("N.N.") == ()
    assert _parse_pmids("doi:10.1007/s13562-013-0225-7") == ()


def test_parse_pmids_todo_no_valido_devuelve_vacio():
    assert _parse_pmids("UniProtKB CARBOHYD;-;N.N.") == ()


def test_split_methylation_lys_arg_y_otro():
    assert _split_methylation("K") == "lys_methylation"
    assert _split_methylation("R") == "arg_methylation"
    assert _split_methylation("S") is None


def test_iter_rows_parsea_fila_real_p53(tmp_path):
    # Fila real verificada contra dbPTM/Phosphorylation.gz descargado 2026-08-09,
    # coincide con GroundTruthSite(15, "S", "phosphorylation", "A", ...) ya
    # existente en biological_panel.py.
    raw = tmp_path / "Phosphorylation.tsv"
    raw.write_text("P53_HUMAN\tP04637\t15\tPhosphorylation\t20123963;16288207\tQSDPSVEPPLSQETFSDLWKL\n")
    rows = list(_iter_rows(raw, "Phosphorylation"))
    assert len(rows) == 1
    row = rows[0]
    assert row == DbptmRow(
        accession="P04637", position=15, dbptm_type="Phosphorylation",
        pmids=(16288207, 20123963), residue="S",
    )


def test_iter_rows_descarta_columnas_de_mas_o_de_menos(tmp_path, capsys):
    raw = tmp_path / "bad.tsv"
    raw.write_text("solo\tcuatro\tcolumnas\tsin_mas\n")
    rows = list(_iter_rows(raw, "bad"))
    assert rows == []
    assert "malformada" in capsys.readouterr().out


def test_iter_rows_descarta_ventana_en_relleno_de_terminal(tmp_path):
    raw = tmp_path / "term.tsv"
    # Un caso real de relleno: el centro de la ventana (indice 10) nunca deberia
    # ser '-' (esa posicion siempre es el residuo real anotado), pero se chequea
    # de todas formas -- defensivo, no se asume ciegamente.
    raw.write_text("X\tP00000\t1\tAcetylation\t123\t---------------------\n")
    assert list(_iter_rows(raw, "Acetylation")) == []


def test_iter_rows_canonicaliza_accession_con_isoforma(tmp_path):
    raw = tmp_path / "iso.tsv"
    raw.write_text("X\tP04637-2\t15\tPhosphorylation\t20123963\tQSDPSVEPPLSQETFSDLWKL\n")
    rows = list(_iter_rows(raw, "Phosphorylation"))
    assert rows[0].accession == "P04637"


def test_build_lookup_agrega_pmids_y_deriva_tier(tmp_path):
    raw_dir = tmp_path / "raw"

    def fake_download(dbptm_type, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dbptm_type == "Phosphorylation":
            # misma (accession, posicion) aparece 2 veces con PMIDs distintos --
            # debe agregarse en un solo LookupSite con la union de PMIDs.
            dest.write_text(
                "A\tP04637\t15\tPhosphorylation\t111\tQSDPSVEPPLSQETFSDLWKL\n"
                "A\tP04637\t15\tPhosphorylation\t222\tQSDPSVEPPLSQETFSDLWKL\n"
                "A\tP99999\t5\tPhosphorylation\t333\tAAAAASAAAAAAAAAAAAAAA\n"
            )
        elif dbptm_type == "Methylation":
            dest.write_text(
                "A\tP04637\t20\tMethylation\t444\tAAAAAAAAAAKAAAAAAAAAA\n"
                "A\tP04637\t21\tMethylation\t555\tAAAAAAAAAARAAAAAAAAAA\n"
            )
        else:
            dest.write_text("")
        return True

    human_accessions = {"P04637"}  # P99999 debe quedar filtrada afuera
    with patch("scripts.import_dbptm_panel._download_type_file", side_effect=fake_download):
        lookup = build_lookup(human_accessions, raw_dir)

    assert lookup[("P04637", 15, "phosphorylation")] == LookupSite(
        residue="S", tier="A", pmids=(111, 222),
    )
    assert lookup[("P04637", 20, "lys_methylation")] == LookupSite(
        residue="K", tier="B", pmids=(444,),
    )
    assert lookup[("P04637", 21, "arg_methylation")] == LookupSite(
        residue="R", tier="B", pmids=(555,),
    )
    assert ("P99999", 5, "phosphorylation") not in lookup


def test_build_lookup_descarta_clave_con_residuo_inconsistente(tmp_path, capsys):
    raw_dir = tmp_path / "raw"

    def fake_download(dbptm_type, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dbptm_type == "Acetylation":
            # misma (accession, posicion, tipo) con residuo distinto -- no deberia
            # pasar con datos reales, pero se descarta la clave entera si pasa.
            dest.write_text(
                "A\tP04637\t10\tAcetylation\t111\tAAAAAAAAAAKAAAAAAAAAA\n"
                "A\tP04637\t10\tAcetylation\t222\tAAAAAAAAAARAAAAAAAAAA\n"
            )
        else:
            dest.write_text("")
        return True

    with patch("scripts.import_dbptm_panel._download_type_file", side_effect=fake_download):
        lookup = build_lookup({"P04637"}, raw_dir)

    assert ("P04637", 10, "acetylation") not in lookup
    assert "inconsistente" in capsys.readouterr().out


def test_write_sqlite_produce_esquema_y_filas_esperadas(tmp_path):
    lookup = {
        ("P04637", 15, "phosphorylation"): LookupSite(residue="S", tier="A", pmids=(111, 222)),
    }
    db_path = tmp_path / "lookup.sqlite3"
    write_sqlite(lookup, db_path)

    conn = sqlite3.connect(str(db_path))
    rows = conn.execute("SELECT accession, position, ptm_type, residue, tier, pmids FROM dbptm_sites").fetchall()
    conn.close()

    assert rows == [("P04637", 15, "phosphorylation", "S", "A", "[111, 222]")]
