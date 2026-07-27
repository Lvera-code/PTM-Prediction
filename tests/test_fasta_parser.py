"""Tests de Fase 1 (src/utils/fasta_parser.py). Logica 100% pura, sin I/O externo."""

import pytest

from src.utils.exceptions import FastaFormatError, InvalidSequenceError
from src.utils.fasta_parser import FastaRecord, load_and_sanitize, parse_fasta, sanitize_sequence, write_fasta


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return p


def test_parse_fasta_separa_multiples_registros(tmp_path):
    p = _write(tmp_path, "multi.fasta", ">A desc1\nMKT\nAYI\n>B desc2\nGGG\n")
    records = parse_fasta(p)
    assert records == [("A desc1", "MKTAYI"), ("B desc2", "GGG")]


def test_parse_fasta_sin_cabecera_lanza_error(tmp_path):
    p = _write(tmp_path, "bad.fasta", "MKTAYI\n")
    with pytest.raises(FastaFormatError):
        parse_fasta(p)


def test_parse_fasta_archivo_inexistente_lanza_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_fasta(tmp_path / "no_existe.fasta")


def test_sanitize_sequence_normaliza_mayusculas():
    upper, invalid = sanitize_sequence("mktayi")
    assert upper == "MKTAYI"
    assert invalid == []


def test_sanitize_sequence_detecta_residuos_no_canonicos():
    upper, invalid = sanitize_sequence("MKTXAYI")
    assert upper == "MKTXAYI"
    assert invalid == ["X"]


def test_load_and_sanitize_registro_valido(tmp_path):
    p = _write(tmp_path, "ok.fasta", ">ACC1 desc\nMKTAYI\n")
    records = load_and_sanitize(p)
    assert records == [FastaRecord(header="ACC1 desc", accession="ACC1", sequence="MKTAYI")]


def test_load_and_sanitize_rechaza_residuo_no_canonico(tmp_path):
    p = _write(tmp_path, "bad.fasta", ">ACC1\nMKTXAYI\n")
    with pytest.raises(FastaFormatError):
        load_and_sanitize(p)


def test_load_and_sanitize_descarta_registro_vacio_sin_detener_lote(tmp_path):
    p = _write(tmp_path, "mixed.fasta", ">EMPTY\n>ACC1\nMKTAYI\n")
    records = load_and_sanitize(p)
    assert [r.accession for r in records] == ["ACC1"]


def test_load_and_sanitize_todos_vacios_lanza_invalid_sequence_error(tmp_path):
    p = _write(tmp_path, "allempty.fasta", ">A\n>B\n")
    with pytest.raises(InvalidSequenceError):
        load_and_sanitize(p)


def test_load_and_sanitize_accession_duplicado_lanza_error(tmp_path):
    p = _write(tmp_path, "dup.fasta", ">ACC1 desc1\nMKTAYI\n>ACC1 desc2\nGGG\n")
    with pytest.raises(FastaFormatError):
        load_and_sanitize(p)


def test_load_and_sanitize_accession_con_slash_se_sanea(tmp_path):
    p = _write(tmp_path, "slash.fasta", ">ACC/1 desc\nMKTAYI\n")
    records = load_and_sanitize(p)
    assert records[0].accession == "ACC_1"


def test_write_fasta_solo_usa_primer_token_de_cabecera(tmp_path):
    records = [FastaRecord(header="ACC1 free text, with comma", accession="ACC1", sequence="MKTAYI")]
    out_path = tmp_path / "out.fasta"
    write_fasta(records, out_path)
    content = out_path.read_text()
    assert content.splitlines()[0] == ">ACC1"
    assert "free text" not in content


def test_write_fasta_envuelve_secuencia_larga(tmp_path):
    long_seq = "A" * 130
    records = [FastaRecord(header="ACC1", accession="ACC1", sequence=long_seq)]
    out_path = tmp_path / "out.fasta"
    write_fasta(records, out_path, line_width=60)
    lines = out_path.read_text().splitlines()
    assert lines[1:] == ["A" * 60, "A" * 60, "A" * 10]
