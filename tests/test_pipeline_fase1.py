"""Tests de integracion ligera del orquestador (pipeline.py) para Fase 1/1.5.

Cubre unicamente el enrutamiento + saneamiento/extraccion de hoy (Fase 3,
DeepMVP/DeepPTMPred, todavia no implementada). Sin mocks de motores porque
todavia no hay ningun motor que mockear.
"""

from pipeline import main

FASTA_CONTENT = ">ACC1 test protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n"

PDB_CONTENT = (
    "HEADER    TEST\n"
    "ATOM      1  N   MET A   1      11.104  13.207   2.100  1.00 20.00           N\n"
    "ATOM      2  CA  MET A   1      12.560  13.207   2.100  1.00 20.00           C\n"
    "ATOM      3  N   GLY A   2      14.500  14.700   2.100  1.00 20.00           N\n"
    "ATOM      4  CA  GLY A   2      15.000  15.700   2.100  1.00 20.00           C\n"
    "END\n"
)


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return p


def test_camino_fasta_produce_clean_fasta(tmp_path):
    input_path = _write(tmp_path, "ACC1.fasta", FASTA_CONTENT)
    output_dir = tmp_path / "out"

    exit_code = main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "ACC1_clean.fasta").is_file()


def test_camino_pdb_produce_fasta_derivado_y_position_mapping(tmp_path):
    input_path = _write(tmp_path, "1abc.pdb", PDB_CONTENT)
    output_dir = tmp_path / "out"

    exit_code = main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "1abc_derived.fasta").is_file()
    assert (output_dir / "1abc_position_mapping.csv").is_file()
    assert (output_dir / "1abc_chain_A.pdb").is_file()


def test_input_invalido_retorna_codigo_de_error(tmp_path):
    input_path = _write(tmp_path, "bad.fasta", "no es fasta valido\n")
    output_dir = tmp_path / "out"

    exit_code = main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 1
