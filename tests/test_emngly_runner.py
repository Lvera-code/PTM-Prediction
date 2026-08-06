"""Tests de src/engines/_emngly_runner.py (solo las piezas testeables sin fair-esm/torch/MIF).

El runner importa fair-esm/torch/MIF/scikit-learn de forma diferida (dentro
de las funciones que las usan), asi que el modulo en si se puede importar y
probar sin esas dependencias pesadas -- mismo patron que
``test_deepptmpred_runner.py``.
"""

import pandas as pd

from src.engines._emngly_runner import _fasta_position_to_pdb_seqid, _sequence_cache_key


def test_sequence_cache_key_es_determinista():
    assert _sequence_cache_key("MKTAYIAKQRQ") == _sequence_cache_key("MKTAYIAKQRQ")


def test_sequence_cache_key_distinta_secuencia_no_colisiona():
    assert _sequence_cache_key("MKTAYIAKQRQ") != _sequence_cache_key("MKTAYIAKQRQAAAA")


def test_fasta_position_to_pdb_seqid_numeracion_continua(tmp_path):
    """Caso comun (AlphaFold2): pdb_seqid == fasta_position, sin huecos."""
    csv_path = tmp_path / "mapping.csv"
    pd.DataFrame({
        "accession": ["ACC1"] * 3, "chain_id": ["A"] * 3, "pdb_seqid": [1, 2, 3],
        "insertion_code": ["", "", ""], "fasta_position": [1, 2, 3], "residue_letter": ["M", "K", "N"],
    }).to_csv(csv_path, index=False)

    mapping = _fasta_position_to_pdb_seqid(csv_path)

    assert mapping == {1: 1, 2: 2, 3: 3}


def test_fasta_position_to_pdb_seqid_numeracion_con_offset_y_huecos(tmp_path):
    """Estructura cristalografica real: arranca en residuo 21, con un hueco (21,22,25 observados)."""
    csv_path = tmp_path / "mapping.csv"
    pd.DataFrame({
        "accession": ["ACC1"] * 3, "chain_id": ["A"] * 3, "pdb_seqid": [21, 22, 25],
        "insertion_code": ["", "", ""], "fasta_position": [1, 2, 3], "residue_letter": ["M", "K", "N"],
    }).to_csv(csv_path, index=False)

    mapping = _fasta_position_to_pdb_seqid(csv_path)

    # fasta_position 3 (3er residuo OBSERVADO, ATMSEQ) -> pdb_seqid real 25, no 3.
    assert mapping == {1: 21, 2: 22, 3: 25}
