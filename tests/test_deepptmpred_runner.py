"""Tests de src/engines/_deepptmpred_runner.py (solo las piezas testeables sin PyRosetta/ESM).

El runner importa torch/esm/pyrosetta de forma diferida (dentro de las
funciones que las usan), asi que el modulo en si se puede importar y probar
sin esas dependencias pesadas -- ver docstring del modulo.
"""

from src.engines._deepptmpred_runner import _esm_cache_path


def test_esm_cache_path_incluye_hash_de_secuencia(tmp_path):
    path_a = _esm_cache_path(tmp_path, "ACC1", "MKTAYIAKQRQ")
    path_b = _esm_cache_path(tmp_path, "ACC1", "MKTAYIAKQRQ")
    assert path_a == path_b  # misma secuencia -> misma cache (determinista)


def test_esm_cache_path_distinta_secuencia_mismo_accession_no_colisiona(tmp_path):
    """Regresion: STATUS.md - auditoria 2026-07-28, item 2.

    Antes la clave de cache era solo el accession -- una secuencia distinta
    bajo el mismo accession (p. ej. PDB actualizado con el mismo nombre de
    archivo) reutilizaba en silencio el embedding ESM viejo.
    """
    path_v1 = _esm_cache_path(tmp_path, "ACC1", "MKTAYIAKQRQ")
    path_v2 = _esm_cache_path(tmp_path, "ACC1", "MKTAYIAKQRQAAAA")

    assert path_v1 != path_v2
    assert path_v1.parent == path_v2.parent == tmp_path


def test_esm_cache_path_distinto_accession_misma_secuencia_no_colisiona(tmp_path):
    path_acc1 = _esm_cache_path(tmp_path, "ACC1", "MKTAYIAKQRQ")
    path_acc2 = _esm_cache_path(tmp_path, "ACC2", "MKTAYIAKQRQ")
    assert path_acc1 != path_acc2
