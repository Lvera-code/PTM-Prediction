"""Tests dedicados de src/structural/ubiquitin_sumo.py (punto 6 del plan de robustez
post-demo-prep). ``_validate_target`` no importa pyrosetta (solo lo hacen ``conjugate``/
``init_pyrosetta``/``_lyx_params_path``, siempre con imports locales) -- es pura logica
Python sobre un objeto ``pose``, testeable con un fake ligero sin pyrosetta instalado. Cubre
las 4 excepciones reales documentadas en el modulo (el docstring las llama "hallazgos
reales": Rosetta mismo NO valida que el residuo sea Lisina, ni que la numeracion PDB
coincida 1:1 con la de pose -- ``_validate_target`` es la unica red de seguridad real).
"""

import sys
import types

import pytest

from src.structural import ubiquitin_sumo


class _FakeResidue:
    def __init__(self, name1):
        self._name1 = name1

    def name1(self):
        return self._name1


class _FakeConformation:
    def __init__(self, num_chains):
        self._num_chains = num_chains

    def num_chains(self):
        return self._num_chains


class _FakePdbInfo:
    def __init__(self, numbering):
        self._numbering = numbering

    def number(self, pose_index):
        return self._numbering[pose_index]


class _FakePose:
    def __init__(self, size=10, num_chains=1, numbering=None, residues=None):
        self._size = size
        self._num_chains = num_chains
        self._numbering = numbering or {i: i for i in range(1, size + 1)}
        self._residues = residues or {}

    def size(self):
        return self._size

    def conformation(self):
        return _FakeConformation(self._num_chains)

    def pdb_info(self):
        return _FakePdbInfo(self._numbering)

    def residue(self, position):
        return _FakeResidue(self._residues.get(position, "A"))


def test_tipo_no_soportado():
    pose = _FakePose(residues={24: "K"})
    with pytest.raises(ubiquitin_sumo.UnsupportedConjugationTypeError):
        ubiquitin_sumo._validate_target(pose, 24, "sumoylation_falsa")


def test_pose_multicadena_rechazada():
    pose = _FakePose(num_chains=2, residues={24: "K"})
    with pytest.raises(ubiquitin_sumo.MultiChainPoseError):
        ubiquitin_sumo._validate_target(pose, 24, "ubiquitination")


def test_numeracion_no_secuencial_rechazada():
    # Residuo 5 en pose tiene numero PDB 999 -- no coincide 1:1, ver docstring del modulo
    # (segundo hallazgo: initialize() hace una doble conversion pose->pdb->pose).
    numbering = {i: i for i in range(1, 11)}
    numbering[5] = 999
    pose = _FakePose(numbering=numbering, residues={24: "K"})
    with pytest.raises(ubiquitin_sumo.NonSequentialNumberingError):
        ubiquitin_sumo._validate_target(pose, 24, "ubiquitination")


def test_residuo_no_lisina_rechazado():
    # Rosetta mismo no verifica esto (runtime_assert comentada en el .cc real, ver docstring).
    pose = _FakePose(residues={24: "A"})
    with pytest.raises(ubiquitin_sumo.ResidueMismatchError):
        ubiquitin_sumo._validate_target(pose, 24, "ubiquitination")


@pytest.mark.parametrize("ptm_type", sorted(ubiquitin_sumo.SUPPORTED_PTM_TYPES))
def test_caso_valido_no_lanza(ptm_type):
    pose = _FakePose(residues={24: "K"})
    ubiquitin_sumo._validate_target(pose, 24, ptm_type)  # no debe lanzar


def test_pdbs_de_referencia_empaquetados_existen():
    """Regresion real: si alguno de estos 2 PDBs (documentados extensamente en el docstring
    del modulo, ver procedencia RCSB 1UBQ/1A5R) se borra o renombra por error, conjugate()
    fallaria en produccion con un FileNotFoundError opaco de PyRosetta en vez de un error
    claro -- mejor detectarlo aqui."""
    for ptm_type, path in ubiquitin_sumo.CONJUGATE_PDB_BY_TYPE.items():
        assert path.is_file(), f"{ptm_type}: falta el PDB de referencia empaquetado '{path}'"


def test_supported_types_coincide_con_conjugate_pdb_by_type():
    assert ubiquitin_sumo.SUPPORTED_PTM_TYPES == frozenset(ubiquitin_sumo.CONJUGATE_PDB_BY_TYPE)


def test_metric_names_tiene_5_angulos():
    # _recompute_metrics hace zip(METRIC_NAMES, torsions) con exactamente 5 tuplas fijas --
    # si METRIC_NAMES cambia de longitud, el zip trunca en silencio en vez de fallar.
    assert len(ubiquitin_sumo.METRIC_NAMES) == 5


def test_lyx_params_path_lanza_si_no_existe(monkeypatch, tmp_path):
    fake_pyrosetta = types.ModuleType("pyrosetta")
    fake_pyrosetta.__file__ = str(tmp_path / "pyrosetta" / "__init__.py")
    monkeypatch.setitem(sys.modules, "pyrosetta", fake_pyrosetta)

    with pytest.raises(FileNotFoundError):
        ubiquitin_sumo._lyx_params_path()


def test_lyx_params_path_resuelve_relativo_al_paquete(monkeypatch, tmp_path):
    lyx_path = (
        tmp_path / "pyrosetta_pkg" / "database" / "chemical" / "residue_type_sets"
        / "fa_standard" / "residue_types" / "sidechain_conjugation" / "LYX.params"
    )
    lyx_path.parent.mkdir(parents=True)
    lyx_path.write_text("# fake LYX params\n")

    fake_pyrosetta = types.ModuleType("pyrosetta")
    fake_pyrosetta.__file__ = str(tmp_path / "pyrosetta_pkg" / "__init__.py")
    monkeypatch.setitem(sys.modules, "pyrosetta", fake_pyrosetta)

    assert ubiquitin_sumo._lyx_params_path() == lyx_path
