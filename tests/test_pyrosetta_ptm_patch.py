"""Tests dedicados de src/structural/pyrosetta_ptm_patch.py (punto 6 del plan de robustez
post-demo-prep). ``apply_ptm_patch`` hace imports incondicionales de pyrosetta.rosetta.*
ANTES de sus propias validaciones (Unsupported/ResidueMismatch) -- a diferencia de
check_glygen_evidence (glycan_patch), no se puede testear la logica de validacion sin que
esos imports resuelvan a algo. Se stubea un arbol minimo de modulos 'pyrosetta.rosetta.*' en
sys.modules (mismo criterio ya usado en test_fase_a_dispatch.py para pose_from_pdb, aqui un
nivel mas profundo del paquete).
"""

import sys
import types
from unittest.mock import MagicMock

import pytest

from src.structural import pyrosetta_ptm_patch


class _FakeResidue:
    def __init__(self, name1):
        self._name1 = name1

    def name1(self):
        return self._name1


class _FakePose:
    def __init__(self, residues=None):
        self._residues = residues or {}

    def residue(self, position):
        return _FakeResidue(self._residues.get(position, "A"))


class _FakeVariantType:
    PHOSPHORYLATION = "PHOSPHORYLATION"
    ACETYLATION = "ACETYLATION"
    HYDROXYLATION = "HYDROXYLATION"
    CARBOXYLATION = "CARBOXYLATION"
    METHYLATION = "METHYLATION"


class _FakeResidueTypeSet:
    def __init__(self):
        self.name_map = MagicMock(side_effect=lambda name: f"residue_type:{name}")


class _FakeChemicalManager:
    _instance = None

    def __init__(self):
        self.residue_type_set = MagicMock(return_value=_FakeResidueTypeSet())

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = _FakeChemicalManager()
        return cls._instance


@pytest.fixture
def fake_pyrosetta_chemical(monkeypatch):
    """Stubea 'pyrosetta.rosetta.core.chemical'/'pyrosetta.rosetta.core.pose' en sys.modules
    para que los 'from pyrosetta.rosetta.core.X import Y' incondicionales de apply_ptm_patch
    resuelvan, sin necesitar pyrosetta real instalado. Incluye tanto la ruta de
    add_variant_type_to_pose_residue (4 tipos) como la de ChemicalManager/
    replace_pose_residue_copying_existing_coordinates (hydroxylation, ver
    PTM_BASE_NAME_MAP -- causa raiz real documentada 2026-08-04 en el modulo)."""
    root = types.ModuleType("pyrosetta")
    rosetta = types.ModuleType("pyrosetta.rosetta")
    core = types.ModuleType("pyrosetta.rosetta.core")
    chemical = types.ModuleType("pyrosetta.rosetta.core.chemical")
    pose_mod = types.ModuleType("pyrosetta.rosetta.core.pose")

    chemical.VariantType = _FakeVariantType
    _FakeChemicalManager._instance = None  # aislar entre tests
    chemical.ChemicalManager = _FakeChemicalManager
    add_variant_mock = MagicMock()
    pose_mod.add_variant_type_to_pose_residue = add_variant_mock
    replace_residue_mock = MagicMock()
    pose_mod.replace_pose_residue_copying_existing_coordinates = replace_residue_mock

    root.rosetta = rosetta
    rosetta.core = core
    core.chemical = chemical
    core.pose = pose_mod

    for name, module in (
        ("pyrosetta", root),
        ("pyrosetta.rosetta", rosetta),
        ("pyrosetta.rosetta.core", core),
        ("pyrosetta.rosetta.core.chemical", chemical),
        ("pyrosetta.rosetta.core.pose", pose_mod),
    ):
        monkeypatch.setitem(sys.modules, name, module)

    return add_variant_mock


def test_tipo_no_soportado(fake_pyrosetta_chemical):
    pose = _FakePose(residues={5: "K"})
    with pytest.raises(pyrosetta_ptm_patch.UnsupportedPTMPatchError):
        pyrosetta_ptm_patch.apply_ptm_patch(pose, 5, "malonylation")
    fake_pyrosetta_chemical.assert_not_called()


def test_residuo_no_coincide_con_el_esperado(fake_pyrosetta_chemical):
    pose = _FakePose(residues={5: "A"})  # acetylation espera K
    with pytest.raises(pyrosetta_ptm_patch.ResidueMismatchError):
        pyrosetta_ptm_patch.apply_ptm_patch(pose, 5, "acetylation")
    fake_pyrosetta_chemical.assert_not_called()


def test_caso_valido_aplica_el_variant_correcto(fake_pyrosetta_chemical):
    pose = _FakePose(residues={24: "K"})

    result = pyrosetta_ptm_patch.apply_ptm_patch(pose, 24, "acetylation")

    assert result is pose
    fake_pyrosetta_chemical.assert_called_once_with(pose, _FakeVariantType.ACETYLATION, 24)


@pytest.mark.parametrize("ptm_type,expected_residue", [
    ("phosphorylation", "S"), ("acetylation", "K"), ("hydroxylation", "P"),
    ("gamma_carboxyglutamic_acid", "E"), ("lys_methylation", "K"),
])
def test_todos_los_tipos_soportados_pasan_con_su_residuo_esperado(
    fake_pyrosetta_chemical, ptm_type, expected_residue,
):
    pose = _FakePose(residues={10: expected_residue})
    pyrosetta_ptm_patch.apply_ptm_patch(pose, 10, ptm_type)  # no debe lanzar


def test_hydroxylation_reemplaza_el_residue_type_en_vez_de_anadir_variant(fake_pyrosetta_chemical):
    """Causa raiz real (2026-08-04, ver docstring del modulo): los patches de hidroxiprolina
    son los unicos de fa_standard con SET_BASE_NAME -> construyen un ResidueType BASE nuevo
    (HYP/0AZ), no una variante de PRO -- add_variant_type_to_pose_residue nunca puede
    resolverlo. apply_ptm_patch usa replace_pose_residue_copying_existing_coordinates en su
    lugar, vía PTM_BASE_NAME_MAP."""
    pose = _FakePose(residues={499: "P"})
    replace_residue_mock = sys.modules["pyrosetta.rosetta.core.pose"].replace_pose_residue_copying_existing_coordinates

    result = pyrosetta_ptm_patch.apply_ptm_patch(pose, 499, "hydroxylation")

    assert result is pose
    fake_pyrosetta_chemical.assert_not_called()  # NO usa add_variant_type_to_pose_residue
    replace_residue_mock.assert_called_once_with(pose, 499, "residue_type:HYP")


def test_variant_map_y_target_residue_cubren_los_mismos_tipos():
    assert set(pyrosetta_ptm_patch.PTM_VARIANT_MAP) == set(pyrosetta_ptm_patch.PTM_TARGET_RESIDUE)
    assert pyrosetta_ptm_patch.SUPPORTED_PTM_TYPES == frozenset(pyrosetta_ptm_patch.PTM_VARIANT_MAP)


def test_base_name_map_solo_cubre_hydroxylation():
    assert dict(pyrosetta_ptm_patch.PTM_BASE_NAME_MAP) == {"hydroxylation": "HYP"}
