"""Tests dedicados de src/structural/fase_a_dispatch.py (punto 6 del plan de robustez
post-demo-prep, STATUS.md/diario del vault 2026-08-03) -- el enrutador de un sitio PTM a su
modulo de modelado estructural (clase 1/2/3) no tenia ninguna cobertura propia hasta ahora,
solo se ejercitaba indirectamente via tests/test_fase_a_engine.py (que mockea FaseAEngine
entero, sin llegar nunca a este modulo).

``fase_a_dispatch.py`` en si NO importa pyrosetta a nivel de modulo (solo los 3 submodulos
que reusa, que tampoco lo hacen -- ver sus propios docstrings), asi que es importable sin
pyrosetta instalado. Las funciones _run_classN SI hacen ``from pyrosetta import
pose_from_pdb`` de forma incondicional (clase 2 y 3) -- se stubea un modulo 'pyrosetta'
minimo en sys.modules para esas, en vez de mockear a nivel de subprocess (mismo criterio
de mocking que el resto del proyecto: mockear en el limite mas bajo posible sin necesitar
el binario/libreria real instalada).
"""

import statistics
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from src.structural import fase_a_dispatch, pyrosetta_glycan_patch, pyrosetta_ptm_patch, ubiquitin_sumo
from src.structural import ddg_estimate as ddg_estimate_module


class _FakePose:
    def __init__(self):
        self.dumped_to = None

    def dump_pdb(self, path):
        self.dumped_to = path


@pytest.fixture
def fake_pyrosetta_pose_from_pdb(monkeypatch):
    """Stubea ``from pyrosetta import pose_from_pdb`` para las clases 2/3 (import
    incondicional dentro de _run_class2/_run_class3 en fase_a_dispatch.py)."""
    fake_module = types.ModuleType("pyrosetta")
    fake_module.pose_from_pdb = MagicMock(return_value=_FakePose())
    monkeypatch.setitem(sys.modules, "pyrosetta", fake_module)
    return fake_module.pose_from_pdb


def test_tipo_sin_soporte_no_inicializa_nada(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(pyrosetta_ptm_patch, "init_pyrosetta", lambda *a, **k: calls.append("ptm"))
    monkeypatch.setattr(pyrosetta_glycan_patch, "init_pyrosetta", lambda *a, **k: calls.append("glycan"))
    monkeypatch.setattr(ubiquitin_sumo, "init_pyrosetta", lambda *a, **k: calls.append("ubq"))

    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 24, "malonylation", tmp_path / "out.pdb"
    )

    assert result["estado"] == "sin_soporte_fase_a"
    assert calls == []  # ningun modulo de PyRosetta se toco para un tipo sin soporte


def test_clase1_rutea_a_ptm_patch_y_propaga_ddg_std(tmp_path, monkeypatch):
    fake_pose = _FakePose()
    monkeypatch.setattr(pyrosetta_ptm_patch, "init_pyrosetta", lambda *a, **k: None)
    monkeypatch.setattr(pyrosetta_ptm_patch, "load_pose", lambda pdb_path: fake_pose)
    monkeypatch.setattr(pyrosetta_ptm_patch, "apply_ptm_patch", lambda pose, pos, t: pose)
    monkeypatch.setattr(pyrosetta_ptm_patch, "relax_neighborhood", lambda pose, pos, radius=6.0: pose)

    wt_scores = [10.0, 10.0, 12.0]
    mut_scores = [8.0, 9.0, 8.5]
    monkeypatch.setattr(
        ddg_estimate_module, "estimate_ddg",
        lambda pdb_path, position, ptm_type, radius, nstruct: (
            min(mut_scores) - min(wt_scores), min(wt_scores), min(mut_scores), wt_scores, mut_scores,
        ),
    )

    out_pdb = tmp_path / "out.pdb"
    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 24, "acetylation", out_pdb, radius=6.0, nstruct=3,
    )

    assert result["estado"] == "modelado"
    assert result["clase"] == "class1_patch_ddg"
    assert result["ddg"] == pytest.approx(min(mut_scores) - min(wt_scores))
    expected_std = (statistics.pstdev(wt_scores) ** 2 + statistics.pstdev(mut_scores) ** 2) ** 0.5
    assert result["ddg_std"] == pytest.approx(expected_std)
    assert result["output_pdb"] == str(out_pdb)
    assert fake_pose.dumped_to == str(out_pdb)


def test_clase2_rutea_a_glycan_patch_con_glygen_opcional(tmp_path, monkeypatch, fake_pyrosetta_pose_from_pdb):
    monkeypatch.setattr(pyrosetta_glycan_patch, "init_pyrosetta", lambda: None)
    monkeypatch.setattr(pyrosetta_glycan_patch, "attach_glycan", lambda pose, pos, t: pose)
    monkeypatch.setattr(pyrosetta_glycan_patch, "refine_glycan", lambda pose, rounds=1: pose)
    monkeypatch.setattr(
        pyrosetta_glycan_patch, "check_glygen_evidence",
        lambda accession, pos, t: f"evidencia para {accession}",
    )

    out_pdb = tmp_path / "out.pdb"
    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 484, "n_linked_glycosylation", out_pdb,
        uniprot_accession="P10636",
    )

    assert result["estado"] == "modelado"
    assert result["clase"] == "class2_glycan"
    assert result["glycan_tree"] == "N-glycan_core"
    assert result["glygen_evidencia"] == "evidencia para P10636"
    fake_pyrosetta_pose_from_pdb.assert_called_once_with(str(tmp_path / "in.pdb"))


def test_clase2_sin_accession_no_consulta_glygen(tmp_path, monkeypatch, fake_pyrosetta_pose_from_pdb):
    monkeypatch.setattr(pyrosetta_glycan_patch, "init_pyrosetta", lambda: None)
    monkeypatch.setattr(pyrosetta_glycan_patch, "attach_glycan", lambda pose, pos, t: pose)
    monkeypatch.setattr(pyrosetta_glycan_patch, "refine_glycan", lambda pose, rounds=1: pose)
    check_glygen = MagicMock()
    monkeypatch.setattr(pyrosetta_glycan_patch, "check_glygen_evidence", check_glygen)

    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 52, "o_linked_glycosylation", tmp_path / "out.pdb",
    )

    assert result["glygen_evidencia"] is None
    check_glygen.assert_not_called()


def test_clase3_rutea_a_ubiquitin_sumo(tmp_path, monkeypatch, fake_pyrosetta_pose_from_pdb):
    monkeypatch.setattr(ubiquitin_sumo, "init_pyrosetta", lambda refine_cycles, refine_repack_cycles: None)
    fake_metrics = {"lysine_chi3_CG-CD-CE-NZ": 12.3}
    monkeypatch.setattr(
        ubiquitin_sumo, "conjugate", lambda pose, pos, t: (pose, fake_metrics)
    )

    out_pdb = tmp_path / "out.pdb"
    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 24, "ubiquitination", out_pdb,
    )

    assert result["estado"] == "modelado"
    assert result["clase"] == "class3_conjugation"
    assert result["conjugation_metrics"] == fake_metrics


def test_excepcion_de_un_submodulo_se_traduce_a_estado_error_sin_propagar(tmp_path, monkeypatch):
    def _boom(*a, **k):
        raise pyrosetta_ptm_patch.ResidueMismatchError("posicion 24 no es K")

    monkeypatch.setattr(pyrosetta_ptm_patch, "init_pyrosetta", lambda *a, **k: None)
    monkeypatch.setattr(ddg_estimate_module, "estimate_ddg", _boom)

    result = fase_a_dispatch.run_fase_a_for_site(
        tmp_path / "in.pdb", 24, "acetylation", tmp_path / "out.pdb",
    )

    assert result["estado"] == "error"
    assert "ResidueMismatchError" in result["error"]
    assert "posicion 24 no es K" in result["error"]


def test_supported_types_coincide_con_settings():
    """Guardarail real que fase_a_dispatch.py ya valida en tiempo de import (ver su propio
    docstring) -- lo re-verificamos explicitamente aqui para que un test falle claro si algun
    dia SUPPORTED_PTM_TYPES y Settings.FASE_A_SUPPORTED_PTM_TYPES divergen, en vez de depender
    solo del RuntimeError de import (que no muestra en pytest cual de las dos listas cambio)."""
    from src.config.settings import Settings

    assert fase_a_dispatch.SUPPORTED_PTM_TYPES == frozenset(Settings.FASE_A_SUPPORTED_PTM_TYPES)
    assert fase_a_dispatch.CLASS1_TYPES | fase_a_dispatch.CLASS2_TYPES | fase_a_dispatch.CLASS3_TYPES == (
        fase_a_dispatch.SUPPORTED_PTM_TYPES
    )
