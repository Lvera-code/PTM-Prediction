"""Tests dedicados de src/structural/ddg_estimate.py (punto 6 del plan de robustez
post-demo-prep) -- la logica de "nstruct trayectorias independientes, se usa la de menor
energia por estado" (practica estandar de los protocolos ddG de Rosetta, ver docstring del
modulo) no tenia ninguna cobertura propia: se ejercitaba solo indirectamente via
test_fase_a_dispatch.py, que mockea ``estimate_ddg`` entero sin probar su propia logica
interna. ``load_pose``/``apply_ptm_patch``/``relax_neighborhood`` se monkeypatchean
directamente en el namespace de ddg_estimate (donde quedaron ligados por el ``from ... import``
del modulo); ``pyrosetta.create_score_function`` (import local dentro de estimate_ddg) se
stubea via sys.modules, mismo criterio que los otros tests de src/structural/.
"""

import sys
import types
from unittest.mock import MagicMock

import pytest

from src.structural import ddg_estimate as ddg_estimate_module


class _FakePose:
    """Distingue WT de la version parcheada solo para que las aserciones puedan verificar
    que apply_ptm_patch se aplico exactamente donde correspondia."""

    def __init__(self):
        self.patched = False


@pytest.fixture
def fake_pyrosetta_scorefxn(monkeypatch):
    """Stubea pyrosetta.create_score_function para que devuelva un scorefxn de mentira cuyo
    resultado por llamada se controla desde el test (via side_effect)."""
    fake_pyrosetta = types.ModuleType("pyrosetta")
    scorefxn = MagicMock()
    fake_pyrosetta.create_score_function = MagicMock(return_value=scorefxn)
    monkeypatch.setitem(sys.modules, "pyrosetta", fake_pyrosetta)
    return scorefxn


def test_best_of_n_relax_devuelve_el_minimo_y_todos_los_scores(monkeypatch):
    poses_creados = []

    def _fake_load_pose(pdb_path):
        pose = _FakePose()
        poses_creados.append(pose)
        return pose

    monkeypatch.setattr(ddg_estimate_module, "load_pose", _fake_load_pose)
    monkeypatch.setattr(ddg_estimate_module, "relax_neighborhood", lambda *a, **k: None)

    scorefxn = MagicMock(side_effect=[10.0, 7.5, 9.0])

    min_score, scores = ddg_estimate_module._best_of_n_relax(
        "fake.pdb", 24, None, scorefxn, nstruct=3, radius=6.0, max_iter=200,
    )

    assert scores == [10.0, 7.5, 9.0]
    assert min_score == 7.5
    assert len(poses_creados) == 3  # una pose fresca por trayectoria, no reusada


def test_best_of_n_relax_aplica_el_parche_solo_si_ptm_type_no_es_none(monkeypatch):
    applied_to = []
    monkeypatch.setattr(ddg_estimate_module, "load_pose", lambda pdb_path: _FakePose())
    monkeypatch.setattr(ddg_estimate_module, "apply_ptm_patch", lambda pose, pos, t: applied_to.append(t))
    monkeypatch.setattr(ddg_estimate_module, "relax_neighborhood", lambda *a, **k: None)
    scorefxn = MagicMock(return_value=5.0)

    ddg_estimate_module._best_of_n_relax("fake.pdb", 24, None, scorefxn, nstruct=2, radius=6.0, max_iter=200)
    assert applied_to == []  # WT: sin parche

    ddg_estimate_module._best_of_n_relax(
        "fake.pdb", 24, "acetylation", scorefxn, nstruct=2, radius=6.0, max_iter=200
    )
    assert applied_to == ["acetylation", "acetylation"]  # parcheado: una vez por trayectoria


def test_estimate_ddg_usa_el_minimo_de_cada_estado(monkeypatch, fake_pyrosetta_scorefxn):
    monkeypatch.setattr(ddg_estimate_module, "load_pose", lambda pdb_path: _FakePose())
    monkeypatch.setattr(ddg_estimate_module, "apply_ptm_patch", lambda pose, pos, t: pose)
    monkeypatch.setattr(ddg_estimate_module, "relax_neighborhood", lambda *a, **k: None)

    # 3 trayectorias WT, luego 3 trayectorias parcheadas (mismo orden que estimate_ddg llama
    # _best_of_n_relax: WT primero, mutante despues).
    fake_pyrosetta_scorefxn.side_effect = [10.0, 11.0, 10.5, 8.0, 9.0, 8.5]

    ddg, wt_score, mut_score, wt_scores, mut_scores = ddg_estimate_module.estimate_ddg(
        "fake.pdb", 24, "acetylation", radius=6.0, max_iter=200, nstruct=3,
    )

    assert wt_scores == [10.0, 11.0, 10.5]
    assert mut_scores == [8.0, 9.0, 8.5]
    assert wt_score == 10.0
    assert mut_score == 8.0
    assert ddg == pytest.approx(8.0 - 10.0)


def test_estimate_ddg_nstruct_1_es_una_sola_trayectoria_por_estado(monkeypatch, fake_pyrosetta_scorefxn):
    monkeypatch.setattr(ddg_estimate_module, "load_pose", lambda pdb_path: _FakePose())
    monkeypatch.setattr(ddg_estimate_module, "apply_ptm_patch", lambda pose, pos, t: pose)
    monkeypatch.setattr(ddg_estimate_module, "relax_neighborhood", lambda *a, **k: None)
    fake_pyrosetta_scorefxn.side_effect = [12.0, 9.0]

    ddg, wt_score, mut_score, wt_scores, mut_scores = ddg_estimate_module.estimate_ddg(
        "fake.pdb", 24, "acetylation", nstruct=1,
    )

    assert wt_scores == [12.0]
    assert mut_scores == [9.0]
    assert ddg == pytest.approx(9.0 - 12.0)
