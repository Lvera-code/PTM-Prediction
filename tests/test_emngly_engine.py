"""Tests de emngly_engine.get_emngly_predictions (src/engines/emngly_engine.py).

Mockea ``subprocess.run`` (nunca invoca fair-esm/torch/MIF/scikit-learn
reales, no instalados en el venv principal -- viven en el venv dedicado
``Settings.EMNGLY_PYTHON_BIN``). Motor REAL de consenso (a diferencia de
MeToken, puramente informativo, ver ``ptm_annotation.py``), pero SIGUE
degradando de forma no-fatal ante cualquier fallo -- estos tests verifican
esa propiedad para cada modo de fallo real (venv ausente, runner ausente,
clon/MIF ausente, pesos MIF/ESM-1b/SVM ausentes, exit code != 0, timeout,
salida malformada).
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.emngly_engine import OUTPUT_COLUMNS, get_emngly_predictions


def _install_fake_resources(tmp_path, monkeypatch, missing=None):
    """Crea todos los recursos que ``_validate_installation`` espera; ``missing`` omite uno."""
    missing = missing or set()

    python_bin = tmp_path / "python3"
    if "python_bin" not in missing:
        python_bin.write_text("#!/bin/sh\n")
        python_bin.chmod(0o755)
    monkeypatch.setattr(Settings, "EMNGLY_PYTHON_BIN", str(python_bin))

    runner_script = tmp_path / "_emngly_runner.py"
    if "runner_script" not in missing:
        runner_script.write_text("# fake\n")
    monkeypatch.setattr(Settings, "EMNGLY_RUNNER_SCRIPT", runner_script)

    emngly_home = tmp_path / "EMNgly"
    mif_init = emngly_home / "model" / "MIF" / "__init__.py"
    if "emngly_home" not in missing:
        mif_init.parent.mkdir(parents=True)
        mif_init.write_text("")
    monkeypatch.setattr(Settings, "EMNGLY_HOME", emngly_home)

    mif_weights = emngly_home / "model" / "MIF" / "weights" / "mif.pt"
    if "mif_weights" not in missing:
        mif_weights.parent.mkdir(parents=True, exist_ok=True)
        mif_weights.write_text("fake weights")
    monkeypatch.setattr(Settings, "EMNGLY_MIF_WEIGHTS", mif_weights)

    esm_checkpoint = tmp_path / "esm1b_t33_650M_UR50S.pt"
    if "esm_checkpoint" not in missing:
        esm_checkpoint.write_text("fake checkpoint")
    monkeypatch.setattr(Settings, "EMNGLY_ESM_CHECKPOINT", esm_checkpoint)

    svm_checkpoint = tmp_path / "N-GlyDE.pickle"
    if "svm_checkpoint" not in missing:
        svm_checkpoint.write_text("fake pickle")
    monkeypatch.setattr(Settings, "EMNGLY_SVM_CHECKPOINT", svm_checkpoint)

    monkeypatch.setattr(Settings, "EMNGLY_CACHE_DIR", tmp_path / "cache")


def _common_args(tmp_path):
    return dict(
        accession="ACC1", sequence="N" * 30, pdb_path=tmp_path / "fake.pdb",
        position_mapping_path=tmp_path / "fake_mapping.csv",
    )


def test_positions_vacias_devuelve_dict_vacio_sin_invocar_subprocess(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = get_emngly_predictions(positions=[], result_dir=tmp_path / "out", **_common_args(tmp_path))

    assert result == {}
    assert called == []


@pytest.mark.parametrize(
    "missing", ["python_bin", "runner_script", "emngly_home", "mif_weights", "esm_checkpoint", "svm_checkpoint"]
)
def test_recurso_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path, missing):
    _install_fake_resources(tmp_path, monkeypatch, missing={missing})

    result = get_emngly_predictions(positions=[5], result_dir=tmp_path / "out", **_common_args(tmp_path))

    assert result == {}


def test_corrida_exitosa_devuelve_dict_por_posicion(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    out_dir = tmp_path / "out"

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame([{"position": 5, "probability": 0.87}], columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_emngly_predictions(positions=[5], result_dir=out_dir, **_common_args(tmp_path))

    assert result == {5: {"emngly_probability": 0.87}}


def test_exit_code_no_cero_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(1, cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_emngly_predictions(positions=[5], result_dir=tmp_path / "out", **_common_args(tmp_path))

    assert result == {}


def test_timeout_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd, 1800)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_emngly_predictions(positions=[5], result_dir=tmp_path / "out", **_common_args(tmp_path))

    assert result == {}


def test_csv_no_generado_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0))

    result = get_emngly_predictions(positions=[5], result_dir=tmp_path / "out", **_common_args(tmp_path))

    assert result == {}


def test_csv_columnas_incompletas_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    out_dir = tmp_path / "out"

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame([{"position": 5}]).to_csv(out_csv, index=False)  # falta 'probability'
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_emngly_predictions(positions=[5], result_dir=out_dir, **_common_args(tmp_path))

    assert result == {}
