"""Tests de stackglyembed_engine.get_nglyco_corroboration (src/engines/stackglyembed_engine.py).

Mockea ``subprocess.run`` (nunca invoca torch/tensorflow/transformers/
ProteinBERT reales, no instalados en el venv principal -- viven en el venv
dedicado del proyecto hermano). Corroboracion PURAMENTE informativa (mismo
patron que MeToken/GlyGen): cualquier fallo debe degradar a un dict vacio,
NUNCA lanzar -- estos tests verifican explicitamente esa propiedad para cada
modo de fallo real (venv ausente, runner ausente, pickles ausentes, exit
code != 0, timeout, salida malformada, archivo de salida ausente).
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.stackglyembed_engine import OUTPUT_COLUMNS, get_nglyco_corroboration


def _install_fake_resources(tmp_path, monkeypatch, with_pickle=True):
    python_bin = tmp_path / "python3"
    python_bin.write_text("#!/bin/sh\n")
    python_bin.chmod(0o755)
    monkeypatch.setattr(Settings, "STACKGLYEMBED_PYTHON_BIN", str(python_bin))

    runner_script = tmp_path / "_stackglyembed_runner.py"
    runner_script.write_text("# fake\n")
    monkeypatch.setattr(Settings, "STACKGLYEMBED_RUNNER_SCRIPT", runner_script)

    models_dir = tmp_path / "prediction"
    pickle_dir = models_dir / "base_layer_pickle_files"
    pickle_dir.mkdir(parents=True)
    if with_pickle:
        (pickle_dir / "SVM_meta_layer.sav").write_text("fake pickle")
    monkeypatch.setattr(Settings, "STACKGLYEMBED_MODELS_DIR", str(models_dir))
    return python_bin, runner_script, models_dir


def test_positions_vacias_devuelve_dict_vacio_sin_invocar_subprocess(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = get_nglyco_corroboration("AAAA", [], result_dir=tmp_path / "out")

    assert result == {}
    assert called == []


def test_python_bin_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "STACKGLYEMBED_PYTHON_BIN", str(tmp_path / "no_existe"))

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_runner_script_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    python_bin = tmp_path / "python3"
    python_bin.write_text("#!/bin/sh\n")
    python_bin.chmod(0o755)
    monkeypatch.setattr(Settings, "STACKGLYEMBED_PYTHON_BIN", str(python_bin))
    monkeypatch.setattr(Settings, "STACKGLYEMBED_RUNNER_SCRIPT", tmp_path / "no_existe.py")

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_pickles_ausentes_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch, with_pickle=False)

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_corrida_exitosa_devuelve_dict_por_posicion(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        rows = [
            {"position": 5, "stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.91},
            {"position": 20, "stackglyembed_veredicto": "No glicosilado", "stackglyembed_score": 0.12},
        ]
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_nglyco_corroboration("N" * 30, [5, 20], result_dir=tmp_path / "out")

    assert result == {
        5: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.91},
        20: {"stackglyembed_veredicto": "No glicosilado", "stackglyembed_score": 0.12},
    }


def test_cmd_incluye_secuencia_y_posiciones(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    captured = []

    def _fake_run(cmd, **kwargs):
        captured.append(cmd)
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    get_nglyco_corroboration("NAT" + "A" * 20, [1, 10], result_dir=tmp_path / "out")

    cmd = captured[0]
    assert cmd[cmd.index("--sequence") + 1] == "NAT" + "A" * 20
    positions_idx = cmd.index("--positions")
    assert cmd[positions_idx + 1 : cmd.index("--models-dir")] == ["1", "10"]


def test_exit_code_distinto_de_cero_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_timeout_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_output_sin_columnas_esperadas_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame({"foo": [1]}).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_output_file_ausente_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""))

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}


def test_python_bin_no_invocable_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise FileNotFoundError("no such file")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_nglyco_corroboration("N" * 20, [5], result_dir=tmp_path / "out")

    assert result == {}
