"""Tests de kinase_library_engine.get_kinase_corroboration (src/engines/kinase_library_engine.py).

Mockea ``subprocess.run`` (nunca invoca ``kinase_library`` real, que vive en
el entorno conda dedicado ``kinase_library``, no en el venv principal --
``numpy``/``pandas`` que fija el paquete son incompatibles con las versiones
fijadas de este venv). Corroboracion PURAMENTE informativa (mismo patron que
MeToken): cualquier fallo debe degradar a un dict vacio, NUNCA
lanzar -- estos tests verifican explicitamente esa propiedad para cada modo
de fallo real (entorno ausente, runner ausente, exit code != 0, timeout,
salida malformada, archivo de salida ausente).
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.kinase_library_engine import OUTPUT_COLUMNS, get_kinase_corroboration


def _install_fake_resources(tmp_path, monkeypatch):
    python_bin = tmp_path / "python3"
    python_bin.write_text("#!/bin/sh\n")
    python_bin.chmod(0o755)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_PYTHON_BIN", str(python_bin))

    runner_script = tmp_path / "_kinase_library_runner.py"
    runner_script.write_text("# fake\n")
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_RUNNER_SCRIPT", runner_script)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", True)
    return python_bin, runner_script


def test_positions_vacias_devuelve_dict_vacio_sin_invocar_subprocess(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = get_kinase_corroboration("AAAA", [], result_dir=tmp_path / "out")

    assert result == {}
    assert called == []


def test_deshabilitado_via_settings_devuelve_dict_vacio_sin_invocar_subprocess(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}
    assert called == []


def test_python_bin_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_PYTHON_BIN", str(tmp_path / "no_existe"))
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", True)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_runner_script_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    python_bin = tmp_path / "python3"
    python_bin.write_text("#!/bin/sh\n")
    python_bin.chmod(0o755)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_PYTHON_BIN", str(python_bin))
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_RUNNER_SCRIPT", tmp_path / "no_existe.py")
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", True)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_corrida_exitosa_devuelve_dict_por_posicion(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        rows = [
            {
                "position": 10, "kinase_library_top_kinase": "ATM", "kinase_library_top_family": "PIKK",
                "kinase_library_percentile": 99.83, "kinase_library_top3_kinases": "ATM,SMG1,ATR",
            },
            {
                "position": 25, "kinase_library_top_kinase": "PKACA", "kinase_library_top_family": "AGC",
                "kinase_library_percentile": 85.4, "kinase_library_top3_kinases": "PKACA,PKACB,RSK2",
            },
        ]
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_kinase_corroboration("S" * 30, [10, 25], result_dir=tmp_path / "out")

    assert result == {
        10: {
            "kinase_library_top_kinase": "ATM", "kinase_library_top_family": "PIKK",
            "kinase_library_percentile": 99.83, "kinase_library_top3_kinases": "ATM,SMG1,ATR",
        },
        25: {
            "kinase_library_top_kinase": "PKACA", "kinase_library_top_family": "AGC",
            "kinase_library_percentile": 85.4, "kinase_library_top3_kinases": "PKACA,PKACB,RSK2",
        },
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
    get_kinase_corroboration("SAT" + "A" * 20, [1, 10], result_dir=tmp_path / "out")

    cmd = captured[0]
    assert cmd[cmd.index("--sequence") + 1] == "SAT" + "A" * 20
    positions_idx = cmd.index("--positions")
    assert cmd[positions_idx + 1 : cmd.index("--out-csv")] == ["1", "10"]


def test_exit_code_distinto_de_cero_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_timeout_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_output_sin_columnas_esperadas_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame({"foo": [1]}).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_output_file_ausente_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""))

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}


def test_python_bin_no_invocable_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_resources(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise FileNotFoundError("no such file")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_kinase_corroboration("S" * 20, [10], result_dir=tmp_path / "out")

    assert result == {}
