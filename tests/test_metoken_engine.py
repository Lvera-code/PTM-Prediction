"""Tests de metoken_engine.get_type_corroboration (src/engines/metoken_engine.py).

Mockea ``subprocess.run`` (nunca invoca torch/torch_scatter/transformers
reales, no instalados en el venv principal). Corroboracion PURAMENTE
informativa: cualquier fallo debe degradar a un dict vacio, NUNCA lanzar --
estos tests verifican explicitamente esa
propiedad para cada modo de fallo real (repo ausente, checkpoint ausente,
exit code != 0, timeout, salida malformada, archivo de salida ausente).
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.metoken_engine import OUTPUT_COLUMNS, get_type_corroboration


def _install_fake_repo(tmp_path, monkeypatch, with_checkpoint=True):
    repo_dir = tmp_path / "MeToken"
    repo_dir.mkdir()
    (repo_dir / "inference.py").write_text("# fake\n")
    monkeypatch.setattr(Settings, "METOKEN_HOME", repo_dir)

    checkpoint = tmp_path / "checkpoint.ckpt"
    if with_checkpoint:
        checkpoint.write_text("fake checkpoint")
    monkeypatch.setattr(Settings, "METOKEN_CHECKPOINT", checkpoint)
    monkeypatch.setattr(Settings, "METOKEN_RUNNER_SCRIPT", tmp_path / "_metoken_runner.py")
    monkeypatch.setattr(Settings, "METOKEN_PYTHON_BIN", "python3")
    return repo_dir, checkpoint


def test_positions_vacias_devuelve_dict_vacio_sin_invocar_subprocess(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = get_type_corroboration(tmp_path / "fake.pdb", [], result_dir=tmp_path / "out")

    assert result == {}
    assert called == []


def test_repo_no_instalado_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "METOKEN_HOME", tmp_path / "no_existe")

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_checkpoint_ausente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch, with_checkpoint=False)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_corrida_exitosa_devuelve_dict_por_posicion(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        rows = [
            {"position": 10, "metoken_type": "Phosphorylation", "metoken_probability": 0.87},
            {"position": 20, "metoken_type": "Ubiquitination", "metoken_probability": 0.61},
        ]
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10, 20], chain_id="A", result_dir=tmp_path / "out")

    assert result == {
        10: {"metoken_type": "Phosphorylation", "metoken_probability": 0.87},
        20: {"metoken_type": "Ubiquitination", "metoken_probability": 0.61},
    }


def test_cmd_incluye_posiciones_y_chain_id(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)
    captured = []

    def _fake_run(cmd, **kwargs):
        captured.append(cmd)
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    get_type_corroboration(tmp_path / "fake.pdb", [5, 8], chain_id="B", result_dir=tmp_path / "out")

    cmd = captured[0]
    assert cmd[cmd.index("--chain-id") + 1] == "B"
    positions_idx = cmd.index("--positions")
    # --positions consume argumentos hasta el siguiente flag (--out-csv)
    assert cmd[positions_idx + 1 : cmd.index("--out-csv")] == ["5", "8"]


def test_exit_code_distinto_de_cero_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_timeout_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_output_sin_columnas_esperadas_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame({"foo": [1]}).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_output_file_ausente_devuelve_dict_vacio(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    monkeypatch.setattr(subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""))

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}


def test_python_bin_inexistente_devuelve_dict_vacio_sin_lanzar(monkeypatch, tmp_path):
    _install_fake_repo(tmp_path, monkeypatch)

    def _fake_run(cmd, **kwargs):
        raise FileNotFoundError("no such file")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = get_type_corroboration(tmp_path / "fake.pdb", [10], result_dir=tmp_path / "out")

    assert result == {}
