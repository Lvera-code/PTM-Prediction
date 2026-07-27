"""Tests de DeepMVPEngine (src/engines/deepmvp_engine.py).

Mockea ``subprocess.run`` en vez de invocar el repo/pesos reales de DeepMVP
(no descargados en CI). El mock escribe un ``site_prediction.tsv`` sintetico
en la carpeta de resultado esperada, igual que haria el subproceso real.
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.deepmvp_engine import DeepMVPEngine, OUTPUT_COLUMNS
from src.utils.exceptions import DeepMVPExecutionError


def _make_engine(tmp_path, deepmvp_home=None, model_dir=None):
    home = deepmvp_home or (tmp_path / "DeepMVP")
    home.mkdir(parents=True, exist_ok=True)
    (home / "DeepMVP.py").write_text("# fake\n")

    models = model_dir or (home / "models")
    models.mkdir(parents=True, exist_ok=True)
    (models / "acetylation_k").mkdir()

    return DeepMVPEngine(deepmvp_home=home, python_bin="python3", model_dir=models)


def _write_fasta(tmp_path, name="ACC1.fasta"):
    p = tmp_path / name
    p.write_text(">ACC1\nMKTAYIAKQRQ\n")
    return p


def _mock_run_writing_tsv(rows):
    """Devuelve un fake de subprocess.run que escribe site_prediction.tsv en -o."""

    def _fake_run(cmd, **kwargs):
        out_dir = cmd[cmd.index("-o") + 1]
        df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
        df.to_csv(f"{out_dir}/site_prediction.tsv", sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _fake_run


def test_run_devuelve_dataframe_con_columnas_esperadas(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)
    rows = [["ACC1", "K", 17, "ASGSCQGCEEDEETLKKLIVRLNNVQEGKQI", 0.99, 0.002, "acetylation_k"]]
    monkeypatch.setattr(subprocess, "run", _mock_run_writing_tsv(rows))

    results = engine.run([str(fasta)], output_dir=tmp_path / "out")

    assert len(results) == 1
    assert list(results[0].columns) == OUTPUT_COLUMNS
    assert results[0].iloc[0]["protein"] == "ACC1"
    assert results[0].iloc[0]["ptm"] == "acetylation_k"


def test_cmd_incluye_tarea_2_y_paths_absolutos(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)
    captured = {}

    def _fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        out_dir = cmd[cmd.index("-o") + 1]
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(f"{out_dir}/site_prediction.tsv", sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    engine.run([str(fasta)], output_dir=tmp_path / "out")

    cmd = captured["cmd"]
    assert cmd[2] == "predict"
    assert "-t" in cmd and cmd[cmd.index("-t") + 1] == "2"
    assert "-d" in cmd and cmd[cmd.index("-d") + 1] == str(fasta.resolve())


def test_exit_code_distinto_de_cero_propaga_deepmvp_execution_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_timeout_propaga_deepmvp_execution_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_output_sin_columnas_esperadas_lanza_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)

    def _fake_run(cmd, **kwargs):
        out_dir = cmd[cmd.index("-o") + 1]
        pd.DataFrame({"foo": [1]}).to_csv(f"{out_dir}/site_prediction.tsv", sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_output_file_ausente_lanza_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path)
    fasta = _write_fasta(tmp_path)

    def _fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_repo_ausente_lanza_error_accionable(tmp_path):
    engine = DeepMVPEngine(
        deepmvp_home=tmp_path / "no_existe", python_bin="python3", model_dir=tmp_path / "models"
    )
    fasta = _write_fasta(tmp_path)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_modelos_ausentes_lanza_error_accionable(tmp_path):
    home = tmp_path / "DeepMVP"
    home.mkdir()
    (home / "DeepMVP.py").write_text("# fake\n")
    empty_models = tmp_path / "models_empty"
    empty_models.mkdir()

    engine = DeepMVPEngine(deepmvp_home=home, python_bin="python3", model_dir=empty_models)
    fasta = _write_fasta(tmp_path)

    with pytest.raises(DeepMVPExecutionError):
        engine.run([str(fasta)], output_dir=tmp_path / "out")


def test_fasta_inexistente_lanza_file_not_found(tmp_path):
    engine = _make_engine(tmp_path)

    with pytest.raises(FileNotFoundError):
        engine.run([str(tmp_path / "no_existe.fasta")], output_dir=tmp_path / "out")
