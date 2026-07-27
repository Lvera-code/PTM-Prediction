"""Tests de DeepPTMPredEngine (src/engines/deepptmpred_engine.py).

Mockea ``subprocess.run`` (nunca invoca PyRosetta/TensorFlow/fair-esm
reales, no instalados en esta maquina). El mock escribe un CSV sintetico en
la ruta ``--out-csv`` esperada, igual que haria el runner real.
"""

import subprocess

import pandas as pd
import pytest

from src.config.settings import Settings
from src.engines.deepptmpred_engine import DeepPTMPredEngine, OUTPUT_COLUMNS
from src.utils.exceptions import DeepPTMPredExecutionError
from src.utils.structure_parser import StructureRecord


def _make_record(tmp_path, accession="ACC1", sequence="MKTAYIAKQRQ"):
    pdb_path = tmp_path / f"{accession}_chain_A.pdb"
    pdb_path.write_text("HEADER  fake\n")
    return StructureRecord(
        accession=accession,
        pdb_path=pdb_path,
        chain_pdb_path=pdb_path,
        fasta_path=tmp_path / f"{accession}_derived.fasta",
        chain_id="A",
        sequence=sequence,
        position_mapping=pd.DataFrame(),
    )


def _make_engine(tmp_path, ptm_types=("phosphorylation", "acetylation")):
    train_ptm_dir = tmp_path / "DeepPTMPred" / "pred" / "train_PTM"
    train_ptm_dir.mkdir(parents=True, exist_ok=True)
    (train_ptm_dir / "predict.py").write_text("# fake\n")

    esm_checkpoint = tmp_path / "esm" / "esm2_t33_650M_UR50D.pt"
    esm_checkpoint.parent.mkdir(parents=True, exist_ok=True)
    esm_checkpoint.write_text("fake checkpoint")

    return DeepPTMPredEngine(
        train_ptm_dir=train_ptm_dir,
        runner_script=tmp_path / "_deepptmpred_runner.py",
        python_bin="python3",
        esm_checkpoint=esm_checkpoint,
        custom_esm_dir=tmp_path / "custom_esm",
        ptm_types=ptm_types,
    )


def _mock_run_writing_csv(row_factory):
    """row_factory(ptm_type) -> list de dicts para ese tipo de PTM."""

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        ptm_type = cmd[cmd.index("--ptm-type") + 1]
        rows = row_factory(ptm_type)
        pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _fake_run


def test_run_concatena_todos_los_tipos_de_ptm(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation", "acetylation"))
    record = _make_record(tmp_path)

    def _rows(ptm_type):
        return [["ACC1", 3, "T", 0.8, ptm_type]]

    monkeypatch.setattr(subprocess, "run", _mock_run_writing_csv(_rows))

    results = engine.run([record], output_dir=tmp_path / "out")

    assert len(results) == 1
    df = results[0]
    assert list(df.columns) == OUTPUT_COLUMNS
    assert sorted(df["ptm_type"].tolist()) == ["acetylation", "phosphorylation"]
    assert len(df) == 2


def test_cmd_incluye_un_tipo_de_ptm_por_invocacion(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record = _make_record(tmp_path)
    captured = []

    def _fake_run(cmd, **kwargs):
        captured.append(cmd)
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame([["ACC1", 3, "T", 0.8, "phosphorylation"]], columns=OUTPUT_COLUMNS).to_csv(
            out_csv, index=False
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)
    engine.run([record], output_dir=tmp_path / "out")

    assert len(captured) == 1
    cmd = captured[0]
    assert cmd[cmd.index("--ptm-type") + 1] == "phosphorylation"
    assert cmd[cmd.index("--protein-id") + 1] == "ACC1"
    assert cmd[cmd.index("--sequence") + 1] == "MKTAYIAKQRQ"


def test_multiples_accessions_invocan_el_runner_por_cada_uno(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record1 = _make_record(tmp_path, accession="ACC1")
    record2 = _make_record(tmp_path, accession="ACC2")

    def _rows(ptm_type):
        return []

    monkeypatch.setattr(subprocess, "run", _mock_run_writing_csv(_rows))

    results = engine.run([record1, record2], output_dir=tmp_path / "out")

    assert len(results) == 2


def test_exit_code_distinto_de_cero_propaga_deepptmpred_execution_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record = _make_record(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")


def test_timeout_propaga_deepptmpred_execution_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record = _make_record(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")


def test_output_sin_columnas_esperadas_lanza_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record = _make_record(tmp_path)

    def _fake_run(cmd, **kwargs):
        out_csv = cmd[cmd.index("--out-csv") + 1]
        pd.DataFrame({"foo": [1]}).to_csv(out_csv, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")


def test_output_file_ausente_lanza_error(monkeypatch, tmp_path):
    engine = _make_engine(tmp_path, ptm_types=("phosphorylation",))
    record = _make_record(tmp_path)

    def _fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")


def test_repo_ausente_lanza_error_accionable(tmp_path):
    engine = DeepPTMPredEngine(
        train_ptm_dir=tmp_path / "no_existe",
        runner_script=tmp_path / "_runner.py",
        python_bin="python3",
        esm_checkpoint=tmp_path / "esm.pt",
        custom_esm_dir=tmp_path / "custom_esm",
        ptm_types=("phosphorylation",),
    )
    record = _make_record(tmp_path)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")


def test_checkpoint_esm_ausente_lanza_error_accionable(tmp_path):
    train_ptm_dir = tmp_path / "DeepPTMPred" / "pred" / "train_PTM"
    train_ptm_dir.mkdir(parents=True)
    (train_ptm_dir / "predict.py").write_text("# fake\n")

    engine = DeepPTMPredEngine(
        train_ptm_dir=train_ptm_dir,
        runner_script=tmp_path / "_runner.py",
        python_bin="python3",
        esm_checkpoint=tmp_path / "no_existe.pt",
        custom_esm_dir=tmp_path / "custom_esm",
        ptm_types=("phosphorylation",),
    )
    record = _make_record(tmp_path)

    with pytest.raises(DeepPTMPredExecutionError):
        engine.run([record], output_dir=tmp_path / "out")
