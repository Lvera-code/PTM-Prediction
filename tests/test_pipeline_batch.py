"""Tests del modo batch (--input apunta a un directorio), punto 5 del plan de robustez
post-demo-prep (STATUS.md, diario del vault 2026-08-03). Mismo criterio de mocking que
test_pipeline_fase1.py: Fase 1/1.5 corren con logica real, los motores se mockean a nivel
de ``Engine.run``.
"""

from pathlib import Path

import pandas as pd

import pipeline
from src.engines.deepmvp_engine import DeepMVPEngine, OUTPUT_COLUMNS as DEEPMVP_COLUMNS

FASTA_ACC1 = ">ACC1 test protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n"
FASTA_ACC2 = ">ACC2 test protein\nMSDEQKLISEEDLGGVKTHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n"


def _write(directory, name, content):
    p = directory / name
    p.write_text(content)
    return p


def _fake_deepmvp(accession):
    return pd.DataFrame(
        [[accession, "K", 5, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )


def test_batch_corre_todos_los_archivos_reconocidos(tmp_path, monkeypatch):
    input_dir = tmp_path / "batch_in"
    input_dir.mkdir()
    _write(input_dir, "ACC1.fasta", FASTA_ACC1)
    _write(input_dir, "ACC2.fasta", FASTA_ACC2)
    _write(input_dir, "notas.txt", "esto no es un input reconocido, debe ignorarse")
    output_dir = tmp_path / "out"

    monkeypatch.setattr(
        DeepMVPEngine, "run",
        lambda self, items, output_dir=None: [_fake_deepmvp(Path(items[0]).stem.replace("_clean", ""))],
    )

    exit_code = pipeline.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "ACC1_ptm_sites.csv").is_file()
    assert (output_dir / "ACC2_ptm_sites.csv").is_file()

    summary = pd.read_csv(output_dir / "batch_summary.csv")
    assert sorted(summary["archivo"]) == ["ACC1.fasta", "ACC2.fasta"]
    assert set(summary["estado"]) == {"ok"}


def test_batch_un_archivo_invalido_no_detiene_el_resto(tmp_path, monkeypatch):
    input_dir = tmp_path / "batch_in"
    input_dir.mkdir()
    _write(input_dir, "ACC1.fasta", FASTA_ACC1)
    _write(input_dir, "bad.fasta", "no es fasta valido\n")
    output_dir = tmp_path / "out"

    monkeypatch.setattr(
        DeepMVPEngine, "run", lambda self, items, output_dir=None: [_fake_deepmvp("ACC1")]
    )

    exit_code = pipeline.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

    assert exit_code == 1
    assert (output_dir / "ACC1_ptm_sites.csv").is_file()

    summary = pd.read_csv(output_dir / "batch_summary.csv").set_index("archivo")
    assert summary.loc["ACC1.fasta", "estado"] == "ok"
    assert summary.loc["bad.fasta", "estado"] == "error"


def test_batch_excepcion_inesperada_no_pipeline_error_no_tumba_el_batch(tmp_path, monkeypatch):
    # Bug real encontrado en auditoria 2026-08-07: _run_batch solo capturaba PipelineError --
    # cualquier otra excepcion real (aqui simulada con un KeyError generico, el tipo de fallo
    # real que un engine con salida malformada produciria) escapaba del bucle entero,
    # tumbando TODO el batch antes de escribir batch_summary.csv y perdiendo el registro de
    # los archivos previos ya procesados con exito.
    input_dir = tmp_path / "batch_in"
    input_dir.mkdir()
    _write(input_dir, "ACC1.fasta", FASTA_ACC1)
    _write(input_dir, "ACC2.fasta", FASTA_ACC2)
    output_dir = tmp_path / "out"

    def _boom_on_acc2(accession):
        if accession == "ACC2":
            raise KeyError("fallo inesperado simulado, no PipelineError")
        return _fake_deepmvp(accession)

    monkeypatch.setattr(
        DeepMVPEngine, "run",
        lambda self, items, output_dir=None: [_boom_on_acc2(Path(items[0]).stem.replace("_clean", ""))],
    )

    exit_code = pipeline.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

    assert exit_code == 1
    assert (output_dir / "ACC1_ptm_sites.csv").is_file()  # el primer archivo SI se completo

    summary = pd.read_csv(output_dir / "batch_summary.csv").set_index("archivo")
    assert summary.loc["ACC1.fasta", "estado"] == "ok"
    assert summary.loc["ACC2.fasta", "estado"] == "error"
    assert "KeyError" in summary.loc["ACC2.fasta", "detalle"]


def test_batch_directorio_sin_archivos_reconocidos_falla(tmp_path):
    input_dir = tmp_path / "batch_in"
    input_dir.mkdir()
    _write(input_dir, "notas.txt", "nada reconocible aqui")
    output_dir = tmp_path / "out"

    exit_code = pipeline.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

    assert exit_code == 1
    assert not (output_dir / "batch_summary.csv").exists()
