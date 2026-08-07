"""Tests de integracion del orquestador (pipeline.py), las 3 fases.

Fase 1/1.5 corren con logica real (sin mocks, sin binarios externos). Los
motores de Fase 3 (DeepMVP/DeepPTMPred) se mockean a nivel de
``Engine.run`` -- no de subprocess (eso ya lo cubren
test_deepmvp_engine.py/test_deepptmpred_engine.py) -- para no depender de si
estan realmente instalados en la maquina que corre los tests (DeepMVP lo
esta desde 2026-07-28; ver STATUS.md) ni de su tiempo de ejecucion real.
"""

import pandas as pd

import pipeline
from src.engines.deepmvp_engine import DeepMVPEngine, OUTPUT_COLUMNS as DEEPMVP_COLUMNS
from src.engines.deepptmpred_engine import DeepPTMPredEngine, OUTPUT_COLUMNS as DEEPPTMPRED_COLUMNS
from src.engines.fase_a_engine import FaseAEngine
from src.utils.exceptions import DeepMVPExecutionError

FASTA_CONTENT = ">ACC1 test protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n"

PDB_CONTENT = (
    "HEADER    TEST\n"
    "ATOM      1  N   MET A   1      11.104  13.207   2.100  1.00 20.00           N\n"
    "ATOM      2  CA  MET A   1      12.560  13.207   2.100  1.00 20.00           C\n"
    "ATOM      3  N   GLY A   2      14.500  14.700   2.100  1.00 20.00           N\n"
    "ATOM      4  CA  GLY A   2      15.000  15.700   2.100  1.00 20.00           C\n"
    "END\n"
)


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return p


def test_camino_fasta_produce_clean_fasta_y_reporte(tmp_path, monkeypatch):
    input_path = _write(tmp_path, "ACC1.fasta", FASTA_CONTENT)
    output_dir = tmp_path / "out"

    fake_deepmvp = pd.DataFrame(
        [["ACC1", "K", 5, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    monkeypatch.setattr(DeepMVPEngine, "run", lambda self, items, output_dir=None: [fake_deepmvp])

    exit_code = pipeline.main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "ACC1_clean.fasta").is_file()
    report = pd.read_csv(output_dir / "ACC1_ptm_sites.csv")
    assert list(report["accession"]) == ["ACC1"]
    assert list(report["tipo_ptm"]) == ["acetylation_k"]


def test_camino_pdb_produce_reporte_con_consenso(tmp_path, monkeypatch):
    input_path = _write(tmp_path, "1abc.pdb", PDB_CONTENT)
    output_dir = tmp_path / "out"

    fake_deepmvp = pd.DataFrame(
        [["1abc", "M", 1, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    fake_deepptmpred = pd.DataFrame(
        [["1abc", 1, "M", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    monkeypatch.setattr(DeepMVPEngine, "run", lambda self, items, output_dir=None: [fake_deepmvp])
    monkeypatch.setattr(DeepPTMPredEngine, "run", lambda self, items, output_dir=None: [fake_deepptmpred])

    exit_code = pipeline.main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "1abc_derived.fasta").is_file()
    assert (output_dir / "1abc_position_mapping.csv").is_file()
    report = pd.read_csv(output_dir / "1abc_ptm_sites.csv")
    assert report.iloc[0]["motor"] == "DeepMVP+DeepPTMPred"
    assert bool(report.iloc[0]["consenso"]) is True


def test_camino_pdb_reporte_incluye_columnas_de_fase_a(tmp_path, monkeypatch):
    """Fase A conectada 2026-08-03: el sitio de consenso aceptado (acetylation, tipo con
    modulo de Fase A real) debe quedar 'modelado' en el reporte final; FaseAEngine se
    mockea (mismo criterio que DeepMVP/DeepPTMPred aqui: PyRosetta no esta instalado en
    el venv que corre los tests, ver tests/test_fase_a_engine.py para el wiring aislado
    de FaseAEngine en si)."""
    input_path = _write(tmp_path, "1abc.pdb", PDB_CONTENT)
    output_dir = tmp_path / "out"

    fake_deepmvp = pd.DataFrame(
        [["1abc", "M", 1, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    fake_deepptmpred = pd.DataFrame(
        [["1abc", 1, "M", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    monkeypatch.setattr(DeepMVPEngine, "run", lambda self, items, output_dir=None: [fake_deepmvp])
    monkeypatch.setattr(DeepPTMPredEngine, "run", lambda self, items, output_dir=None: [fake_deepptmpred])

    def _fake_fase_a_run(self, items, output_dir=None):
        return [
            {
                "estado": "modelado", "clase": "class1_patch_ddg", "ddg": -1.5,
                "wt_score": 10.0, "mut_score": 8.5, "glycan_tree": None,
                "glygen_evidencia": None, "conjugation_metrics": None,
                "output_pdb": "/tmp/fake.pdb", "error": None,
            }
            for _ in items
        ]

    monkeypatch.setattr(FaseAEngine, "run", _fake_fase_a_run)

    exit_code = pipeline.main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 0
    report = pd.read_csv(output_dir / "1abc_ptm_sites.csv")
    assert report.iloc[0]["fase_a_estado"] == "modelado"
    assert report.iloc[0]["fase_a_ddg"] == -1.5


def test_input_invalido_retorna_codigo_de_error(tmp_path):
    input_path = _write(tmp_path, "bad.fasta", "no es fasta valido\n")
    output_dir = tmp_path / "out"

    exit_code = pipeline.main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 1


def test_motor_no_instalado_falla_con_error_accionable_no_en_silencio(tmp_path, monkeypatch):
    # DeepMVP esta realmente instalado en esta maquina desde 2026-07-28 (ver
    # STATUS.md), asi que el fallo de "motor no instalado" se simula en el
    # motor en vez de depender del estado real de instalacion de la maquina
    # que corre los tests. Lo que se verifica es el comportamiento del
    # pipeline ante un DeepMVPExecutionError: Fase 1 debe completarse igual
    # (no falla en silencio a mitad de camino).
    def _fake_run_no_instalado(self, items, output_dir=None):
        raise DeepMVPExecutionError("No se encontro la instalacion local de DeepMVP.")

    monkeypatch.setattr(DeepMVPEngine, "run", _fake_run_no_instalado)

    input_path = _write(tmp_path, "ACC1.fasta", FASTA_CONTENT)
    output_dir = tmp_path / "out"

    exit_code = pipeline.main(["--input", str(input_path), "--output-dir", str(output_dir)])

    assert exit_code == 1
    assert (output_dir / "ACC1_clean.fasta").is_file()
