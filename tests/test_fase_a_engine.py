"""Tests de FaseAEngine (src/engines/fase_a_engine.py).

Mockea ``subprocess.run`` (nunca invoca PyRosetta real, no instalado en el
venv principal -- vive en el conda env dedicado ``deepptmpred``). Fase A
NUNCA es un motor de consenso -- cualquier fallo (interprete/runner
ausentes, exit code != 0, timeout, JSON ausente o malformado) debe degradar
a ``estado="error"``/``"no_disponible"``, NUNCA lanzar -- mismo patron que
``test_stackglyembed_engine.py``/``test_metoken_engine.py``.

``FaseAEngine.__init__`` toma ``python_bin``/``runner_script`` como
parametros con default ``Settings.FASE_A_*`` (evaluados UNA vez al definir
la clase, mismo patron que ``DeepPTMPredEngine``) -- por eso estos tests
construyen el engine con argumentos EXPLICITOS en vez de monkeypatchear
``Settings`` y confiar en el default (que no se re-evalua por instancia, ver
``tests/test_deepptmpred_engine.py`` para el mismo criterio ya establecido).
``Settings.FASE_A_ENABLED`` si se lee en vivo dentro de ``run()``, ese si se
puede monkeypatchear.
"""

import json
import subprocess

from src.config.settings import Settings
from src.engines.fase_a_engine import FaseAEngine, FaseASiteRequest

_EMPTY_EXTRA = {
    "clase": None, "ddg": None, "ddg_std": None,
    "wt_score": None, "wt_score_std": None, "mut_score": None, "mut_score_std": None,
    "glycan_tree": None, "glygen_evidencia": None, "conjugation_metrics": None,
    "cadena_tipo_aviso": None,
    "output_pdb": None, "error": None,
}


def _make_engine(tmp_path, python_bin_ok=True, runner_ok=True):
    if python_bin_ok:
        python_bin = tmp_path / "python3"
        python_bin.write_text("#!/bin/sh\n")
        python_bin.chmod(0o755)
    else:
        python_bin = tmp_path / "no_existe_python"

    if runner_ok:
        runner_script = tmp_path / "_fase_a_runner.py"
        runner_script.write_text("# fake\n")
    else:
        runner_script = tmp_path / "no_existe_runner.py"

    return FaseAEngine(python_bin=str(python_bin), runner_script=runner_script)


def _request(**overrides):
    defaults = dict(
        accession="ACC1", pdb_path="/tmp/ACC1_chain_A.pdb", position=24, ptm_type="acetylation",
    )
    defaults.update(overrides)
    return FaseASiteRequest(**defaults)


def test_fase_a_deshabilitado_devuelve_no_disponible_sin_invocar_subprocess(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", False)
    engine = _make_engine(tmp_path)
    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append(1))

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result == [{"estado": "no_disponible", **_EMPTY_EXTRA}]
    assert called == []


def test_items_vacios_devuelve_lista_vacia(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)

    result = engine.run([], output_dir=tmp_path / "out")

    assert result == []


def test_python_bin_ausente_devuelve_no_disponible_para_todos(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path, python_bin_ok=False)

    result = engine.run([_request(), _request(position=50)], output_dir=tmp_path / "out")

    assert result == [{"estado": "no_disponible", **_EMPTY_EXTRA}] * 2


def test_runner_script_ausente_devuelve_no_disponible(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path, runner_ok=False)

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result == [{"estado": "no_disponible", **_EMPTY_EXTRA}]


def test_corrida_exitosa_lee_json_de_salida(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)

    def _fake_run(cmd, **kwargs):
        out_json = cmd[cmd.index("--out-json") + 1]
        payload = {
            "estado": "modelado", "clase": "class1_patch_ddg", "ddg": -3.21,
            "wt_score": 10.0, "mut_score": 6.79, "glycan_tree": None,
            "glygen_evidencia": None, "conjugation_metrics": None,
            "output_pdb": cmd[cmd.index("--out-pdb") + 1], "error": None,
            "position": 24, "ptm_type": "acetylation",
        }
        with open(out_json, "w") as f:
            json.dump(payload, f)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result[0]["estado"] == "modelado"
    assert result[0]["ddg"] == -3.21


def test_cmd_incluye_pdb_position_tipo_y_uniprot_opcional(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)
    captured = []

    def _fake_run(cmd, **kwargs):
        captured.append(cmd)
        out_json = cmd[cmd.index("--out-json") + 1]
        with open(out_json, "w") as f:
            json.dump({"estado": "modelado"}, f)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    engine.run(
        [_request(ptm_type="n_linked_glycosylation", uniprot_accession="P10636")],
        output_dir=tmp_path / "out",
    )

    cmd = captured[0]
    assert cmd[cmd.index("--position") + 1] == "24"
    assert cmd[cmd.index("--ptm-type") + 1] == "n_linked_glycosylation"
    assert cmd[cmd.index("--uniprot-accession") + 1] == "P10636"


def test_uniprot_accession_omitido_si_no_se_da(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)
    captured = []

    def _fake_run(cmd, **kwargs):
        captured.append(cmd)
        out_json = cmd[cmd.index("--out-json") + 1]
        with open(out_json, "w") as f:
            json.dump({"estado": "modelado"}, f)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    engine.run([_request()], output_dir=tmp_path / "out")

    assert "--uniprot-accession" not in captured[0]


def test_exit_code_distinto_de_cero_devuelve_error(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result[0]["estado"] == "error"
    assert "exit code 1" in result[0]["error"]


def test_timeout_devuelve_error(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)

    def _fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=1)

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result[0]["estado"] == "error"
    assert "timeout" in result[0]["error"]


def test_json_ausente_devuelve_error(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)
    monkeypatch.setattr(
        subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
    )

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result[0]["estado"] == "error"
    assert "no generado" in result[0]["error"]


def test_json_malformado_devuelve_error(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "FASE_A_ENABLED", True)
    engine = _make_engine(tmp_path)

    def _fake_run(cmd, **kwargs):
        out_json = cmd[cmd.index("--out-json") + 1]
        with open(out_json, "w") as f:
            f.write("no es json valido")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    result = engine.run([_request()], output_dir=tmp_path / "out")

    assert result[0]["estado"] == "error"
    assert "JSON invalido" in result[0]["error"]
