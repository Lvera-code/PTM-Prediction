"""Tests de src/engines/base_engine.py (contrato abstracto Strategy)."""

import pytest

from src.engines.base_engine import BaseEngine


def test_no_se_puede_instanciar_directamente():
    with pytest.raises(TypeError):
        BaseEngine()


def test_subclase_debe_implementar_run():
    class MotorIncompleto(BaseEngine):
        pass

    with pytest.raises(TypeError):
        MotorIncompleto()


def test_subclase_que_implementa_run_es_instanciable_y_acepta_output_dir():
    class MotorFalso(BaseEngine):
        def run(self, items, output_dir=None):
            return [f"{item}:{output_dir}" for item in items]

    motor = MotorFalso()
    assert motor.run(["a", "b"], output_dir="/tmp/out") == ["a:/tmp/out", "b:/tmp/out"]
    assert motor.run(["a"]) == ["a:None"]  # output_dir es opcional, default None
