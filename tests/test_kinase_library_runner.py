"""Tests de src/engines/_kinase_library_runner.py (solo las piezas testeables sin kinase-library real).

``kl.Substrate``/``kl.get_kinase`` se importan de forma diferida (dentro de
``main()``, nunca a nivel de modulo) -- mismo patron que
``test_emngly_runner.py``/``test_deepptmpred_runner.py``: el modulo en si se
puede importar y probar sin el paquete pesado instalado (vive en el entorno
conda dedicado ``kinase_library``, ver ``Settings.KINASE_LIBRARY_PYTHON_BIN``).
``_score_position`` recibe el modulo ``kl`` como parametro (inyeccion de
dependencia), asi que se puede probar con un doble ligero que reproduce
exactamente la forma real verificada 2026-08-07 contra el paquete real
(``kl.Substrate(seq, phos_pos=pos).predict()`` -> DataFrame indexado por
nombre de quinasa con columnas Score/Score Rank/Percentile/Percentile Rank;
``kl.get_kinase(name).family`` -> string), sin invocar el paquete real.
"""

import pandas as pd
import pytest

from src.engines._kinase_library_runner import OUTPUT_COLUMNS, _score_position, main


class _FakeSubstrate:
    def __init__(self, sequence, phos_pos):
        self.sequence = sequence
        self.phos_pos = phos_pos

    def predict(self):
        if self.sequence[self.phos_pos - 1].upper() not in ("S", "T", "Y"):
            raise Exception("Invalid phosphoacceptor")
        return pd.DataFrame(
            {
                "Score": [5.0385, 4.2377, 3.5045],
                "Score Rank": [1, 2, 4],
                "Percentile": [99.83, 99.77, 99.69],
                "Percentile Rank": [1, 2, 3],
            },
            index=["ATM", "SMG1", "ATR"],
        )


class _FakeKinase:
    def __init__(self, name):
        self._family = {"ATM": "PIKK", "SMG1": "PIKK", "ATR": "PIKK"}[name]

    @property
    def family(self):
        return self._family


class _FakeKinaseLibraryModule:
    Substrate = _FakeSubstrate

    @staticmethod
    def get_kinase(name):
        return _FakeKinase(name)


def test_score_position_devuelve_top_kinasa_familia_percentil_top3():
    scored = _score_position(_FakeKinaseLibraryModule, "A" * 20 + "S" + "A" * 20, 21)

    assert scored == {
        "kinase_library_top_kinase": "ATM",
        "kinase_library_top_family": "PIKK",
        "kinase_library_percentile": 99.83,
        "kinase_library_top3_kinases": "ATM,SMG1,ATR",
    }


def test_score_position_residuo_invalido_lanza():
    with pytest.raises(Exception):
        _score_position(_FakeKinaseLibraryModule, "A" * 40, 21)  # 'A' no es S/T/Y


def test_main_omite_posicion_invalida_y_continua_con_el_resto(monkeypatch, tmp_path):
    # 'kinase_library' se importa con 'import kinase_library as kl' DENTRO de
    # main() -- import mira sys.modules primero, asi que pre-poblarlo con el
    # doble ligero evita necesitar el paquete real instalado en este venv.
    monkeypatch.setitem(__import__("sys").modules, "kinase_library", _FakeKinaseLibraryModule)

    sequence = "A" * 20 + "S" + "A" * 20
    out_csv = tmp_path / "out.csv"
    argv = [
        "_kinase_library_runner.py",
        "--sequence", sequence,
        "--positions", "1", "21",
        "--out-csv", str(out_csv),
    ]
    monkeypatch.setattr("sys.argv", argv)

    exit_code = main()

    assert exit_code == 0
    df = pd.read_csv(out_csv)
    assert list(df.columns) == OUTPUT_COLUMNS
    assert list(df["position"]) == [21]
    assert df.iloc[0]["kinase_library_top_kinase"] == "ATM"
