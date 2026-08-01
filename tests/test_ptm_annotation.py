"""Tests del nucleo de Fase 3 (src/engines/ptm_annotation.py). Logica 100% pura.

Los tests de la corroboracion opcional MeToken mockean
``src.engines.ptm_annotation.get_type_corroboration`` (nunca invocan el
subproceso real, ver ``tests/test_metoken_engine.py`` para esos) --
verifican unicamente el WIRING (que columnas se agregan/no, en que filas,
que nunca cambia pasa_umbral/consenso, que un fallo no tumba la anotacion).
"""

import pandas as pd
import pytest

import src.engines.ptm_annotation as ptm_annotation
from src.config.settings import Settings
from src.engines.ptm_annotation import (
    OUTPUT_COLUMNS,
    annotate_fasta_path,
    annotate_pdb_path,
    apply_workflow_filter,
)

DEEPMVP_COLUMNS = ["protein", "aa", "pos", "x", "y_pred", "fpr", "ptm"]
DEEPPTMPRED_COLUMNS = ["protein_id", "position", "residue", "probability", "ptm_type"]


def test_annotate_fasta_path_columnas_y_valores_basicos():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert list(result.columns) == OUTPUT_COLUMNS
    row = result.iloc[0]
    assert row["accession"] == "ACC1"
    assert row["posicion"] == 17
    assert row["tipo_ptm"] == "acetylation_k"
    assert row["motor"] == "DeepMVP"
    assert row["score_deepmvp"] == 0.9
    assert pd.isna(row["score_deepptmpred"])
    assert bool(row["consenso"]) is False
    assert row["camino"] == "FASTA"


def test_annotate_fasta_path_pasa_umbral_usa_fpr():
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],  # pasa (fpr bajo)
            ["ACC1", "K", 30, "xxx", 0.6, 0.5, "acetylation_k"],   # no pasa (fpr alto)
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 40, deepmvp_df)
    assert result.set_index("posicion")["pasa_umbral"].to_dict() == {17: True, 30: False}


def test_annotate_fasta_path_ventana_solo_para_glicosilacion():
    sequence = "AAAANKSAAAAA"  # N en pos 5, X=K, S en pos 7 -> sequon valido
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"],
            ["ACC1", "K", 8, "xxx", 0.9, 0.01, "acetylation_k"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)
    by_pos = result.set_index("posicion")["ventana"].to_dict()
    assert by_pos[5] == "NKS"
    assert by_pos[8] is None


def test_annotate_pdb_path_fusiona_consenso_cuando_ambos_motores_coinciden():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 17, "K", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["motor"] == "DeepMVP+DeepPTMPred"
    assert row["tipo_ptm"] == "acetylation"  # nombre canonico (DeepPTMPred)
    assert row["score_deepmvp"] == 0.9
    assert row["score_deepptmpred"] == 0.8
    assert bool(row["consenso"]) is True


def test_annotate_pdb_path_consenso_false_si_uno_no_pasa_umbral():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(
        # 0.2 bajo el umbral calibrado de acetylation (0.6350621, ver Settings)
        [["ACC1", 17, "K", 0.2, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)

    row = result.iloc[0]
    assert bool(row["consenso"]) is False
    # pasa_umbral sigue True porque DeepMVP si paso (union, no interseccion)
    assert bool(row["pasa_umbral"]) is True


def test_annotate_pdb_path_phosphorylation_y_sin_equivalente_queda_deepmvp_solo():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "Y", 40, "xxx", 0.9, 0.01, "phosphorylation_y"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path("ACC1", "A" * 50, deepmvp_df, deepptmpred_df)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["motor"] == "DeepMVP"
    assert row["tipo_ptm"] == "phosphorylation_y"
    assert bool(row["consenso"]) is False


def test_annotate_pdb_path_n_linked_glycosylation_nunca_hace_consenso():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 25, "N", 0.9, "n_linked_glycosylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "N" * 30, deepmvp_df, deepptmpred_df)

    # Ambos motores reportan la misma posicion, pero NUNCA fusionados en
    # una fila de consenso (decision 2026-08-01, DeepPTMPred no tiene poder
    # discriminativo real para este tipo -- ver STATUS.md).
    assert len(result) == 2
    assert set(result["motor"]) == {"DeepMVP", "DeepPTMPred"}
    assert (result["consenso"] == False).all()


def test_annotate_pdb_path_tipo_exclusivo_deepptmpred_incluido_marcado():
    deepmvp_df = pd.DataFrame(columns=DEEPMVP_COLUMNS)
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 12, "K", 0.7, "malonylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["motor"] == "DeepPTMPred"
    assert row["tipo_ptm"] == "malonylation"
    assert pd.isna(row["score_deepmvp"])
    assert bool(row["consenso"]) is False
    assert bool(row["pasa_umbral"]) is True  # 0.7 >= umbral calibrado de malonylation (0.41699925)


def test_apply_workflow_filter_mantiene_solo_pasa_umbral():
    df = pd.DataFrame(
        [
            {"accession": "A", "posicion": 1, "residuo_wt": "K", "tipo_ptm": "acetylation",
             "motor": "DeepMVP", "score_deepmvp": 0.9, "score_deepptmpred": None,
             "consenso": False, "ventana": None, "camino": "FASTA", "pasa_umbral": True},
            {"accession": "A", "posicion": 2, "residuo_wt": "K", "tipo_ptm": "acetylation",
             "motor": "DeepMVP", "score_deepmvp": 0.3, "score_deepptmpred": None,
             "consenso": False, "ventana": None, "camino": "FASTA", "pasa_umbral": False},
        ],
        columns=OUTPUT_COLUMNS,
    )
    filtered = apply_workflow_filter(df)
    assert filtered["posicion"].tolist() == [1]


def test_annotate_fasta_path_vacio_devuelve_dataframe_vacio_con_columnas():
    result = annotate_fasta_path("ACC1", "AAAA", pd.DataFrame(columns=DEEPMVP_COLUMNS))
    assert list(result.columns) == OUTPUT_COLUMNS
    assert len(result) == 0


# --- Corroboracion opcional de tipo (MeToken, decision 2026-08-01) ---


def test_annotate_pdb_path_sin_pdb_path_no_agrega_columnas_metoken(monkeypatch):
    called = []
    monkeypatch.setattr(ptm_annotation, "get_type_corroboration", lambda *a, **k: called.append(1))

    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 17, "K", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)  # pdb_path=None (default)

    assert list(result.columns) == OUTPUT_COLUMNS  # ninguna columna metoken_* agregada
    assert called == []  # nunca se invoca si no hay pdb_path


def test_annotate_pdb_path_con_pdb_path_agrega_columnas_solo_en_filas_pasa_umbral(monkeypatch, tmp_path):
    def _fake_corroboration(pdb_path, positions, chain_id="A", result_dir=None):
        assert sorted(positions) == [17]  # solo la fila con pasa_umbral=True
        return {17: {"metoken_type": "Acetylation", "metoken_probability": 0.77}}

    monkeypatch.setattr(ptm_annotation, "get_type_corroboration", _fake_corroboration)

    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],  # pasa_umbral=True (fpr bajo)
            ["ACC1", "K", 30, "xxx", 0.9, 0.5, "acetylation_k"],   # pasa_umbral=False (fpr alto, sin consenso)
        ],
        columns=DEEPMVP_COLUMNS,
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "A" * 40, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb", chain_id="A",
    )

    by_pos = result.set_index("posicion")
    assert by_pos.loc[17, "metoken_type"] == "Acetylation"
    assert by_pos.loc[17, "metoken_probability"] == 0.77
    assert by_pos.loc[17, "metoken_type_coincide"] is True  # tipo_ptm='acetylation_k' -> canonico 'acetylation'
    assert by_pos.loc[30, "metoken_type"] is None  # no eligible (pasa_umbral=False), nunca se evalua
    # pasa_umbral/consenso identicos a los que calcularia sin MeToken:
    assert bool(by_pos.loc[17, "pasa_umbral"]) is True
    assert bool(by_pos.loc[30, "pasa_umbral"]) is False


def test_annotate_pdb_path_metoken_type_coincide_false_si_discrepa(monkeypatch, tmp_path):
    monkeypatch.setattr(
        ptm_annotation, "get_type_corroboration",
        lambda *a, **k: {17: {"metoken_type": "Ubiquitination", "metoken_probability": 0.6}},
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "A" * 20, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb",
    )

    assert result.iloc[0]["metoken_type_coincide"] is False


def test_annotate_pdb_path_metoken_type_coincide_none_si_tipo_sin_equivalente(monkeypatch, tmp_path):
    # 'crotonylation' no tiene equivalente en MeToken (ver CANONICAL_TO_METOKEN_TYPE)
    monkeypatch.setattr(
        ptm_annotation, "get_type_corroboration",
        lambda *a, **k: {12: {"metoken_type": "Phosphorylation", "metoken_probability": 0.5}},
    )
    deepmvp_df = pd.DataFrame(columns=DEEPMVP_COLUMNS)
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 12, "K", 0.9, "crotonylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path(
        "ACC1", "A" * 20, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb",
    )

    assert result.iloc[0]["metoken_type"] == "Phosphorylation"  # se reporta igual (informativo)
    assert result.iloc[0]["metoken_type_coincide"] is None  # pero no se puede evaluar coincidencia


def test_annotate_pdb_path_metoken_vacio_deja_columnas_en_none(monkeypatch, tmp_path):
    monkeypatch.setattr(ptm_annotation, "get_type_corroboration", lambda *a, **k: {})  # degradado (no instalado, etc.)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "A" * 20, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb",
    )

    assert result.iloc[0]["metoken_type"] is None
    assert result.iloc[0]["metoken_probability"] is None
    assert result.iloc[0]["metoken_type_coincide"] is None
    assert bool(result.iloc[0]["pasa_umbral"]) is True  # el resto de la anotacion no se ve afectado


def test_annotate_pdb_path_metoken_deshabilitado_via_settings_ignora_pdb_path(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "METOKEN_ENABLED", False)
    called = []
    monkeypatch.setattr(ptm_annotation, "get_type_corroboration", lambda *a, **k: called.append(1))

    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "A" * 20, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb",
    )

    assert list(result.columns) == OUTPUT_COLUMNS  # ninguna columna metoken_* agregada
    assert called == []


def test_annotate_pdb_path_metoken_excepcion_inesperada_no_tumba_la_anotacion(monkeypatch, tmp_path):
    def _boom(*a, **k):
        raise RuntimeError("fallo inesperado")

    monkeypatch.setattr(ptm_annotation, "get_type_corroboration", _boom)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "A" * 20, deepmvp_df, deepptmpred_df, pdb_path=tmp_path / "fake.pdb",
    )

    assert bool(result.iloc[0]["pasa_umbral"]) is True  # la anotacion sigue completa pese al fallo


# --- Corroboracion opcional de N-glicosilacion (StackGlyEmbed, decision 2026-08-01) ---


def test_annotate_fasta_path_enable_stackglyembed_false_no_agrega_columnas(monkeypatch):
    called = []
    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", lambda *a, **k: called.append(1))

    sequence = "AAAANKSAAAAA"  # N en pos 5, sequon NKS valido
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)  # enable_stackglyembed=False (default)

    assert list(result.columns) == OUTPUT_COLUMNS  # ninguna columna stackglyembed_* agregada
    assert called == []


def test_annotate_fasta_path_enable_stackglyembed_true_agrega_columnas_solo_en_nglyco_pasa_umbral(monkeypatch):
    def _fake_corroboration(sequence, positions, **kwargs):
        assert sorted(positions) == [5]  # solo la fila n-glyco con pasa_umbral=True
        return {5: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.93}}

    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", _fake_corroboration)

    sequence = "AAAANKSAAAAAAAAAAAAAAAAAA"
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"],   # pasa_umbral=True, n-glyco
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],    # pasa_umbral=True, NO n-glyco
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df, enable_stackglyembed=True)

    by_pos = result.set_index("posicion")
    assert by_pos.loc[5, "stackglyembed_veredicto"] == "Glicosilado"
    assert by_pos.loc[5, "stackglyembed_score"] == 0.93
    assert by_pos.loc[5, "stackglyembed_coincide"] is True
    assert by_pos.loc[17, "stackglyembed_veredicto"] is None  # no es n-glyco, nunca se evalua
    # pasa_umbral no se ve afectado:
    assert bool(by_pos.loc[5, "pasa_umbral"]) is True
    assert bool(by_pos.loc[17, "pasa_umbral"]) is True


def test_annotate_pdb_path_stackglyembed_no_requiere_pdb_path(monkeypatch):
    """A diferencia de MeToken, StackGlyEmbed solo necesita 'sequence' -- no requiere pdb_path."""
    def _fake_corroboration(sequence, positions, **kwargs):
        return {25: {"stackglyembed_veredicto": "No glicosilado", "stackglyembed_score": 0.08}}

    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", _fake_corroboration)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=None, enable_stackglyembed=True,
    )

    assert result.iloc[0]["stackglyembed_veredicto"] == "No glicosilado"
    assert result.iloc[0]["stackglyembed_score"] == 0.08
    assert result.iloc[0]["stackglyembed_coincide"] is False


def test_annotate_pdb_path_stackglyembed_cubre_ambos_nombres_de_tipo(monkeypatch):
    """tipo_ptm puede ser 'glycosylation_n' (DeepMVP crudo) o 'n_linked_glycosylation' (DeepPTMPred)."""
    def _fake_corroboration(sequence, positions, **kwargs):
        assert sorted(positions) == [10, 25]
        return {
            10: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.88},
            25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.77},
        }

    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", _fake_corroboration)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(
        # 0.999 >= umbral calibrado de n_linked_glycosylation (0.99802846, ver Settings)
        [["ACC1", 10, "N", 0.999, "n_linked_glycosylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
    )

    by_pos = result.set_index("posicion")
    assert by_pos.loc[10, "stackglyembed_coincide"] is True
    assert by_pos.loc[25, "stackglyembed_coincide"] is True


def test_annotate_pdb_path_stackglyembed_vacio_deja_columnas_en_none(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", lambda *a, **k: {})  # degradado

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
    )

    assert result.iloc[0]["stackglyembed_veredicto"] is None
    assert result.iloc[0]["stackglyembed_score"] is None
    assert result.iloc[0]["stackglyembed_coincide"] is None
    assert bool(result.iloc[0]["pasa_umbral"]) is True  # el resto de la anotacion no se ve afectado


def test_annotate_pdb_path_stackglyembed_deshabilitado_via_settings(monkeypatch):
    monkeypatch.setattr(Settings, "STACKGLYEMBED_ENABLED", False)
    called = []
    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", lambda *a, **k: called.append(1))

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
    )

    assert list(result.columns) == OUTPUT_COLUMNS  # ninguna columna stackglyembed_* agregada
    assert called == []


def test_annotate_pdb_path_stackglyembed_excepcion_inesperada_no_tumba_la_anotacion(monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("fallo inesperado")

    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", _boom)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
    )

    assert bool(result.iloc[0]["pasa_umbral"]) is True  # la anotacion sigue completa pese al fallo
