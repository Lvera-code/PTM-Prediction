"""Tests del nucleo de Fase 3 (src/engines/ptm_annotation.py). Logica 100% pura.

Los tests de la corroboracion opcional MeToken mockean
``src.engines.ptm_annotation.get_type_corroboration`` (nunca invocan el
subproceso real, ver ``tests/test_metoken_engine.py`` para esos) --
verifican unicamente el WIRING (que columnas se agregan/no, en que filas,
que nunca cambia pasa_umbral/consenso, que un fallo no tumba la anotacion).

``_mock_uniprot_lookup`` (autouse): muchos tests de este archivo usan
``glycosylation_n``/``pasa_umbral=True`` (secciones de StackGlyEmbed/EMNGly)
-- sin este fixture, ``_add_secretory_pathway_evidence`` (2026-08-07) haria
una llamada de red REAL a UniProt en cada uno de ellos (accession='ACC1',
inexistente mas alla del contexto del test, pero igual un round-trip de red
real durante la suite principal -- viola la convencion de este proyecto de
nunca golpear una API real en tests, ver docstring de test_glygen_client.py).
Mockea a ``None`` por defecto (mismo significado que "no se pudo verificar"
-- nunca cambia pasa_umbral/consenso); los tests dedicados de esta funcion
mas abajo sobreescriben el mock explicitamente para probar True/False/error.

``_mock_kinase_library`` (autouse, analisis 2026-08-07 punto 5): mismo
motivo -- varios tests usan ``phosphorylation``/``phosphorylation_st``/
``pasa_umbral=True``, y ``KINASE_LIBRARY_PYTHON_BIN`` apunta a un entorno
conda REAL ya instalado en esta maquina (ver Settings), asi que sin este
fixture ``_add_kinase_library_corroboration`` lanzaria un subprocess REAL en
cada uno de ellos durante la suite principal. Mockea a ``{}`` por defecto
(dict vacio, mismo significado que "no disponible"); los tests dedicados mas
abajo sobreescriben el mock explicitamente.
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
    select_fase_a_candidates,
)

DEEPMVP_COLUMNS = ["protein", "aa", "pos", "x", "y_pred", "fpr", "ptm"]
DEEPPTMPRED_COLUMNS = ["protein_id", "position", "residue", "probability", "ptm_type"]


@pytest.fixture(autouse=True)
def _mock_uniprot_lookup(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: None)


@pytest.fixture(autouse=True)
def _mock_kinase_library(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "get_kinase_corroboration", lambda sequence, positions: {})


def test_annotate_fasta_path_columnas_y_valores_basicos(monkeypatch):
    # SECRETORY_PATHWAY_CHECK_ENABLED/PTM_CROSSTALK_CHECK_ENABLED deshabilitados aqui:
    # este test verifica el esquema BASE de columnas, ortogonal a esos 2 avisos
    # informativos (2026-08-07) -- tienen su propio test dedicado mas abajo.
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
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


def _fase_a_row(tipo_ptm, posicion, score_deepmvp=None, score_deepptmpred=None):
    return {
        "accession": "A", "posicion": posicion, "residuo_wt": "K", "tipo_ptm": tipo_ptm,
        "motor": "DeepMVP+DeepPTMPred", "score_deepmvp": score_deepmvp,
        "score_deepptmpred": score_deepptmpred, "consenso": True, "ventana": None,
        "camino": "PDB", "pasa_umbral": True,
    }


def test_select_fase_a_candidates_vacio_si_filtered_vacio():
    df = pd.DataFrame(columns=OUTPUT_COLUMNS)
    result = select_fase_a_candidates(df)
    assert result.empty


def test_select_fase_a_candidates_excluye_tipos_sin_soporte_fase_a():
    df = pd.DataFrame(
        [_fase_a_row("crotonylation", 1, score_deepptmpred=0.9)],  # sin modulo de Fase A
        columns=OUTPUT_COLUMNS,
    )
    result = select_fase_a_candidates(df)
    assert result.empty


def test_select_fase_a_candidates_top1_prioriza_score_deepptmpred():
    df = pd.DataFrame(
        [
            _fase_a_row("acetylation", 10, score_deepmvp=0.99, score_deepptmpred=0.5),
            _fase_a_row("acetylation", 20, score_deepmvp=0.1, score_deepptmpred=0.9),
        ],
        columns=OUTPUT_COLUMNS,
    )
    result = select_fase_a_candidates(df, top_n_per_type=1)
    assert result["posicion"].tolist() == [20]


def test_select_fase_a_candidates_usa_score_deepmvp_si_deepptmpred_ausente():
    df = pd.DataFrame(
        [
            _fase_a_row("phosphorylation", 5, score_deepmvp=0.95, score_deepptmpred=None),
            _fase_a_row("phosphorylation", 6, score_deepmvp=0.2, score_deepptmpred=None),
        ],
        columns=OUTPUT_COLUMNS,
    )
    result = select_fase_a_candidates(df, top_n_per_type=1)
    assert result["posicion"].tolist() == [5]


def test_select_fase_a_candidates_top_n_mayor_a_uno():
    df = pd.DataFrame(
        [
            _fase_a_row("ubiquitination", 1, score_deepptmpred=0.9),
            _fase_a_row("ubiquitination", 2, score_deepptmpred=0.8),
            _fase_a_row("ubiquitination", 3, score_deepptmpred=0.1),
        ],
        columns=OUTPUT_COLUMNS,
    )
    result = select_fase_a_candidates(df, top_n_per_type=2)
    assert sorted(result["posicion"].tolist()) == [1, 2]


def test_select_fase_a_candidates_multiples_tipos_ordenado_por_tipo():
    df = pd.DataFrame(
        [
            _fase_a_row("ubiquitination", 1, score_deepptmpred=0.9),
            _fase_a_row("acetylation", 2, score_deepptmpred=0.9),
        ],
        columns=OUTPUT_COLUMNS,
    )
    result = select_fase_a_candidates(df, top_n_per_type=1)
    assert result["tipo_ptm"].tolist() == ["acetylation", "ubiquitination"]


def test_annotate_fasta_path_vacio_devuelve_dataframe_vacio_con_columnas(monkeypatch):
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
    result = annotate_fasta_path("ACC1", "AAAA", pd.DataFrame(columns=DEEPMVP_COLUMNS))
    assert list(result.columns) == OUTPUT_COLUMNS
    assert len(result) == 0


# --- Corroboracion opcional de tipo (MeToken, decision 2026-08-01) ---


def test_annotate_pdb_path_sin_pdb_path_no_agrega_columnas_metoken(monkeypatch):
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
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


def test_annotate_pdb_path_metoken_corrobora_fila_nglyco_promovida_por_consenso(monkeypatch, tmp_path):
    # Bug real encontrado en auditoria 2026-08-07: MeToken corria ANTES que el consenso
    # de N-glicosilacion (_apply_nglyco_consensus), asi que calculaba sus filas elegibles
    # (pasa_umbral=True) sobre el estado VIEJO -- una fila que DeepMVP solo NO pasaba
    # (fpr alto) pero que EMNGly/StackGlyEmbed promovian a pasa_umbral=True nunca recibia
    # corroboracion MeToken. Este test usa fpr=0.5 (> Settings.DEEPMVP_MAX_FPR, DeepMVP
    # solo = pasa_umbral False) + EMNGly con probabilidad alta (promueve a True) y verifica
    # que, con el orden corregido, MeToken SI corrobora esa fila.
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions",
        lambda *a, **k: {25: {"emngly_probability": 0.95}},
    )
    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", lambda *a, **k: {})
    monkeypatch.setattr(
        ptm_annotation, "get_type_corroboration",
        lambda *a, **k: {25: {"metoken_type": "N-linked Glycosylation", "metoken_probability": 0.8}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.3, 0.5, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    row = result.iloc[0]
    assert bool(row["pasa_umbral"]) is True  # promovida por EMNGly, no por DeepMVP solo
    assert row["metoken_type"] == "N-linked Glycosylation"  # corroboracion SI llego


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
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
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
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
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
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
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


# --- Consenso real de N-glicosilacion en Camino PDB (EMNGly + StackGlyEmbed, decision 2026-08-06) ---


def _position_mapping_df(accession="ACC1", positions=(25,)):
    return pd.DataFrame({
        "accession": [accession] * len(positions),
        "chain_id": ["A"] * len(positions),
        "pdb_seqid": list(positions),
        "insertion_code": [""] * len(positions),
        "fasta_position": list(positions),
        "residue_letter": ["N"] * len(positions),
    })


def test_annotate_pdb_path_nglyco_consenso_3_motores_todos_pasan(monkeypatch, tmp_path):
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions",
        lambda *a, **k: {25: {"emngly_probability": 0.9}},
    )
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    row = result.iloc[0]
    assert row["motor"] == "DeepMVP+EMNGly+StackGlyEmbed"
    assert row["score_emngly"] == 0.9
    assert row["stackglyembed_veredicto"] == "Glicosilado"
    assert bool(row["pasa_umbral"]) is True
    assert bool(row["consenso"]) is True  # 3/3 pasan >= NGLYCO_CONSENSUS_MIN_ENGINES (2)


def test_annotate_pdb_path_nglyco_consenso_canoniza_tipo_ptm_para_fase_a(monkeypatch, tmp_path):
    # Bug real encontrado en auditoria 2026-08-07: _apply_nglyco_consensus dejaba
    # tipo_ptm='glycosylation_n' (nombre crudo de DeepMVP, unico origen posible) sin
    # canonizar a 'n_linked_glycosylation' -- select_fase_a_candidates filtra por
    # Settings.FASE_A_SUPPORTED_PTM_TYPES, que solo contiene el nombre canonico, asi
    # que los sitios de MAYOR confianza (confirmados por 3 motores) quedaban
    # silenciosamente excluidos de Fase A. Este test verifica el fix end-to-end: desde
    # el consenso hasta que select_fase_a_candidates realmente selecciona la fila.
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions",
        lambda *a, **k: {25: {"emngly_probability": 0.9}},
    )
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    assert result.iloc[0]["tipo_ptm"] == "n_linked_glycosylation"

    filtered = apply_workflow_filter(result)
    candidates = select_fase_a_candidates(filtered)
    assert len(candidates) == 1
    assert candidates.iloc[0]["posicion"] == 25


def test_annotate_pdb_path_nglyco_consenso_false_si_solo_1_de_3_pasa(monkeypatch, tmp_path):
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions",
        lambda *a, **k: {25: {"emngly_probability": 0.1}},  # bajo el umbral (0.5)
    )
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "No glicosilado", "stackglyembed_score": 0.2}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    row = result.iloc[0]
    assert row["motor"] == "DeepMVP+EMNGly+StackGlyEmbed"  # los 3 lograron evaluar
    assert bool(row["pasa_umbral"]) is True  # DeepMVP solo ya basta (regla OR)
    assert bool(row["consenso"]) is False  # solo 1/3 pasa su propio umbral


def test_annotate_pdb_path_nglyco_degrada_a_2_motores_si_emngly_no_disponible(monkeypatch, tmp_path):
    monkeypatch.setattr(ptm_annotation, "get_emngly_predictions", lambda *a, **k: {})  # degradado
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    row = result.iloc[0]
    assert row["motor"] == "DeepMVP+StackGlyEmbed"
    assert pd.isna(row["score_emngly"])
    assert bool(row["consenso"]) is True  # DeepMVP + StackGlyEmbed = 2/2 >= minimo


def test_annotate_pdb_path_nglyco_sin_position_mapping_omite_emngly(monkeypatch, tmp_path):
    called = []
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions", lambda *a, **k: called.append(1)
    )
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=None,
    )

    assert called == []  # EMNGly nunca se invoca sin position_mapping
    assert result.iloc[0]["motor"] == "DeepMVP+StackGlyEmbed"


def test_annotate_pdb_path_nglyco_consenso_no_afecta_otros_tipos(monkeypatch, tmp_path):
    monkeypatch.setattr(ptm_annotation, "get_emngly_predictions", lambda *a, **k: {})
    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", lambda *a, **k: {})

    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 17, "K", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path(
        "ACC1", "A" * 10 + "N" * 20, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    acetyl_row = result[result["tipo_ptm"] == "acetylation"].iloc[0]
    assert acetyl_row["motor"] == "DeepMVP+DeepPTMPred"  # fusion normal, sin tocar
    assert bool(acetyl_row["consenso"]) is True


def test_annotate_pdb_path_nglyco_consenso_omite_pathway_generico_stackglyembed(monkeypatch, tmp_path):
    """Con el consenso real activo, NUNCA se invoca el pathway informativo generico (doble subproceso)."""
    calls = []
    monkeypatch.setattr(ptm_annotation, "get_emngly_predictions", lambda *a, **k: {})

    def _spy(*a, **k):
        calls.append(1)
        return {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}}

    monkeypatch.setattr(ptm_annotation, "get_nglyco_corroboration", _spy)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    assert len(calls) == 1  # no 2 -- el pathway generico se omite cuando el consenso real esta activo


def test_annotate_pdb_path_nglyco_emngly_deshabilitado_via_settings_usa_pathway_generico(monkeypatch, tmp_path):
    monkeypatch.setattr(Settings, "EMNGLY_ENABLED", False)
    emngly_called = []
    monkeypatch.setattr(
        ptm_annotation, "get_emngly_predictions", lambda *a, **k: emngly_called.append(1)
    )
    monkeypatch.setattr(
        ptm_annotation, "get_nglyco_corroboration",
        lambda *a, **k: {25: {"stackglyembed_veredicto": "Glicosilado", "stackglyembed_score": 0.85}},
    )

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df, enable_stackglyembed=True,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    assert emngly_called == []
    # Pathway generico (informativo, no decide motor/pasa_umbral/consenso):
    assert result.iloc[0]["motor"] == "DeepMVP"
    assert result.iloc[0]["stackglyembed_veredicto"] == "Glicosilado"
    assert bool(result.iloc[0]["consenso"]) is False


def test_annotate_pdb_path_nglyco_consenso_excepcion_inesperada_no_tumba_la_anotacion(monkeypatch, tmp_path):
    def _boom(*a, **k):
        raise RuntimeError("fallo inesperado")

    monkeypatch.setattr(ptm_annotation, "get_emngly_predictions", _boom)

    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 25, "xxx", 0.95, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path(
        "ACC1", "N" * 30, deepmvp_df, deepptmpred_df,
        pdb_path=tmp_path / "fake.pdb", position_mapping=_position_mapping_df(),
    )

    row = result.iloc[0]
    assert row["motor"] == "DeepMVP"  # sin cambios, el fallo no tumba la anotacion
    assert bool(row["pasa_umbral"]) is True


# --- Corroboracion opcional de VIA SECRETORA (UniProt, analisis 2026-08-07) ---


def test_annotate_fasta_path_via_secretora_evidencia_true(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: True)
    sequence = "AAAANKSAAAAA"  # N en pos 5, sequon NKS valido
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    assert result.iloc[0]["via_secretora_evidencia"] is True


def test_annotate_fasta_path_via_secretora_evidencia_false(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: False)
    sequence = "AAAANKSAAAAA"
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    assert result.iloc[0]["via_secretora_evidencia"] is False


def test_annotate_fasta_path_via_secretora_evidencia_none_por_defecto():
    # Fixture autouse ya mockea a None (equivalente a "no se pudo verificar").
    sequence = "AAAANKSAAAAA"
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    assert "via_secretora_evidencia" in result.columns
    assert result.iloc[0]["via_secretora_evidencia"] is None


def test_annotate_fasta_path_via_secretora_evidencia_no_afecta_filas_no_elegibles(monkeypatch):
    monkeypatch.setattr(ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: True)
    sequence = "AAAKAANKSAAAAA"  # K en pos 4, N en pos 8 (sequon NKS)
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 4, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "N", 8, "xxx", 0.9, 0.01, "glycosylation_n"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    by_pos = result.set_index("posicion")
    assert by_pos.loc[4, "via_secretora_evidencia"] is None  # no es N-glicosilacion
    assert by_pos.loc[8, "via_secretora_evidencia"] is True


def test_annotate_fasta_path_via_secretora_evidencia_sin_filas_elegibles_no_consulta_uniprot(monkeypatch):
    called = []
    monkeypatch.setattr(
        ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: called.append(accession)
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 4, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert called == []
    assert result.iloc[0]["via_secretora_evidencia"] is None


def test_annotate_pdb_path_via_secretora_evidencia_o_linked_fuera_de_alcance(monkeypatch):
    # o_linked_glycosylation deliberadamente fuera de alcance (ver docstring del
    # modulo): dos vias biologicas distintas, este cliente no las distingue.
    called = []
    monkeypatch.setattr(
        ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: called.append(accession)
    )
    deepmvp_df = pd.DataFrame(columns=DEEPMVP_COLUMNS)
    deepptmpred_df = pd.DataFrame(
        [["ACC1", 10, "S", 0.9, "o_linked_glycosylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)  # pdb_path=None

    assert called == []
    assert result.iloc[0]["via_secretora_evidencia"] is None


def test_annotate_fasta_path_via_secretora_evidencia_error_no_tumba_la_anotacion(monkeypatch):
    def _boom(accession):
        raise ptm_annotation.UniProtLookupError("fallo de red simulado")

    monkeypatch.setattr(ptm_annotation, "lookup_secretory_pathway_evidence", _boom)
    sequence = "AAAANKSAAAAA"
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    assert result.iloc[0]["via_secretora_evidencia"] is None  # degrada, no lanza
    assert bool(result.iloc[0]["pasa_umbral"]) is True  # sin cambios en el resto de la fila


def test_annotate_fasta_path_via_secretora_evidencia_deshabilitado_via_settings(monkeypatch):
    monkeypatch.setattr(Settings, "SECRETORY_PATHWAY_CHECK_ENABLED", False)
    called = []
    monkeypatch.setattr(
        ptm_annotation, "lookup_secretory_pathway_evidence", lambda accession: called.append(accession)
    )
    sequence = "AAAANKSAAAAA"
    deepmvp_df = pd.DataFrame(
        [["ACC1", "N", 5, "xxx", 0.9, 0.01, "glycosylation_n"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", sequence, deepmvp_df)

    assert called == []
    assert "via_secretora_evidencia" not in result.columns


# --- Corroboracion opcional de ESPECIFICIDAD DE QUINASA (Kinase Library, analisis 2026-08-07) ---


def test_annotate_fasta_path_kinase_library_agrega_columnas(monkeypatch):
    monkeypatch.setattr(
        ptm_annotation, "get_kinase_corroboration",
        lambda sequence, positions: {
            10: {
                "kinase_library_top_kinase": "ATM", "kinase_library_top_family": "PIKK",
                "kinase_library_percentile": 99.83, "kinase_library_top3_kinases": "ATM,SMG1,ATR",
            },
        },
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "S", 10, "xxx", 0.9, 0.01, "phosphorylation_st"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    row = result.iloc[0]
    assert row["kinase_library_top_kinase"] == "ATM"
    assert row["kinase_library_top_family"] == "PIKK"
    assert row["kinase_library_percentile"] == 99.83
    assert row["kinase_library_top3_kinases"] == "ATM,SMG1,ATR"


def test_annotate_fasta_path_kinase_library_vacio_deja_columnas_en_none():
    # Fixture autouse ya mockea a {} (equivalente a "no disponible").
    deepmvp_df = pd.DataFrame(
        [["ACC1", "S", 10, "xxx", 0.9, 0.01, "phosphorylation_st"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert "kinase_library_top_kinase" in result.columns
    assert result.iloc[0]["kinase_library_top_kinase"] is None


def test_annotate_fasta_path_kinase_library_no_afecta_filas_no_elegibles(monkeypatch):
    monkeypatch.setattr(
        ptm_annotation, "get_kinase_corroboration",
        lambda sequence, positions: {
            8: {
                "kinase_library_top_kinase": "ATM", "kinase_library_top_family": "PIKK",
                "kinase_library_percentile": 99.83, "kinase_library_top3_kinases": "ATM,SMG1,ATR",
            },
        },
    )
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 4, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "S", 8, "xxx", 0.9, 0.01, "phosphorylation_st"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    by_pos = result.set_index("posicion")
    assert by_pos.loc[4, "kinase_library_top_kinase"] is None  # no es fosforilacion
    assert by_pos.loc[8, "kinase_library_top_kinase"] == "ATM"


def test_annotate_fasta_path_kinase_library_sin_filas_elegibles_no_invoca(monkeypatch):
    called = []
    monkeypatch.setattr(
        ptm_annotation, "get_kinase_corroboration",
        lambda sequence, positions: called.append(positions),
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 4, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert called == []
    assert result.iloc[0]["kinase_library_top_kinase"] is None


def test_annotate_pdb_path_kinase_library_no_requiere_pdb_path(monkeypatch):
    # Igual que StackGlyEmbed (ver seccion arriba): solo necesita la secuencia,
    # aplica aunque pdb_path sea None.
    monkeypatch.setattr(
        ptm_annotation, "get_kinase_corroboration",
        lambda sequence, positions: {
            10: {
                "kinase_library_top_kinase": "ATM", "kinase_library_top_family": "PIKK",
                "kinase_library_percentile": 99.83, "kinase_library_top3_kinases": "ATM,SMG1,ATR",
            },
        },
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "S", 10, "xxx", 0.9, 0.01, "phosphorylation_st"]], columns=DEEPMVP_COLUMNS
    )
    deepptmpred_df = pd.DataFrame(columns=DEEPPTMPRED_COLUMNS)
    result = annotate_pdb_path("ACC1", "A" * 20, deepmvp_df, deepptmpred_df)  # pdb_path=None

    assert result.iloc[0]["kinase_library_top_kinase"] == "ATM"


def test_annotate_fasta_path_kinase_library_deshabilitado_via_settings(monkeypatch):
    monkeypatch.setattr(Settings, "KINASE_LIBRARY_ENABLED", False)
    called = []
    monkeypatch.setattr(
        ptm_annotation, "get_kinase_corroboration",
        lambda sequence, positions: called.append(positions),
    )
    deepmvp_df = pd.DataFrame(
        [["ACC1", "S", 10, "xxx", 0.9, 0.01, "phosphorylation_st"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert called == []
    assert "kinase_library_top_kinase" not in result.columns


def test_annotate_fasta_path_kinase_library_excepcion_inesperada_no_tumba_la_anotacion(monkeypatch):
    def _boom(sequence, positions):
        raise RuntimeError("fallo inesperado simulado")

    monkeypatch.setattr(ptm_annotation, "get_kinase_corroboration", _boom)
    deepmvp_df = pd.DataFrame(
        [["ACC1", "S", 10, "xxx", 0.9, 0.01, "phosphorylation_st"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 20, deepmvp_df)

    assert bool(result.iloc[0]["pasa_umbral"]) is True  # sin cambios en el resto de la fila


# --- Aviso de competencia/crosstalk entre PTMs del mismo residuo (analisis 2026-08-07) ---


def test_annotate_fasta_path_crosstalk_aviso_dos_tipos_compiten_mismo_residuo():
    # Acetilacion y ubiquitinacion compiten por el mismo grupo epsilon-amino
    # de la Lys17 -- ambas reales, DeepMVP puede proponer ambas de forma
    # independiente (2 clasificadores separados por tipo).
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "ubiquitination_k"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    avisos = sorted(result["ptm_crosstalk_aviso"].tolist())
    assert all(a is not None for a in avisos)
    assert "acetylation" in avisos[0] or "ubiquitination" in avisos[0]
    tipos_en_avisos = set(result["tipo_ptm"]) 
    assert tipos_en_avisos == {"acetylation_k", "ubiquitination_k"}


def test_annotate_fasta_path_crosstalk_sin_aviso_un_solo_tipo():
    deepmvp_df = pd.DataFrame(
        [["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"]], columns=DEEPMVP_COLUMNS
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    assert result.iloc[0]["ptm_crosstalk_aviso"] is None


def test_annotate_fasta_path_crosstalk_sin_aviso_si_uno_no_pasa_umbral():
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],   # pasa (fpr bajo)
            ["ACC1", "K", 17, "xxx", 0.9, 0.99, "ubiquitination_k"],  # no pasa (fpr alto)
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    assert result["ptm_crosstalk_aviso"].isna().all()


def test_annotate_fasta_path_crosstalk_sin_aviso_posiciones_distintas():
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 5, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "K", 20, "xxx", 0.9, 0.01, "ubiquitination_k"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    assert result["ptm_crosstalk_aviso"].isna().all()


def test_annotate_fasta_path_crosstalk_valida_residuo_wt():
    # 'acetylation_k' con residuo_wt distinto de 'K' (no deberia pasar en la
    # practica -- DeepMVP reporta el residuo real -- pero la funcion debe
    # ignorar la fila si el residuo no coincide, nunca asumir por el tipo).
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "X", 17, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "ubiquitination_k"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    assert result["ptm_crosstalk_aviso"].isna().all()  # solo 1 miembro real del grupo (K)


def test_annotate_fasta_path_crosstalk_deshabilitado_via_settings(monkeypatch):
    monkeypatch.setattr(Settings, "PTM_CROSSTALK_CHECK_ENABLED", False)
    deepmvp_df = pd.DataFrame(
        [
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "acetylation_k"],
            ["ACC1", "K", 17, "xxx", 0.9, 0.01, "ubiquitination_k"],
        ],
        columns=DEEPMVP_COLUMNS,
    )
    result = annotate_fasta_path("ACC1", "A" * 30, deepmvp_df)

    assert "ptm_crosstalk_aviso" not in result.columns
