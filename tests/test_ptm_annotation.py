"""Tests del nucleo de Fase 3 (src/engines/ptm_annotation.py). Logica 100% pura."""

import pandas as pd
import pytest

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
        # 0.2 bajo el umbral calibrado de acetylation (0.6299973, ver Settings)
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
    assert bool(row["pasa_umbral"]) is True  # 0.7 >= umbral calibrado de malonylation (0.4223403)


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
