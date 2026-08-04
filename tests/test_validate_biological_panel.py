"""Tests del runner de validacion biologica (scripts/validate_biological_panel.py). Las
funciones puras (_recall_by_tier/_negative_control_report) se testean con un PanelEntry
sintetico; main() end-to-end se testea con un PDB real pequeño del panel pero motores
mockeados (mismo criterio que tests/test_pipeline_fase1.py -- ninguno de los engines esta
mockeado a nivel de subprocess, DeepMVPEngine/DeepPTMPredEngine.run se reemplazan enteros).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.validate_biological_panel import _negative_control_report, _recall_by_tier, main
from src.engines.deepmvp_engine import DeepMVPEngine, OUTPUT_COLUMNS as DEEPMVP_COLUMNS
from src.engines.deepptmpred_engine import DeepPTMPredEngine, OUTPUT_COLUMNS as DEEPPTMPRED_COLUMNS
from src.validation.biological_panel import PANEL, GroundTruthSite, PanelEntry

_FAKE_ENTRY = PanelEntry(
    name="fake", uniprot_accession="X00000", pdb_filename="fake.pdb", length=10,
    sites=(
        GroundTruthSite(1, "K", "acetylation", "A", (1,)),
        GroundTruthSite(2, "K", "acetylation", "A", (2,)),
        GroundTruthSite(3, "S", "phosphorylation", "B", (3,)),
        GroundTruthSite(4, "N", "n_linked_glycosylation", "A", (4,), is_negative=True),
    ),
)


def test_recall_by_tier_cuenta_aciertos_por_separado():
    accepted = {(1, "acetylation"), (3, "phosphorylation")}  # acierta 1/2 de A, 1/1 de B

    recall = _recall_by_tier(_FAKE_ENTRY, accepted)

    assert recall["A"] == (1, 2)
    assert recall["B"] == (1, 1)


def test_negative_control_report_marca_falso_positivo():
    lines_ok = _negative_control_report(_FAKE_ENTRY, accepted=set())
    assert "correctamente NO aceptado" in lines_ok[0]

    lines_bad = _negative_control_report(_FAKE_ENTRY, accepted={(4, "n_linked_glycosylation")})
    assert "FALSO POSITIVO" in lines_bad[0]


def test_only_filtra_por_nombre_y_falla_si_no_hay_coincidencias(capsys):
    exit_code = main(["--only", "nombre_que_no_existe"])
    assert exit_code == 1


def test_main_end_to_end_con_motores_mockeados(tmp_path, monkeypatch):
    """Corrida real de Fase 1.5/2/3 (sin mock de subprocess) sobre un PDB real del panel
    (histone_h4, el mas pequeño), motores mockeados para no depender de que
    DeepMVP/DeepPTMPred esten instalados en el entorno que corre pytest."""
    entry = next(e for e in PANEL if e.name == "histone_h4")

    fake_deepmvp = pd.DataFrame(
        [[entry.uniprot_accession, "K", 6, "xxx", 0.9, 0.01, "acetylation_k"]],
        columns=DEEPMVP_COLUMNS,
    )
    fake_deepptmpred = pd.DataFrame(
        [[entry.uniprot_accession, 6, "K", 0.8, "acetylation"]], columns=DEEPPTMPRED_COLUMNS
    )
    monkeypatch.setattr(DeepMVPEngine, "run", lambda self, items, output_dir=None: [fake_deepmvp])
    monkeypatch.setattr(DeepPTMPredEngine, "run", lambda self, items, output_dir=None: [fake_deepptmpred])

    exit_code = main(["--only", "histone_h4", "--output-dir", str(tmp_path)])

    assert exit_code == 0
