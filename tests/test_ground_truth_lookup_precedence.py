"""Test de precedencia de pipeline.py::_ground_truth_lookup (decision 2026-08-09):
PANEL (subconjunto curado, verificado contra PDB real) siempre gana sobre el
fallback de dbPTM (base grande, solo lookup) para las proteinas que curra.
"""

from pathlib import Path
from unittest.mock import patch

import pipeline
from src.validation.biological_panel import PANEL


def test_panel_tiene_precedencia_sobre_el_fallback_dbptm():
    entry = PANEL[0]  # p53, PANEL[0] == _P53 en biological_panel.py
    # _ground_truth_lookup matchea contra Path(pdb_filename).stem, no contra el
    # accession UniProt pelado -- "p53_P04637", no "P04637".
    accession = Path(entry.pdb_filename).stem

    with patch("src.validation.dbptm_lookup.lookup_ground_truth") as mock_lookup:
        result = pipeline._ground_truth_lookup(accession)

    mock_lookup.assert_not_called()
    expected = {
        (site.position, site.ptm_type): (site.tier, "curado", site.pmids)
        for site in entry.sites
        if not site.is_negative
    }
    assert result == expected


def test_accession_fuera_del_panel_cae_al_fallback_dbptm():
    with patch(
        "src.validation.dbptm_lookup.lookup_ground_truth",
        return_value={(42, "phosphorylation"): ("B", (12345678,))},
    ) as mock_lookup:
        result = pipeline._ground_truth_lookup("Q99999_no_esta_en_panel")

    mock_lookup.assert_called_once_with("Q99999_no_esta_en_panel")
    assert result == {(42, "phosphorylation"): ("B", "dbPTM", (12345678,))}
