"""Tests de src/structural/glygen_client.py.

A diferencia del resto de src/structural/ (requiere pyrosetta), este modulo
solo usa urllib/json de la stdlib -- se puede testear en la suite principal
mockeando la red (nunca golpea api.glygen.org de verdad en CI/tests).
"""

import json
import urllib.error
from unittest.mock import MagicMock, patch

import pytest

from src.structural.glygen_client import (
    GlyGenLookupError,
    fetch_glycosylation_sites,
    lookup_site,
)

# Subconjunto real de la respuesta de GlyGen para P10636 (Tau), verificado
# 2026-07-29 con una consulta en vivo contra api.glygen.org.
FAKE_GLYGEN_RESPONSE = {
    "glycosylation": [
        {
            "type": "O-linked",
            "start_pos": 111,
            "site_category": "predicted",
        },
        {
            "glytoucan_ac": "G49108TO",
            "type": "O-linked",
            "start_pos": 123,
            "site_category": "reported_with_glycan",
        },
        {
            "type": "N-linked",
            "start_pos": 484,
            "site_category": "predicted",
        },
    ]
}


def _mock_response(body: dict, status: int = 200):
    mock = MagicMock()
    mock.read.return_value = json.dumps(body).encode("utf-8")
    mock.__enter__.return_value = mock
    mock.__exit__.return_value = False
    return mock


def test_fetch_glycosylation_sites_devuelve_lista_real():
    with patch("src.structural.glygen_client.urllib.request.urlopen", return_value=_mock_response(FAKE_GLYGEN_RESPONSE)):
        sites = fetch_glycosylation_sites("P10636")
    assert len(sites) == 3
    assert sites[1]["glytoucan_ac"] == "G49108TO"


def test_fetch_glycosylation_sites_accession_no_reconocido_lanza(monkeypatch):
    error_body = {"error_list": [{"error_code": "non-existent-record"}], "reason": {"description": "Invalid accession"}}
    with patch("src.structural.glygen_client.urllib.request.urlopen", return_value=_mock_response(error_body)):
        with pytest.raises(GlyGenLookupError, match="Invalid accession"):
            fetch_glycosylation_sites("P10636-2")


def test_fetch_glycosylation_sites_http_error_se_propaga_como_glygen_error():
    http_error = urllib.error.HTTPError(
        url="https://api.glygen.org/protein/detail/XXXX/", code=500, msg="Internal Error",
        hdrs=None, fp=MagicMock(read=lambda: b"boom"),
    )
    with patch("src.structural.glygen_client.urllib.request.urlopen", side_effect=http_error):
        with pytest.raises(GlyGenLookupError, match="HTTP 500"):
            fetch_glycosylation_sites("XXXX")


def test_fetch_glycosylation_sites_url_error_se_propaga_como_glygen_error():
    url_error = urllib.error.URLError("timed out")
    with patch("src.structural.glygen_client.urllib.request.urlopen", side_effect=url_error):
        with pytest.raises(GlyGenLookupError, match="No se pudo contactar"):
            fetch_glycosylation_sites("P10636")


def test_lookup_site_encuentra_coincidencia_n_linked():
    with patch("src.structural.glygen_client.urllib.request.urlopen", return_value=_mock_response(FAKE_GLYGEN_RESPONSE)):
        site = lookup_site("P10636", 484, "n_linked_glycosylation")
    assert site is not None
    assert site["start_pos"] == 484
    assert site["type"] == "N-linked"


def test_lookup_site_prioriza_glytoucan_ac_si_hay_varios_matches():
    dup_response = {
        "glycosylation": [
            {"type": "O-linked", "start_pos": 123, "site_category": "predicted"},
            {"type": "O-linked", "start_pos": 123, "site_category": "reported_with_glycan", "glytoucan_ac": "G49108TO"},
        ]
    }
    with patch("src.structural.glygen_client.urllib.request.urlopen", return_value=_mock_response(dup_response)):
        site = lookup_site("P10636", 123, "o_linked_glycosylation")
    assert site["glytoucan_ac"] == "G49108TO"


def test_lookup_site_sin_coincidencia_devuelve_none():
    with patch("src.structural.glygen_client.urllib.request.urlopen", return_value=_mock_response(FAKE_GLYGEN_RESPONSE)):
        site = lookup_site("P10636", 999, "n_linked_glycosylation")
    assert site is None


def test_lookup_site_tipo_no_soportado_lanza_value_error():
    with pytest.raises(ValueError, match="no es un tipo de glicosilacion soportado"):
        lookup_site("P10636", 1, "phosphorylation")
