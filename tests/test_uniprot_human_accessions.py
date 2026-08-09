"""Tests de src/validation/uniprot_human_accessions.py.

Mismo patron que test_uniprot_localization_client.py: mockea la red, nunca
golpea rest.uniprot.org de verdad en CI/tests (la verificacion contra la API
real -- 20431 accessions, P04637 incluido -- se hizo a mano al implementar,
no se repite en cada corrida de pytest).
"""

import urllib.error
from unittest.mock import MagicMock, patch

import pytest

from src.validation.uniprot_human_accessions import (
    UniProtLookupError,
    fetch_human_reviewed_accessions,
    load_or_fetch_human_reviewed_accessions,
)


def _mock_response(body: str):
    mock = MagicMock()
    mock.read.return_value = body.encode("utf-8")
    mock.__enter__.return_value = mock
    mock.__exit__.return_value = False
    return mock


def test_fetch_parsea_una_accession_por_linea():
    body = "P04637\nP68431\nP62805\n"
    with patch(
        "src.validation.uniprot_human_accessions.urllib.request.urlopen",
        return_value=_mock_response(body),
    ):
        accessions = fetch_human_reviewed_accessions()
    assert accessions == {"P04637", "P68431", "P62805"}


def test_fetch_ignora_lineas_vacias():
    body = "P04637\n\n\nP68431\n"
    with patch(
        "src.validation.uniprot_human_accessions.urllib.request.urlopen",
        return_value=_mock_response(body),
    ):
        accessions = fetch_human_reviewed_accessions()
    assert accessions == {"P04637", "P68431"}


def test_fetch_lanza_en_respuesta_vacia():
    with patch(
        "src.validation.uniprot_human_accessions.urllib.request.urlopen",
        return_value=_mock_response(""),
    ):
        with pytest.raises(UniProtLookupError):
            fetch_human_reviewed_accessions()


def test_fetch_lanza_en_error_http():
    error = urllib.error.HTTPError(
        url="https://rest.uniprot.org/uniprotkb/stream", code=503, msg="err",
        hdrs=None, fp=MagicMock(read=lambda: b"error"),
    )
    with patch(
        "src.validation.uniprot_human_accessions.urllib.request.urlopen",
        side_effect=error,
    ):
        with pytest.raises(UniProtLookupError):
            fetch_human_reviewed_accessions()


def test_fetch_lanza_en_error_de_red():
    error = urllib.error.URLError("DNS failure")
    with patch(
        "src.validation.uniprot_human_accessions.urllib.request.urlopen",
        side_effect=error,
    ):
        with pytest.raises(UniProtLookupError):
            fetch_human_reviewed_accessions()


def test_load_or_fetch_usa_cache_si_existe(tmp_path):
    cache_path = tmp_path / "human_reviewed_accessions.txt"
    cache_path.write_text("P04637\nP68431\n")
    with patch(
        "src.validation.uniprot_human_accessions.fetch_human_reviewed_accessions"
    ) as mock_fetch:
        accessions = load_or_fetch_human_reviewed_accessions(cache_path)
    mock_fetch.assert_not_called()
    assert accessions == {"P04637", "P68431"}


def test_load_or_fetch_descarga_y_cachea_si_no_existe(tmp_path):
    cache_path = tmp_path / "human_reviewed_accessions.txt"
    with patch(
        "src.validation.uniprot_human_accessions.fetch_human_reviewed_accessions",
        return_value={"P04637", "P68431"},
    ) as mock_fetch:
        accessions = load_or_fetch_human_reviewed_accessions(cache_path)
    mock_fetch.assert_called_once()
    assert accessions == {"P04637", "P68431"}
    assert cache_path.is_file()
    assert set(cache_path.read_text().split()) == {"P04637", "P68431"}


def test_load_or_fetch_force_refresh_ignora_cache_existente(tmp_path):
    cache_path = tmp_path / "human_reviewed_accessions.txt"
    cache_path.write_text("P00000\n")
    with patch(
        "src.validation.uniprot_human_accessions.fetch_human_reviewed_accessions",
        return_value={"P04637"},
    ) as mock_fetch:
        accessions = load_or_fetch_human_reviewed_accessions(cache_path, force_refresh=True)
    mock_fetch.assert_called_once()
    assert accessions == {"P04637"}
