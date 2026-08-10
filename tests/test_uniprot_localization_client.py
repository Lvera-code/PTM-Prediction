"""Tests de src/structural/uniprot_localization_client.py.

Solo usa urllib/json de la stdlib -- se puede testear en la suite principal
mockeando la red (nunca golpea rest.uniprot.org de verdad en CI/tests).
"""

import json
import urllib.error
from unittest.mock import MagicMock, patch

import pytest

from src.structural.uniprot_localization_client import (
    UniProtLookupError,
    lookup_secretory_pathway_evidence,
)


def _mock_response(body: dict):
    mock = MagicMock()
    mock.read.return_value = json.dumps(body).encode("utf-8")
    mock.__enter__.return_value = mock
    mock.__exit__.return_value = False
    return mock


def _http_error(code: int, body: bytes = b"error"):
    return urllib.error.HTTPError(
        url="https://rest.uniprot.org/uniprotkb/X.json", code=code, msg="err",
        hdrs=None, fp=MagicMock(read=lambda: body),
    )


def test_evidencia_via_keyword_secretorio_real():
    # Subconjunto real de la respuesta de UniProt para P01588 (EPO), verificado
    # 2026-08-07 con una consulta en vivo contra rest.uniprot.org.
    body = {
        "keywords": [
            {"name": "Secreted", "category": "Cellular component"},
            {"name": "Signal", "category": "Domain"},
        ],
        "comments": [],
    }
    with patch("src.structural.uniprot_localization_client.urllib.request.urlopen", return_value=_mock_response(body)):
        assert lookup_secretory_pathway_evidence("P01588") is True


def test_evidencia_via_comentario_subcellular_location_sin_keyword():
    body = {
        "keywords": [],
        "comments": [
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [{"location": {"value": "Golgi apparatus membrane"}}],
            }
        ],
    }
    with patch("src.structural.uniprot_localization_client.urllib.request.urlopen", return_value=_mock_response(body)):
        assert lookup_secretory_pathway_evidence("P00000") is True


def test_localizacion_conocida_pero_no_secretora_devuelve_false():
    # Subconjunto real de la respuesta de UniProt para P04406 (GAPDH,
    # citoplasmatica/nuclear), verificado 2026-08-07 -- caso real que
    # descarto la keyword 'Membrane' sola del set (ver docstring del modulo).
    body = {
        "keywords": [
            {"name": "Cytoplasm", "category": "Cellular component"},
            {"name": "Nucleus", "category": "Cellular component"},
            {"name": "Membrane", "category": "Cellular component"},
        ],
        "comments": [
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [{"location": {"value": "Cytoplasm"}}],
            }
        ],
    }
    with patch("src.structural.uniprot_localization_client.urllib.request.urlopen", return_value=_mock_response(body)):
        assert lookup_secretory_pathway_evidence("P04406") is False


def test_sin_datos_de_localizacion_devuelve_none():
    # 'category' real de '3D-structure' es 'Technical term', no 'Cellular
    # component' -- bug real encontrado 2026-08-07 escribiendo este test:
    # una version anterior contaba CUALQUIER keyword (sin filtrar por
    # categoria) como "hay datos de localizacion", dando un falso 'False'
    # en vez del 'None' correcto ("no se pudo determinar nada").
    body = {
        "keywords": [{"name": "3D-structure", "category": "Technical term"}],
        "comments": [{"commentType": "FUNCTION"}],
    }
    with patch("src.structural.uniprot_localization_client.urllib.request.urlopen", return_value=_mock_response(body)):
        assert lookup_secretory_pathway_evidence("P99999") is None


def test_http_400_accession_con_formato_invalido_devuelve_none():
    # Bug real confirmado 2026-08-07: UniProt devuelve HTTP 400 (no 404) para
    # accessions con formato invalido -- el caso MAS COMUN en este proyecto,
    # ya que 'accession' normalmente viene del stem del archivo PDB/FASTA de
    # entrada, casi nunca un ID UniProt real (p.ej. '1qlp'). NO debe lanzar.
    with patch(
        "src.structural.uniprot_localization_client.urllib.request.urlopen",
        side_effect=_http_error(400),
    ):
        assert lookup_secretory_pathway_evidence("1qlp") is None


def test_http_404_accession_no_encontrado_devuelve_none():
    with patch(
        "src.structural.uniprot_localization_client.urllib.request.urlopen",
        side_effect=_http_error(404),
    ):
        assert lookup_secretory_pathway_evidence("NOTREAL999") is None


def test_http_500_lanza_uniprot_lookup_error():
    with patch(
        "src.structural.uniprot_localization_client.urllib.request.urlopen",
        side_effect=_http_error(500),
    ):
        with pytest.raises(UniProtLookupError, match="HTTP 500"):
            lookup_secretory_pathway_evidence("P01588")


def test_url_error_lanza_uniprot_lookup_error():
    with patch(
        "src.structural.uniprot_localization_client.urllib.request.urlopen",
        side_effect=urllib.error.URLError("timed out"),
    ):
        with pytest.raises(UniProtLookupError, match="No se pudo contactar"):
            lookup_secretory_pathway_evidence("P01588")


def test_json_invalido_lanza_uniprot_lookup_error():
    mock = MagicMock()
    mock.read.return_value = b"esto no es json valido {{{"
    mock.__enter__.return_value = mock
    mock.__exit__.return_value = False
    with patch("src.structural.uniprot_localization_client.urllib.request.urlopen", return_value=mock):
        with pytest.raises(UniProtLookupError, match="no parseable"):
            lookup_secretory_pathway_evidence("P01588")
