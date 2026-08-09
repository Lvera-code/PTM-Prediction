"""Lista de accessions UniProt humanas *reviewed* (Swiss-Prot).

## Por que existe

El panel de validacion biologica migra de curacion manual a dbPTM (decision
2026-08-09, vault) -- dbPTM es multi-especie, sin filtro de organismo en sus
archivos de descarga (verificado 2026-08-09 leyendo la pagina de descarga real,
no asumido). Este modulo resuelve ese filtro consultando UniProt directamente,
mismo patron de cliente minimo via ``urllib`` que ``uniprot_localization_client.py``.

## Endpoint real (verificado 2026-08-09 con una consulta real, no asumido de
   la documentacion)

``https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=list``
(GET, sin auth). El endpoint ``/stream`` devuelve el resultado completo en una
sola respuesta sin paginar -- confirmado real: 20431 accessions, un accession
por linea, texto plano. No usar ``/search`` (pagina de a 500 por defecto).
"""

import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Set

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_QUERY = "organism_id:9606 AND reviewed:true"
UNIPROT_TIMEOUT_SECONDS = 60


class UniProtLookupError(Exception):
    """Fallo al consultar la API de UniProt: red, timeout, o respuesta vacia."""


def fetch_human_reviewed_accessions() -> Set[str]:
    """Descarga la lista completa de accessions humanas *reviewed* de UniProt.

    Returns:
        Set de accessions UniProt (p.ej. ``{"P04637", "P68431", ...}``),
        ~20000 elementos (confirmado real 2026-08-09).

    Raises:
        UniProtLookupError: fallo de red/timeout, o respuesta vacia/no
            parseable -- nunca devuelve un set vacio silenciosamente, ya que
            eso haria que el importador de dbPTM descarte TODO por "no
            humano", un fallo silencioso peor que abortar.
    """
    url = f"{UNIPROT_STREAM_URL}?query={urllib.parse.quote(UNIPROT_QUERY)}&format=list"
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})

    try:
        with urllib.request.urlopen(request, timeout=UNIPROT_TIMEOUT_SECONDS) as response:
            body = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raise UniProtLookupError(
            f"UniProt devolvio HTTP {exc.code} consultando accessions humanas reviewed."
        ) from exc
    except urllib.error.URLError as exc:
        raise UniProtLookupError(
            f"No se pudo contactar a UniProt para la lista de accessions humanas: {exc.reason}"
        ) from exc

    accessions = {line.strip() for line in body.splitlines() if line.strip()}
    if not accessions:
        raise UniProtLookupError(
            "UniProt devolvio una respuesta vacia para la lista de accessions humanas reviewed."
        )
    return accessions


def load_or_fetch_human_reviewed_accessions(
    cache_path: Path,
    force_refresh: bool = False,
) -> Set[str]:
    """Lee ``cache_path`` si existe (una accession por linea); si no, la descarga y cachea.

    Sin expiracion automatica -- es un fetch de setup unico, se refresca a
    mano con ``force_refresh=True`` si se quiere una foto mas reciente del
    proteoma humano reviewed.
    """
    if cache_path.is_file() and not force_refresh:
        return {line.strip() for line in cache_path.read_text().splitlines() if line.strip()}

    accessions = fetch_human_reviewed_accessions()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("\n".join(sorted(accessions)) + "\n")
    return accessions
