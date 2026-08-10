"""Cliente minimo para UniProt REST (evidencia real de localizacion subcelular).

## Por que existe

La N-glicosilacion ocurre quimicamente en el lumen del reticulo
endoplasmatico/aparato de Golgi -- una proteina puramente citoplasmatica o
nuclear (sin senal de secrecion/membrana real) practicamente nunca se
glicosila en un sequon N-X-[S/T], sin importar cuan fuerte sea el score de
consenso del pipeline (DeepMVP/DeepPTMPred/EMNGly son todos predictores de
secuencia/estructura, ninguno modela la via biosintetica real del
sustrato). Hallazgo del analisis de coherencia biologica 2026-08-07: el
pipeline no comprobaba esto en ningun punto. Este cliente consulta la
anotacion REAL de localizacion subcelular de UniProt (curada, con evidencia
citada) para dar una senal informativa -- NUNCA decide ``pasa_umbral``/
``consenso``, solo avisa cuando el consenso ya acepto un sitio N-glico sin
evidencia conocida de via secretora.

## Endpoint real (verificado 2026-08-07 con una consulta real contra P01588/
   EPO, proteina secretada conocida -- no asumido de la documentacion)

``https://rest.uniprot.org/uniprotkb/{accession}.json`` (GET, sin auth).
Dos senales reales en la respuesta, ambas comprobadas con la consulta real:
``keywords[].name`` (vocabulario controlado, p.ej. ``"Secreted"``/``"Signal"``
-- confirmado presente para P01588) y
``comments[commentType="SUBCELLULAR LOCATION"].subcellularLocations[].location.value``
(nota: el ``commentType`` real tiene un ESPACIO, ``"SUBCELLULAR LOCATION"``,
no guion bajo -- el primer intento asumiendolo con guion bajo devolvio 0
resultados contra un caso real conocido, corregido tras verificar la
respuesta cruda).

Nunca aborta el flujo principal: un accession no reconocido por UniProt (muy
comun aqui, ya que ``accession`` en este proyecto suele derivarse del nombre
de archivo PDB/FASTA de entrada, no siempre un ID UniProt real) o un fallo de
red se reportan como "no se pudo verificar" (``None``), NUNCA como "no
secretado" -- ausencia de evidencia no es evidencia de ausencia.
"""

import json
import urllib.error
import urllib.request
from typing import Optional

UNIPROT_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_TIMEOUT_SECONDS = 20

# Vocabulario controlado REAL de UniProt (keywords de categoria "Cellular
# component" mas "Signal", confirmado contra la respuesta real de P01588)
# consistente con que la proteina transite la via secretora (RE/Golgi ->
# superficie/extracelular/lisosoma), donde la N-glicosilacion es
# quimicamente posible. No es una lista exhaustiva de doctrina de biologia
# celular -- es deliberadamente conservadora (mejor un falso "sin evidencia"
# que nunca ocurre aqui, ya que solo se activa esta funcion si UniProt SI
# tiene datos de localizacion, ver ``lookup_secretory_pathway_evidence``).
#
# 'Membrane' (sin calificar) se probo y DESCARTO deliberadamente 2026-08-07:
# contra un caso real (GAPDH/P04406, proteina citoplasmatica/nuclear con
# roles moonlighting documentados de asociacion perifierica de membrana)
# esa keyword sola daba un falso positivo real -- demasiado amplia, cubre
# asociacion periferica ademas de topologia de membrana integral real.
SECRETORY_PATHWAY_KEYWORDS = {
    "Secreted", "Signal", "Signal-anchor", "Cell membrane",
    "Transmembrane", "Golgi apparatus", "Endoplasmic reticulum",
    "Extracellular matrix", "Lysosome", "Cell surface", "Cell projection",
    "Postsynaptic cell membrane", "Presynaptic cell membrane",
}


class UniProtLookupError(Exception):
    """Fallo al consultar la API de UniProt: red, timeout, o respuesta no parseable."""


def lookup_secretory_pathway_evidence(accession: str) -> Optional[bool]:
    """Evidencia real de UniProt de que ``accession`` transita la via secretora.

    Args:
        accession: Accession UniProt (canonico, sin sufijo de isoforma). Si
            no es un ID UniProt real (comun en este proyecto, ver docstring
            del modulo), UniProt devuelve HTTP 404 -- se traduce a ``None``,
            no a ``False``.

    Returns:
        ``True`` si UniProt reporta evidencia de localizacion consistente
        con la via secretora (keyword o comentario de localizacion
        subcelular). ``False`` si UniProt SI tiene datos reales de
        localizacion pero NINGUNO es consistente con la via secretora
        (evidencia real en contra, no solo ausencia). ``None`` si UniProt no
        tiene ningun dato de localizacion para este accession, o si no se
        pudo determinar si el accession es valido (ver ``UniProtLookupError``
        mas abajo, que el llamador debe capturar por separado).

    Raises:
        UniProtLookupError: Si la peticion de red falla (timeout, DNS) o la
            respuesta no es JSON valido. Un accession no reconocido NO lanza
            -- se trata como "sin datos", devuelve ``None`` (ver arriba).
            Verificado real 2026-08-07: UniProt devuelve HTTP 404 para un
            accession con formato valido pero inexistente, Y TAMBIEN HTTP
            400 para un accession con formato invalido (el caso MAS COMUN
            en este proyecto, ya que ``accession`` normalmente se deriva del
            nombre del archivo PDB/FASTA de entrada, casi nunca un ID
            UniProt real) -- ambos codigos se tratan igual aqui, ninguno
            lanza. Solo errores de servidor/red genuinos (5xx, timeout, DNS)
            lanzan.
    """
    url = UNIPROT_URL_TEMPLATE.format(accession=accession)
    request = urllib.request.Request(url, headers={"Accept": "application/json"})

    try:
        with urllib.request.urlopen(request, timeout=UNIPROT_TIMEOUT_SECONDS) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 404):
            return None
        raise UniProtLookupError(
            f"UniProt devolvio HTTP {exc.code} para '{accession}'."
        ) from exc
    except urllib.error.URLError as exc:
        raise UniProtLookupError(
            f"No se pudo contactar a UniProt para '{accession}': {exc.reason}"
        ) from exc
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise UniProtLookupError(
            f"Respuesta de UniProt no parseable como JSON para '{accession}': {exc}"
        ) from exc

    keyword_names = {kw.get("name") for kw in data.get("keywords", [])}
    if keyword_names & SECRETORY_PATHWAY_KEYWORDS:
        return True

    # Bug real encontrado 2026-08-07 escribiendo el test: CUALQUIER keyword
    # (p.ej. '3D-structure', sin relacion con localizacion) contaba como
    # "hay datos de localizacion" en una version anterior de este codigo --
    # 'category' real de UniProt distingue keywords de localizacion
    # ('Cellular component') de las demas (PTM/Disease/Technical term/etc,
    # confirmado en la respuesta real de P01588/EPO).
    has_localization_data = any(
        kw.get("category") == "Cellular component" for kw in data.get("keywords", [])
    )
    for comment in data.get("comments", []):
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc in comment.get("subcellularLocations", []):
            value = loc.get("location", {}).get("value", "")
            has_localization_data = True
            if any(kw in value for kw in SECRETORY_PATHWAY_KEYWORDS):
                return True

    if has_localization_data:
        return False
    return None
