"""Cliente minimo para la API real de GlyGen (api.glygen.org).

## Por que existe

Ni DeepMVP ni DeepPTMPred predicen la COMPOSICION real del glicano
(glicoforma) en un sitio de glicosilacion, solo el SITIO -- por eso
``pyrosetta_glycan_patch.py`` usa el nucleo biosintetico conservado
(``N-glycan_core`` / ``core_1_O-glycan``) como default DOCUMENTADO, no una
prediccion real (ver su docstring). Mejora identificada 2026-07-28,
implementada aqui 2026-07-29: antes de caer a ese default generico, consultar
GlyGen -- si la proteina ya tiene evidencia EXPERIMENTAL real (no solo
predicha por otra herramienta) de glicosilacion en ese sitio exacto, vale la
pena reportarlo como corroboracion, aunque el glicano especifico que reporta
GlyGen (identificado por su propio ID GlyTouCan) no se traduzca aqui a un
arbol PyRosetta concreto -- esa traduccion GlyTouCan -> nucleo IUPAC
construible por PyRosetta es un problema real de mapeo de nomenclatura de
glicanos, no resuelto por este cliente, fuera de alcance de esta mejora.

## Endpoint real (verificado 2026-07-29, no asumido)

Descubierto leyendo ``https://api.glygen.org/swagger.json`` directamente (no
hay un endpoint dedicado ``/glycosylation/...`` como se asumia originalmente
-- la info vive dentro de la respuesta de ``/protein/detail/{accession}/``,
un POST con body ``{"uniprot_canonical_ac": accession}``). Verificado con una
consulta real en vivo contra P10636 (Tau, mismo caso de prueba que el resto
de Fase A): 14 sitios reales devueltos, con ``site_category`` distinguiendo
``"predicted"`` (solo herramienta computacional, p. ej. ISOGlyP) de
``"reported_with_glycan"`` (evidencia real con ``glytoucan_ac`` asociado,
referencias a PubMed/O-GlcNAc Atlas/GlyConnect). El accession debe ser el
UniProt canonico SIN sufijo de isoforma (``"P10636"``, no ``"P10636-2"``) --
verificado empiricamente: la version con isoforma devuelve
``{"error_list": [{"error_code": "non-existent-record"}]}`` (HTTP 500) para
este caso real, mientras que la version sin sufijo si funciona (HTTP 200).

100% consulta de red (GlyGen no tiene un dump local practico para esto) --
a diferencia del resto del pipeline, este modulo SI depende de conectividad
externa. Fallos de red se propagan como :class:`GlyGenLookupError`, nunca
como una excepcion generica -- quien llama decide si es fatal o solo se
reporta como "sin corroboracion disponible" (ver
``pyrosetta_glycan_patch.py::check_glygen_evidence``, que degrada
graciosamente).
"""

import json
import urllib.error
import urllib.request
from typing import Dict, List, Optional

GLYGEN_DETAIL_URL_TEMPLATE = "https://api.glygen.org/protein/detail/{accession}/"
GLYGEN_TIMEOUT_SECONDS = 20

# Mapeo tipo interno del pipeline -> valor real del campo 'type' de GlyGen
# (confirmado en la respuesta real, no inventado: "N-linked" / "O-linked").
PTM_TYPE_TO_GLYGEN_TYPE = {
    "n_linked_glycosylation": "N-linked",
    "o_linked_glycosylation": "O-linked",
}


class GlyGenLookupError(Exception):
    """Fallo al consultar la API de GlyGen: red, HTTP no-200, o accession no reconocido."""


def fetch_glycosylation_sites(uniprot_accession: str) -> List[Dict]:
    """Consulta GlyGen y devuelve la lista cruda de sitios de glicosilacion reportados.

    Args:
        uniprot_accession: Accession UniProt CANONICO, sin sufijo de isoforma
            (p. ej. ``"P10636"``, no ``"P10636-2"`` -- ver docstring del
            modulo, verificado empiricamente que el sufijo produce un error
            "non-existent-record" real).

    Returns:
        Lista de dicts, uno por sitio reportado (puede estar vacia si la
        proteina no tiene ningun sitio de glicosilacion en GlyGen). Cada dict
        conserva el esquema real de la API (``type``, ``start_pos``,
        ``site_category``, ``glytoucan_ac`` opcional, ``evidence``, etc.).

    Raises:
        GlyGenLookupError: Si la peticion de red falla (timeout, DNS, HTTP
            no-200) o el accession no es reconocido por GlyGen.
    """
    url = GLYGEN_DETAIL_URL_TEMPLATE.format(accession=uniprot_accession)
    payload = json.dumps({"uniprot_canonical_ac": uniprot_accession}).encode("utf-8")
    request = urllib.request.Request(
        url, data=payload, headers={"Content-Type": "application/json"}, method="POST",
    )

    try:
        with urllib.request.urlopen(request, timeout=GLYGEN_TIMEOUT_SECONDS) as response:
            body = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[:500]
        raise GlyGenLookupError(
            f"GlyGen devolvio HTTP {exc.code} para '{uniprot_accession}': {detail}"
        ) from exc
    except urllib.error.URLError as exc:
        raise GlyGenLookupError(
            f"No se pudo contactar a GlyGen para '{uniprot_accession}': {exc.reason}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise GlyGenLookupError(
            f"Respuesta de GlyGen para '{uniprot_accession}' no es JSON valido: {exc}"
        ) from exc

    if "error_list" in body:
        reason = body.get("reason", {}).get("description", body.get("error_list"))
        raise GlyGenLookupError(f"GlyGen no reconoce el accession '{uniprot_accession}': {reason}")

    return body.get("glycosylation", [])


def lookup_site(uniprot_accession: str, position: int, ptm_type: str) -> Optional[Dict]:
    """Busca evidencia GlyGen para un sitio (posicion 1-based) y tipo de glicosilacion dados.

    Args:
        uniprot_accession: Ver :func:`fetch_glycosylation_sites`.
        position: Posicion 1-based en la secuencia (mismo sistema de
            numeracion que ``pyrosetta_glycan_patch.attach_glycan``).
        ptm_type: ``"n_linked_glycosylation"`` u ``"o_linked_glycosylation"``
            (mismos valores que ``GLYCAN_TREE_BY_TYPE`` en
            ``pyrosetta_glycan_patch.py``).

    Returns:
        El registro de GlyGen (dict con el esquema real de la API) que mejor
        evidencia tiene entre los que coinciden en posicion+tipo (prioriza
        uno con ``glytoucan_ac`` real -- evidencia experimental con glicano
        identificado -- sobre uno solo ``"predicted"``), o ``None`` si GlyGen
        no reporta ningun sitio en esa posicion+tipo exactos.

    Raises:
        ValueError: Si ``ptm_type`` no es un tipo de glicosilacion soportado.
        GlyGenLookupError: Ver :func:`fetch_glycosylation_sites`.
    """
    glygen_type = PTM_TYPE_TO_GLYGEN_TYPE.get(ptm_type)
    if glygen_type is None:
        raise ValueError(
            f"ptm_type '{ptm_type}' no es un tipo de glicosilacion soportado "
            f"(soportados: {sorted(PTM_TYPE_TO_GLYGEN_TYPE)})."
        )

    sites = fetch_glycosylation_sites(uniprot_accession)
    matches = [s for s in sites if s.get("type") == glygen_type and s.get("start_pos") == position]
    if not matches:
        return None

    return max(matches, key=lambda site: 1 if site.get("glytoucan_ac") else 0)
