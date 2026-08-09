"""Lookup de ground truth PTM contra la base derivada de dbPTM (decision 2026-08-09).

## Por que existe

`biological_panel.py::PANEL` solo cubre las ~20-30 proteinas curadas que ademas
corren el pipeline completo para medir recall (ver `scripts/select_recall_subset.py`).
Este modulo consulta la base GRANDE (`scripts/import_dbptm_panel.py`, miles de
proteinas humanas reviewed) para poblar la columna "Literatura" del CLI en
CUALQUIER proteina humana -- puramente informativo, nunca decide
`pasa_umbral`/consenso, mismo patron ya establecido por
`src/structural/glygen_client.py`/`uniprot_localization_client.py`.
"""

import json
import sqlite3
from pathlib import Path
from typing import Dict, Tuple

from src.config.settings import Settings

DBPTM_LOOKUP_PATH: Path = Settings.DBPTM_DATA_DIR / "lookup.sqlite3"


def lookup_ground_truth(
    accession: str, db_path: Path = DBPTM_LOOKUP_PATH,
) -> Dict[Tuple[int, str], Tuple[str, Tuple[int, ...]]]:
    """``{(posicion, tipo_ptm): (tier, pmids)}`` para ``accession``, o ``{}`` si no hay datos.

    ``pmids`` es la tupla cruda de dbPTM para ese sitio -- a diferencia de
    ``biological_panel.py::PANEL``, estos PMIDs NUNCA se verificaron
    individualmente contra NCBI eutils (ver decision 2026-08-09: la
    verificacion real del subconjunto curado encontro PMIDs inexistentes o
    mal atribuidos en el borrador dbPTM). El llamador debe tratar esto como
    una señal mas debil que la del panel curado.

    Nunca lanza: ``db_path`` puede no existir todavia (setup opcional de
    ``scripts/import_dbptm_panel.py`` no corrido) -- se trata igual que "sin
    datos", mismo criterio no-decisorio que el resto de clientes de
    ``src/validation``/``src/structural``.
    """
    if not db_path.is_file():
        return {}
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT position, ptm_type, tier, pmids FROM dbptm_sites WHERE accession = ?",
            (accession,),
        ).fetchall()
    finally:
        conn.close()
    return {
        (position, ptm_type): (tier, tuple(json.loads(pmids_json)))
        for position, ptm_type, tier, pmids_json in rows
    }
