"""PTM ground truth lookup against the dbPTM-derived database (decision 2026-08-09).

## Why it exists

`biological_panel.py::PANEL` only covers the ~20-30 curated proteins that also
run the complete pipeline to measure recall (see `scripts/select_recall_subset.py`).
This module queries the LARGE database (`scripts/import_dbptm_panel.py`, thousands of
reviewed human proteins) in order to populate the CLI's "Literatura" column for
ANY human protein. It is purely informative and never decides
`pasa_umbral`/consensus, the same pattern already established by
`src/structural/uniprot_localization_client.py`.
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
