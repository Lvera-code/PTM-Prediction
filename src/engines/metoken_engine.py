"""Corroboracion opcional de TIPO (MeToken), Camino PDB unicamente.

NUNCA es un motor de consenso -- Decision 2 (segundo motor completo del
Camino PDB) sigue pausada, ver STATUS.md. Este modulo invoca
``src/engines/_metoken_runner.py`` (vendorizado, ver su docstring para el
detalle completo de por que MeToken existe aqui, los 2 bugs reales
parcheados, y la deteccion de la clase null/rare) via subprocess, sobre el
venv dedicado ``Settings.METOKEN_PYTHON_BIN``.

Mismo patron no-decisorio que ``src/structural/uniprot_localization_client.py``:
cualquier fallo (repo no instalado, checkpoint ausente, subproceso que
revienta, timeout, salida malformada) degrada a un aviso en el log y
devuelve ``{}`` -- NUNCA
propaga una excepcion que tumbe el pipeline principal. Quien llama
(``src/engines/ptm_annotation.py::annotate_pdb_path``) decide que hacer con
un resultado vacio (no anota corroboracion en ninguna fila, el resto del
pipeline sigue exactamente igual).
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd

from src.config.settings import Settings
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

OUTPUT_COLUMNS = ["position", "metoken_type", "metoken_probability"]


def get_type_corroboration(
    pdb_path: Path,
    positions: Sequence[int],
    chain_id: str = "A",
    result_dir: Path = None,
) -> Dict[int, Dict]:
    """Corre MeToken sobre ``pdb_path`` para ``positions`` (1-based) y devuelve la corroboracion por posicion.

    Args:
        pdb_path: PDB de una sola cadena (mismo criterio que DeepPTMPred --
            se espera ``record.chain_pdb_path`` de Fase 1.5, para que las
            posiciones 1-based coincidan exactamente con ``record.sequence``,
            no ``record.pdb_path`` original que puede tener mas de una
            cadena).
        positions: Posiciones 1-based ya aceptadas por el consenso
            (``pasa_umbral=true``) a corroborar -- ver
            ``ptm_annotation.py::annotate_pdb_path``.
        chain_id: Cadena a leer del PDB (default ``"A"``, mismo default que
            ``inference.py::get_seq_str`` del repo).
        result_dir: Carpeta donde guardar el CSV crudo. Si es ``None``, usa
            ``Settings.FASTA_OUTPUT_DIR``.

    Returns:
        Dict ``{posicion: {"metoken_type": str, "metoken_probability": float}}``,
        una entrada por posicion que MeToken pudo predecir (puede faltar
        alguna si cae fuera de la secuencia que ``get_seq_str`` logro leer
        del PDB, ver runner). Dict vacio si MeToken no esta instalado, el
        subproceso falla, excede el timeout, o la lista de posiciones esta
        vacia -- NUNCA lanza.
    """
    if not positions:
        return {}

    if not _validate_installation():
        return {}

    base_output_dir = Path(result_dir) if result_dir else Settings.FASTA_OUTPUT_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = base_output_dir / f"{Path(pdb_path).stem}_metoken.csv"

    cmd = [
        Settings.METOKEN_PYTHON_BIN, str(Settings.METOKEN_RUNNER_SCRIPT.resolve()),
        "--repo-dir", str(Settings.METOKEN_HOME.resolve()),
        "--checkpoint-path", str(Settings.METOKEN_CHECKPOINT.resolve()),
        "--pdb-path", str(Path(pdb_path).resolve()),
        "--chain-id", chain_id,
        "--positions", *[str(p) for p in positions],
        "--out-csv", str(out_csv.resolve()),
    ]

    logger.info("Ejecutando MeToken (corroboracion de tipo) sobre '%s': %s", pdb_path, " ".join(cmd))
    try:
        # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
        # backend inline que no existe en el conda env aislado) -- ver
        # deepptmpred_engine.py para el caso real que disparo esto.
        env = {**os.environ, "MPLBACKEND": "Agg"}
        subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=Settings.METOKEN_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "MeToken termino con exit code %d para '%s' -- se omite la corroboracion de tipo "
            "(no fatal, no afecta consenso/pasa_umbral). stderr: %s",
            exc.returncode, pdb_path, (exc.stderr or "<vacio>")[-6000:],
        )
        return {}
    except subprocess.TimeoutExpired:
        logger.warning(
            "MeToken excedio el tiempo limite de %ds para '%s' -- se omite la corroboracion de tipo "
            "(no fatal).", Settings.METOKEN_TIMEOUT_SECONDS, pdb_path,
        )
        return {}
    except FileNotFoundError as exc:
        logger.warning(
            "No se pudo invocar METOKEN_PYTHON_BIN='%s' -- se omite la corroboracion de tipo "
            "(no fatal): %s", Settings.METOKEN_PYTHON_BIN, exc,
        )
        return {}

    if not out_csv.is_file():
        logger.warning(
            "MeToken reporto exit code 0 pero no genero '%s' -- se omite la corroboracion de tipo "
            "(no fatal).", out_csv,
        )
        return {}

    try:
        df = pd.read_csv(out_csv)
    except Exception as exc:  # noqa: BLE001 -- degradacion no fatal, cualquier fallo de lectura cuenta
        logger.warning(
            "'%s' no se pudo leer como CSV valido -- se omite la corroboracion de tipo (no fatal): %s",
            out_csv, exc,
        )
        return {}

    missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
    if missing:
        logger.warning(
            "'%s' no tiene las columnas esperadas %s (faltan: %s) -- se omite la corroboracion de "
            "tipo (no fatal).", out_csv, OUTPUT_COLUMNS, missing,
        )
        return {}

    return {
        int(r["position"]): {"metoken_type": r["metoken_type"], "metoken_probability": float(r["metoken_probability"])}
        for _, r in df.iterrows()
    }


def _validate_installation() -> bool:
    """Comprueba que el repo clonado, la runner script y el checkpoint esten presentes.

    A diferencia de ``DeepMVPEngine``/``DeepPTMPredEngine`` (motores
    OBLIGATORIOS del consenso, su ausencia lanza una excepcion fatal
    accionable), aqui una instalacion faltante NUNCA es fatal -- MeToken es
    puramente informativo. Se registra un aviso una unica vez por llamada
    con instrucciones de instalacion.
    """
    repo_predict_entry = Settings.METOKEN_HOME / "inference.py"
    if not repo_predict_entry.is_file():
        logger.warning(
            "MeToken no esta instalado en '%s' -- se omite la corroboracion de tipo (no fatal). "
            "Clona 'git clone https://github.com/A4Bio/MeToken' en la raiz del proyecto (o apunta "
            "METOKEN_HOME a su ubicacion) y descarga los pesos (release 1.0, pretrained_model.zip) "
            "para habilitarla. Ver README.md.",
            Settings.METOKEN_HOME,
        )
        return False
    if not Settings.METOKEN_CHECKPOINT.is_file():
        logger.warning(
            "No se encontro el checkpoint de MeToken en '%s' -- se omite la corroboracion de tipo "
            "(no fatal). Descarga 'pretrained_model.zip' desde "
            "https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip y "
            "descomprimelo en METOKEN_HOME.",
            Settings.METOKEN_CHECKPOINT,
        )
        return False
    return True
