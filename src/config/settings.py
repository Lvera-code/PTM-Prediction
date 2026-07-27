"""Configuracion centralizada del pipeline: rutas, umbrales y credenciales externas.

Todos los parametros ajustables se resuelven desde variables de entorno con
valores por defecto conservadores, para permitir reconfiguracion sin tocar
codigo fuente ni comprometer credenciales al subir el repositorio a GitHub.
"""

import os
import sys
from pathlib import Path


def _env_str(key: str, default: str) -> str:
    return os.environ.get(key, default)


def _env_float(key: str, default: float) -> float:
    raw = os.environ.get(key)
    if raw is None:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def _env_int(key: str, default: int) -> int:
    raw = os.environ.get(key)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _env_bool(key: str, default: bool) -> bool:
    raw = os.environ.get(key)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


class Settings:
    """Punto unico de verdad para toda configuracion del pipeline."""

    # --- Fase 1 / Orquestador: carpetas de entrada y salida del pipeline ---
    FASTA_INPUT_DIR: Path = Path(_env_str("FASTA_INPUT_DIR", "fasta_inputs"))
    FASTA_OUTPUT_DIR: Path = Path(_env_str("FASTA_OUTPUT_DIR", "fasta_outputs"))

    # --- Fase 1.5: Extraccion de estructura (PDB/mmCIF via gemmi, LOCAL) ---
    # Estrategia de seleccion de cadena cuando el archivo de entrada tiene mas
    # de una cadena proteica (ver `src/utils/structure_parser.py`). Nunca
    # implicita: la cadena elegida siempre se loggea con su motivo.
    #   'longest'  -> se elige la cadena con mas residuos en su polimero.
    #   'explicit' -> se usa PDB_EXPLICIT_CHAIN_ID (obligatorio en ese caso).
    PDB_CHAIN_SELECTION_STRATEGY: str = _env_str("PDB_CHAIN_SELECTION_STRATEGY", "longest")
    PDB_EXPLICIT_CHAIN_ID: str = _env_str("PDB_EXPLICIT_CHAIN_ID", "")

    # --- Fase 3a: DeepMVP (motor unico Camino FASTA, motor 1/2 Camino PDB) ---
    # Verificado leyendo github.com/bzhanglab/DeepMVP directamente (README.md,
    # DeepMVP.py, lib/PTModels.py, lib/Metrics.py) el 2026-07-27. Repo real,
    # empaquetado (requirements.txt + environment.yml), GPL-3.0.
    DEEPMVP_HOME: Path = Path(_env_str("DEEPMVP_HOME", "DeepMVP"))
    DEEPMVP_SCRIPT_NAME: str = _env_str("DEEPMVP_SCRIPT_NAME", "DeepMVP.py")
    # Python 3.7.10 + TensorFlow==2.4.2 (ver environment.yml del repo): stack
    # antiguo, se recomienda venv dedicado (mismo patron que BepiPred-3.0 en
    # BCell-Epitope-Prediction). Por defecto usa el mismo interprete que
    # corre pipeline.py; apunta esto a un venv dedicado en produccion.
    DEEPMVP_PYTHON_BIN: str = _env_str("DEEPMVP_PYTHON_BIN", sys.executable)
    # Pesos pre-entrenados NO incluidos en el repo (confirmado en README):
    # descarga manual desde DEEPMVP_DOWNLOAD_URL (un .tar.gz por version),
    # descomprimir dentro de esta carpeta. Contiene un subdirectorio por tipo
    # de PTM (acetylation_k/, glycosylation_n/, methylation_k/, methylation_r/,
    # phosphorylation_st/, phosphorylation_y/, sumoylation_k/, ubiquitination_k/),
    # cada uno con su propio model.json + ensemble de modelos .h5 -- confirmado
    # leyendo lib/PTModels.py::ptm_prediction_for_multiple_ptms.
    DEEPMVP_MODEL_DIR: Path = Path(_env_str("DEEPMVP_MODEL_DIR", "DeepMVP/models"))
    DEEPMVP_DOWNLOAD_URL: str = "https://deepmvp.ptmax.org/"
    DEEPMVP_TIMEOUT_SECONDS: int = _env_int("DEEPMVP_TIMEOUT_SECONDS", 1800)
    # Tarea 2 = "PTM site prediction" en el CLI de DeepMVP (tarea 1 es
    # prediccion de impacto de mutacion, no usada aqui). Confirmado leyendo
    # DeepMVP.py directamente, no configurable por Settings: es un valor fijo
    # del contrato con el CLI real, no un parametro ajustable del pipeline.
    DEEPMVP_TASK_PTM_SITE_PREDICTION: str = "2"
    # Filtro de confianza (Fase 3, nucleo -- decision 2026-07-27): DeepMVP no
    # publica un cutoff de probabilidad fijo por tipo de PTM. En su lugar,
    # cada carpeta de pesos trae su propio 'site_prediction.tsv' de
    # validacion, y la columna 'fpr' que devuelve por fila es el FPR real que
    # tendria ese modelo especifico si se usara esa probabilidad como corte
    # (confirmado leyendo lib/Metrics.py::add_confidence_metrics/calc_FPR).
    # Es un umbral MEJOR calibrado que un score de probabilidad fijo: el
    # filtro de Fase 3 usa 'fpr <= DEEPMVP_MAX_FPR', no un umbral de
    # 'y_pred'. El score crudo (y_pred) se conserva siempre igual.
    DEEPMVP_MAX_FPR: float = _env_float("DEEPMVP_MAX_FPR", 0.05)

    @classmethod
    def ensure_dirs(cls) -> None:
        """Crea todas las carpetas de entrada/salida configuradas si no existen."""
        cls.FASTA_INPUT_DIR.mkdir(parents=True, exist_ok=True)
        cls.FASTA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
