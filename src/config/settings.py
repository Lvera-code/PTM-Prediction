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

    # --- Fase 3a: DeepPTMPred (motor 2/2 del consenso, Camino PDB unicamente) ---
    # Verificado leyendo github.com/kuikui-wang/DeepPTMPred directamente
    # (README.md, pred/train_PTM/predict.py, pred/train_PTM/e2_single_data.py,
    # pred/train_PTM/environment.yml) el 2026-07-27. Repo real, pesos SI
    # incluidos en el repo (.h5 por tipo, ~19MB c/u), pero SIN licencia
    # declarada (a diferencia de DeepMVP, GPL-3.0) -- verificar con Carlos
    # antes de cualquier uso mas alla de investigacion/TFG.
    #
    # HALLAZGO IMPORTANTE: a diferencia de DeepMVP, ni predict.py ni
    # e2_single_data.py tienen CLI real -- ptm_type/pdb_path/protein_id/ruta
    # del checkpoint ESM estan hardcodeados dentro de su bloque
    # `if __name__ == "__main__":`. Por eso este proyecto NO invoca esos
    # scripts directamente por subprocess (a diferencia de DeepMVP): usa un
    # runner propio (`src/engines/_deepptmpred_runner.py`) que importa sus
    # clases parametrizables (PredictConfig, PTMPredictor) y REIMPLEMENTA la
    # extraccion de features ESM-2 -- e2_single_data.py::extract_full_sequence_esm
    # tiene un bug real confirmado: redefine 'custom_checkpoint_path' como
    # variable LOCAL con una ruta absoluta hardcodeada de AutoDL
    # (/root/autodl-tmp/...), ignorando cualquier valor pasado o de modulo.
    DEEPPTMPRED_HOME: Path = Path(_env_str("DEEPPTMPRED_HOME", "DeepPTMPred"))
    DEEPPTMPRED_TRAIN_PTM_DIR: Path = Path(
        _env_str("DEEPPTMPRED_TRAIN_PTM_DIR", "DeepPTMPred/pred/train_PTM")
    )
    DEEPPTMPRED_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "DEEPPTMPRED_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_deepptmpred_runner.py"),
        )
    )
    # Python 3.10 + TensorFlow==2.15 + PyTorch 2.0 + PyRosetta + fair-esm (ver
    # environment.yml del repo): stack pesado y DISTINTO del de DeepMVP
    # (Python 3.7 + TF 2.4), requiere su propio venv/conda dedicado, nunca
    # compartido con DEEPMVP_PYTHON_BIN. tensorflow-addons (usado por
    # predict.py para su loss function) esta archivado/deprecado por Google
    # desde 2024 -- riesgo real de instalacion contra TF 2.15 sin verificar
    # todavia, no asumido como resuelto.
    DEEPPTMPRED_PYTHON_BIN: str = _env_str("DEEPPTMPRED_PYTHON_BIN", sys.executable)
    # Checkpoint ESM-2 (fair-esm, ~2.5GB), descarga separada (no incluido en
    # el repo, a diferencia de los pesos .h5 de PTM que si vienen incluidos).
    DEEPPTMPRED_ESM_CHECKPOINT: Path = Path(
        _env_str("DEEPPTMPRED_ESM_CHECKPOINT", "DeepPTMPred/esm/checkpoints/esm2_t33_650M_UR50D.pt")
    )
    # Cache de features ESM por accession (.npz), reutilizado entre corridas
    # y entre tipos de PTM del mismo accession (evita recalcular el embedding
    # completo 17 veces, una por tipo de PTM).
    DEEPPTMPRED_CUSTOM_ESM_DIR: Path = Path(
        _env_str("DEEPPTMPRED_CUSTOM_ESM_DIR", "DeepPTMPred/pred/custom_esm")
    )
    DEEPPTMPRED_TIMEOUT_SECONDS: int = _env_int("DEEPPTMPRED_TIMEOUT_SECONDS", 3600)
    # Los 17 tipos de PTM soportados, confirmado leyendo
    # predict.py::PredictConfig.ptm_aa_map (nombres exactos del repo, no
    # traducidos): a diferencia de DeepMVP (task=2 predice todos los tipos en
    # una sola invocacion), DeepPTMPred predice UN tipo por invocacion -- el
    # engine debe invocar el runner una vez POR TIPO, por accession.
    # Filtro de confianza (Fase 3, nucleo): a diferencia de DeepMVP,
    # DeepPTMPred no expone ningun mecanismo de calibracion tipo 'fpr' (ni en
    # predict.py ni en su README) -- solo una probabilidad cruda 0-1. Su
    # propio script usa un cutoff de 0.5 hardcodeado sin documentar de donde
    # sale (no calibrado contra ningun validation set publicado, a
    # diferencia del 'fpr' de DeepMVP). Este default de 0.5 es PROVISIONAL,
    # heredado sin verificacion adicional -- revisar si aparece evidencia de
    # un umbral mejor calibrado.
    DEEPPTMPRED_MIN_PROBABILITY: float = _env_float("DEEPPTMPRED_MIN_PROBABILITY", 0.5)

    DEEPPTMPRED_PTM_TYPES: tuple = (
        "phosphorylation", "acetylation", "ubiquitination", "hydroxylation",
        "gamma_carboxyglutamic_acid", "lys_methylation", "malonylation",
        "arg_methylation", "crotonylation", "succinylation", "glutathionylation",
        "sumoylation", "s_nitrosylation", "glutarylation", "citrullination",
        "o_linked_glycosylation", "n_linked_glycosylation",
    )

    @classmethod
    def ensure_dirs(cls) -> None:
        """Crea todas las carpetas de entrada/salida configuradas si no existen."""
        cls.FASTA_INPUT_DIR.mkdir(parents=True, exist_ok=True)
        cls.FASTA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
