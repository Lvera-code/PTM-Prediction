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

    # --- Fase 2: DeepMVP (motor unico Camino FASTA, motor 1/2 Camino PDB) ---
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

    # --- Fase 2: DeepPTMPred (motor 2/2 del consenso, Camino PDB unicamente) ---
    # Verificado leyendo github.com/kuikui-wang/DeepPTMPred directamente
    # (README.md, pred/train_PTM/predict.py, pred/train_PTM/e2_single_data.py,
    # pred/train_PTM/environment.yml) el 2026-07-27. Repo real, pesos SI
    # incluidos en el repo (.h5 por tipo, ~19MB c/u). Licencia RESUELTA
    # 2026-07-29: el repo no tiene LICENSE declarado, pero Junwen Wang
    # (autor de correspondencia del paper, CC BY-NC 4.0) confirmo
    # directamente por email que el codigo del repo sigue esos mismos
    # terminos -- uso no comercial (investigacion/TFG, Scipion/CNB) cubierto.
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
    # desde 2024 -- RIESGO DESCARTADO 2026-07-27: instalado e importado de
    # verdad en un venv aislado contra TF 2.15, sin error (solo un
    # UserWarning de deprecacion esperado). Ver STATUS.md.
    DEEPPTMPRED_PYTHON_BIN: str = _env_str("DEEPPTMPRED_PYTHON_BIN", sys.executable)
    # Checkpoint ESM-2 (fair-esm, ~2.5GB), descarga separada (no incluido en
    # el repo, a diferencia de los pesos .h5 de PTM que si vienen incluidos).
    # IMPORTANTE (verificado 2026-07-27 leyendo esm/pretrained.py directamente):
    # fair-esm exige un archivo COMPANERO junto al checkpoint principal,
    # '<checkpoint>-contact-regression.pt' (mismo directorio, mismo nombre
    # base) -- su heuristica interna (_has_regression_weights) no excluye
    # los modelos esm2_*, asi que SIEMPRE lo busca salvo variantes esm1v/
    # esm_if/270K/500K. Si falta, torch.load() falla con FileNotFoundError
    # LOCAL (no es un intento de red -- ver docstring de
    # _deepptmpred_runner.py::_extract_esm_features para el detalle
    # completo de por que esto SI es 100% local pese a este archivo extra).
    # Verificar al descargar el checkpoint que este archivo compañero
    # tambien este disponible (mismo origen que el checkpoint principal).
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
    # diferencia del 'fpr' de DeepMVP).
    #
    # CALIBRADO 2026-07-31: umbrales reales por tipo, calculados por
    # `scripts/generate_deepptmpred_calibration.py` + `_rebuild_calibration_summary.py`
    # sobre `github.com/meilerlab/PTMPrediction/data/ptm_data.csv.gz`
    # (376,557 filas reales, n=75 proteinas/clase), CON el fix de
    # 2026-07-31 ya aplicado (phi/psi=0.0 forzado en inferencia para igualar
    # la distribucion de entrenamiento -- ver STATUS.md, seccion "Calibracion
    # real de DeepPTMPred + bug de distribucion phi/psi"). Regla real
    # (bug-fixed respecto al paper, no literal): FPR<=0.18 TPR-max para
    # phosphorylation, FPR<=0.20 TPR-max para el resto; fallback a 0.5 si
    # ningun punto de operacion cumple el FPR maximo en la muestra
    # (arg_methylation, lys_methylation cayeron en este caso).
    # Fuente exacta: `DeepPTMPred/data/calibration/summary.tsv` (gitignored,
    # regenerable con los dos scripts de arriba). AUROC de referencia y
    # salvedades por tipo en STATUS.md -- 14/17 tipos en rango solido-a-
    # excelente; `n_linked_glycosylation` (AUROC 0.49, peor que azar) y los
    # 4 tipos mediocres (s_nitrosylation/glutarylation/citrullination/
    # crotonylation, 0.61-0.68) usan igual su umbral calibrado aqui pese al
    # AUROC bajo -- son limitaciones de PODER DISCRIMINATIVO del modelo, no
    # de calibracion del umbral; investigacion de causa raiz pendiente.
    DEEPPTMPRED_CALIBRATED_THRESHOLDS: dict = {
        "acetylation": 0.6299973,
        "arg_methylation": 0.34319976,
        "citrullination": 0.39169243,
        "crotonylation": 0.93576247,
        "gamma_carboxyglutamic_acid": 0.21486656,
        "glutarylation": 0.4618006,
        "glutathionylation": 0.47251946,
        "hydroxylation": 0.3752605,
        "lys_methylation": 0.41040188,
        "malonylation": 0.4223403,
        "n_linked_glycosylation": 0.99741435,
        "o_linked_glycosylation": 0.3105993,
        "phosphorylation": 0.24348031,
        "s_nitrosylation": 0.500678,
        "succinylation": 0.48264492,
        "sumoylation": 0.4091335,
        "ubiquitination": 0.5434938,
    }
    # Fallback si DEEPPTMPRED_CALIBRATED_THRESHOLDS no tiene el tipo (no
    # deberia ocurrir con los 17 tipos de DEEPPTMPRED_PTM_TYPES, pero
    # mantiene el pipeline funcional si se agrega un tipo nuevo sin
    # recalibrar). Tambien override manual via env var para experimentacion
    # puntual sin editar codigo.
    DEEPPTMPRED_MIN_PROBABILITY: float = _env_float("DEEPPTMPRED_MIN_PROBABILITY", 0.5)

    @classmethod
    def deepptmpred_threshold_for(cls, ptm_type: str) -> float:
        """Umbral calibrado para ``ptm_type``, o el fallback generico si no existe."""
        return cls.DEEPPTMPRED_CALIBRATED_THRESHOLDS.get(ptm_type, cls.DEEPPTMPRED_MIN_PROBABILITY)

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
