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
    # CALIBRADO 2026-08-01 (reemplaza la calibracion 2026-07-31, invalidada
    # por el bug de plDDT -- ver abajo): umbrales reales por tipo, calculados
    # por `scripts/generate_deepptmpred_calibration.py` +
    # `_rebuild_calibration_summary.py` sobre
    # `github.com/meilerlab/PTMPrediction/data/ptm_data.csv.gz`
    # (n=75 proteinas/clase), CON dos fixes de train/inferencia ya aplicados
    # en `_deepptmpred_runner.py::_load_predict_module`: (1) phi/psi=0.0
    # forzado en inferencia (2026-07-31, igualaba la distribucion de
    # entrenamiento) y (2) local_plDDT recalculado desde el B-factor CA real
    # del pose de PyRosetta en vez del array `plDDT_values` de `predict.py`
    # (mal indexado por atomo en vez de por residuo, 2026-08-01) -- ver
    # STATUS.md, secciones "Calibracion real de DeepPTMPred + bug de
    # distribucion phi/psi" e "investigacion de n_linked_glycosylation y los
    # 4 tipos mediocres". Regla real (bug-fixed respecto al paper, no
    # literal): FPR<=0.18 TPR-max para phosphorylation, FPR<=0.20 TPR-max
    # para el resto; fallback a 0.5 si ningun punto de operacion cumple el
    # FPR maximo en la muestra.
    # Fuente exacta: `DeepPTMPred/data/calibration/summary.tsv` (gitignored,
    # regenerable con los dos scripts de arriba). Con el fix de plDDT, 16/17
    # tipos quedan en rango solido-a-excelente (AUROC 0.77-0.99) -- incluye
    # los 4 tipos que antes eran "mediocres" (crotonylation 0.61->0.82,
    # citrullination 0.66->0.78, glutarylation 0.67->0.77, s_nitrosylation
    # 0.68->0.77), confirmando que el bug de plDDT (no un techo del modelo)
    # era la causa real. `n_linked_glycosylation` sigue AUROC 0.51 (peor que
    # azar, limite real de los pesos publicados -- ver STATUS.md) -- por eso
    # esta EXCLUIDO del consenso en `ptm_annotation.py` pese a tener un
    # umbral calibrado aqui (el umbral no es el problema, es el poder
    # discriminativo del modelo).
    DEEPPTMPRED_CALIBRATED_THRESHOLDS: dict = {
        "acetylation": 0.6350621,
        "arg_methylation": 0.34068727,
        "citrullination": 0.36854228,
        "crotonylation": 0.86312497,
        "gamma_carboxyglutamic_acid": 0.24807824,
        "glutarylation": 0.4747097,
        "glutathionylation": 0.466462,
        "hydroxylation": 0.35899624,
        "lys_methylation": 0.43064716,
        "malonylation": 0.41699925,
        "n_linked_glycosylation": 0.99802846,
        "o_linked_glycosylation": 0.2619363,
        "phosphorylation": 0.24020174,
        "s_nitrosylation": 0.5140331,
        "succinylation": 0.50403893,
        "sumoylation": 0.37326753,
        "ubiquitination": 0.5321394,
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

    # --- Corroboracion opcional: MeToken (Camino PDB unicamente, NUNCA un motor
    # de consenso -- Decision 2 sigue pausada) ---
    # Verificado leyendo github.com/A4Bio/MeToken (ICLR 2025) el 2026-08-01: el
    # checkpoint publicado es un CLASIFICADOR DE TIPO en sitios YA CONOCIDOS
    # (confirmado en model_interface.py:40, la clase "no-PTM" queda enmascarada
    # cuando with_null_ptm=0, como viene el checkpoint publicado -- y
    # verificado empiricamente: predice tipos con alta confianza tambien en
    # posiciones sin PTM real). Por eso NUNCA decide pasa_umbral/consenso --
    # mismo patron puramente informativo que GlyGen
    # (src/structural/glygen_client.py): solo corrobora el TIPO en sitios que
    # el consenso YA acepto (pasa_umbral=true), ver
    # src/engines/ptm_annotation.py::annotate_pdb_path y
    # src/engines/metoken_engine.py.
    METOKEN_ENABLED: bool = _env_bool("METOKEN_ENABLED", True)
    METOKEN_HOME: Path = Path(_env_str("METOKEN_HOME", "MeToken"))
    METOKEN_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "METOKEN_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_metoken_runner.py"),
        )
    )
    # Python 3.10 + torch (CPU) + torch_scatter (compilado desde fuente, sin
    # wheel prebuilt para esta combinacion torch/CPU) + transformers +
    # biopython>=1.80 (ver bug 1 parcheado en _metoken_runner.py) -- venv
    # dedicado, nunca compartido con DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN.
    METOKEN_PYTHON_BIN: str = _env_str("METOKEN_PYTHON_BIN", sys.executable)
    # Pesos reales (release 1.0, pretrained_model.zip, ~88MB, verificado
    # descargable el 2026-08-01): checkpoint.ckpt es un state_dict plano
    # (confirmado con torch.load + model.load_state_dict -> "<All keys
    # matched successfully>"), no un checkpoint completo de pytorch-lightning
    # (ese es lightning_checkpoint.ckpt, sin usar aqui).
    METOKEN_CHECKPOINT: Path = Path(
        _env_str("METOKEN_CHECKPOINT", "MeToken/pretrained_model/checkpoint.ckpt")
    )
    METOKEN_TIMEOUT_SECONDS: int = _env_int("METOKEN_TIMEOUT_SECONDS", 600)

    # --- Corroboracion opcional: StackGlyEmbed (N-glicosilacion, ambos caminos --
    # FASTA y PDB, NUNCA un motor de consenso) ---
    # Motivo (decision 2026-08-01, ver STATUS.md): 'n_linked_glycosylation' en
    # DeepPTMPred esta CONFIRMADO como modelo muerto (AUROC ~0.51, ya excluido
    # del consenso via CONSENSUS_EXCLUDED_TYPES en ptm_annotation.py) -- no
    # arreglable reentrenando, un problema real del dataset/modelo publicado.
    # StackGlyEmbed (github.com/GaryChan-lab/StackGlyEmbed) es un tercer motor
    # INDEPENDIENTE de arquitectura (ProteinBERT + ESM-2 650M + ProtT5 apilados
    # -> meta-clasificador SVM), especializado solo en N-glicosilacion -- ya
    # instalado y verificado funcionando de verdad en el proyecto HERMANO
    # 'B-Cell-Epitope-Prediction' (proyecto 1, independiente de este por
    # decision explicita 2026-07-26 -- NO se importa codigo de un proyecto al
    # otro, pero SI se reusa su venv/pickles ya instalados como recurso
    # externo, mismo criterio que cualquier otro motor externo de este
    # proyecto). A diferencia de MeToken (requiere PDB, Camino PDB unicamente),
    # StackGlyEmbed solo necesita la secuencia completa -- aplica a ambos
    # caminos, ver src/engines/ptm_annotation.py.
    STACKGLYEMBED_ENABLED: bool = _env_bool("STACKGLYEMBED_ENABLED", True)
    # Venv REAL del proyecto hermano (torch/tensorflow/transformers/sklearn +
    # ProteinBERT ya instalados, ~pesado -- nunca reinstalado aqui). Verificado
    # 2026-08-01: existe de verdad en esta maquina (python3.10).
    STACKGLYEMBED_PYTHON_BIN: str = _env_str(
        "STACKGLYEMBED_PYTHON_BIN",
        "/home/enzo/DiffSBDD/B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed/bin/python",
    )
    # Runner PROPIO de este proyecto (vendorizado, adaptacion de la logica real
    # ya verificada en stackglyembed_predict_local.py del proyecto hermano --
    # NUNCA se importa ese archivo directamente). Vive dentro de este repo,
    # mismo patron que DEEPPTMPRED_RUNNER_SCRIPT/METOKEN_RUNNER_SCRIPT.
    STACKGLYEMBED_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "STACKGLYEMBED_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_stackglyembed_runner.py"),
        )
    )
    # Carpeta 'prediction/' del clon externo de StackGlyEmbed (proyecto
    # hermano): aqui viven los pickles del clasificador ya entrenado
    # (power_transformer_*.sav, base_layer_pickle_files/). Verificado 2026-08-01
    # con 'ls' real, no asumido.
    STACKGLYEMBED_MODELS_DIR: str = _env_str(
        "STACKGLYEMBED_MODELS_DIR",
        "/home/enzo/DiffSBDD/B-Cell-Epitope-Prediction/StackGlyEmbed/prediction",
    )
    # Pesos de ProtT5 (~3GB) reusados de scipion-chem-tmbed (mismo encoder que
    # usa el proyecto hermano para su propio StackGlyEmbed, ver su
    # settings.py) -- nunca descargados de nuevo aqui.
    STACKGLYEMBED_T5_MODEL_PATH: str = _env_str(
        "STACKGLYEMBED_T5_MODEL_PATH",
        "/home/enzo/DiffSBDD/scipion-chem-tmbed/tmbed_src/tmbed/models/t5",
    )
    STACKGLYEMBED_ESM_MODEL_NAME: str = _env_str("STACKGLYEMBED_ESM_MODEL_NAME", "facebook/esm2_t33_650M_UR50D")
    # Generoso por defecto: carga en frio de 3 modelos (ProteinBERT + ESM-2
    # 650M + ProtT5) sobre CPU antes de procesar el primer sitio -- mismo valor
    # que el proyecto hermano usa para la misma carga en frio.
    STACKGLYEMBED_TIMEOUT_SECONDS: int = _env_int("STACKGLYEMBED_TIMEOUT_SECONDS", 900)

    # --- Fase A / Extension 3: modelado estructural real (PyRosetta), Camino PDB
    # unicamente -- conectado al pipeline principal 2026-08-03 ---
    # Revierte la decision 2026-07-27 de que D (apply_workflow_filter) no ruta a
    # Extension 3/Fase A porque esas fases no existian todavia -- ya existen y
    # estan verificadas con corridas reales (ver STATUS.md). Cubre solo 9/17
    # tipos de PTM (ver src/structural/fase_a_dispatch.py::SUPPORTED_PTM_TYPES);
    # los otros 8 no tienen modulo de Fase A (requeririan construir un residuo
    # no-canonico propio, cheminformatica real, no implementado).
    #
    # Modelar TODOS los sitios que pasan el umbral es computacionalmente
    # inviable (un caso real como Tau acepta ~572 sitios; cada modelado tarda de
    # minutos, no segundos) -- FASE_A_TOP_N_PER_TYPE acota el barrido al sitio
    # de mayor score por cada uno de los 9 tipos soportados (representativo,
    # demostrable en una corrida real), no un recorte arbitrario.
    FASE_A_ENABLED: bool = _env_bool("FASE_A_ENABLED", True)
    FASE_A_TOP_N_PER_TYPE: int = _env_int("FASE_A_TOP_N_PER_TYPE", 1)
    # Constante de datos pura (sin pyrosetta) para que src/engines/ptm_annotation.py
    # (proceso principal, sin PyRosetta instalado) pueda seleccionar candidatos
    # sin importar src/structural/*. Fuente de verdad real: la union de los 3
    # SUPPORTED_PTM_TYPES de pyrosetta_ptm_patch/pyrosetta_glycan_patch/
    # ubiquitin_sumo (ver src/structural/fase_a_dispatch.py, que valida en tiempo
    # de import que esta lista coincide exactamente -- si alguna vez difieren,
    # falla alto en vez de seleccionar candidatos para un tipo sin soporte real).
    FASE_A_SUPPORTED_PTM_TYPES: tuple = (
        "phosphorylation", "acetylation", "hydroxylation", "gamma_carboxyglutamic_acid",
        "lys_methylation", "n_linked_glycosylation", "o_linked_glycosylation",
        "ubiquitination", "sumoylation",
    )
    # Plantilla (sin pyrosetta) de todas las claves que puede tener un
    # resultado de Fase A -- unica fuente de verdad compartida entre
    # fase_a_engine.py (proceso principal), fase_a_dispatch.py y
    # _fase_a_runner.py (subprocess con PyRosetta). Evita que las distintas
    # rutas de salida (modelado/sin_soporte/error/no_disponible) queden
    # desalineadas entre si -- el bug real 2026-08-03 de "ddg_std" calculado
    # pero nunca incluido en el dict de resultado (se perdia antes de llegar
    # al reporte final) fue exactamente este tipo de desalineacion.
    FASE_A_RESULT_TEMPLATE: dict = {
        "estado": None, "clase": None, "ddg": None, "ddg_std": None,
        "wt_score": None, "wt_score_std": None, "mut_score": None, "mut_score_std": None,
        "glycan_tree": None, "glygen_evidencia": None, "conjugation_metrics": None,
        "output_pdb": None, "error": None,
    }
    FASE_A_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "FASE_A_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_fase_a_runner.py"),
        )
    )
    # Requiere PyRosetta -- reusa el MISMO conda env 'deepptmpred' ya instalado
    # y verificado para DeepPTMPred (ver DEEPPTMPRED_PYTHON_BIN arriba), nunca
    # un env separado: ambos necesitan exactamente el mismo PyRosetta real.
    FASE_A_PYTHON_BIN: str = _env_str("FASE_A_PYTHON_BIN", sys.executable)
    # Generoso por defecto: clase 1 (ddG) corre nstruct*2 relax independientes
    # (~1 min cada uno sobre una proteina de tamano medio, ver STATUS.md) mas
    # un relax final para el PDB de salida -- varios minutos en el peor caso.
    FASE_A_TIMEOUT_SECONDS: int = _env_int("FASE_A_TIMEOUT_SECONDS", 900)

    @classmethod
    def ensure_dirs(cls) -> None:
        """Crea todas las carpetas de entrada/salida configuradas si no existen."""
        cls.FASTA_INPUT_DIR.mkdir(parents=True, exist_ok=True)
        cls.FASTA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
