"""Corroboracion opcional de N-glicosilacion (StackGlyEmbed), ambos caminos (FASTA y PDB).

NUNCA es un motor de consenso -- StackGlyEmbed no decide ``pasa_umbral``/
``consenso``, solo corrobora candidatos de ``n_linked_glycosylation`` que ya
propuso DeepMVP y/o DeepPTMPred. Este modulo invoca
``src/engines/_stackglyembed_runner.py`` (vendorizado, adaptacion propia de
la logica real ya verificada en el proyecto hermano
``B-Cell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py`` --
ver el docstring del runner) via subprocess, sobre el venv dedicado del
proyecto hermano (``Settings.STACKGLYEMBED_PYTHON_BIN``) -- reusado como
recurso externo, nunca importado directamente (independencia entre
proyectos, decision 2026-07-26).

Motivo (decision 2026-08-01, ver STATUS.md): ``n_linked_glycosylation`` en
DeepPTMPred esta CONFIRMADO como modelo muerto (AUROC ~0.51, ya excluido del
consenso). StackGlyEmbed es un tercer motor INDEPENDIENTE de arquitectura
(ProteinBERT + ESM-2 650M + ProtT5 -> meta-clasificador SVM), especializado
solo en N-glicosilacion -- util para corroborar esos candidatos sin
pretender arreglar el modelo muerto.

A diferencia de MeToken (requiere un PDB, Camino PDB unicamente), StackGlyEmbed
solo necesita la SECUENCIA COMPLETA de la proteina -- aplica igual a Camino
FASTA (``annotate_fasta_path``) y Camino PDB (``annotate_pdb_path``).

Mismo patron no-decisorio que GlyGen/MeToken: cualquier fallo (venv/pickles
no instalados, subproceso que revienta, timeout, salida malformada) degrada
a un aviso en el log y devuelve ``{}`` -- NUNCA propaga una excepcion que
tumbe el pipeline principal. Quien llama
(``src/engines/ptm_annotation.py``) decide que hacer con un resultado vacio
(no anota corroboracion en ninguna fila, el resto del pipeline sigue
exactamente igual).
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd

from src.config.settings import Settings
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

OUTPUT_COLUMNS = ["position", "stackglyembed_veredicto", "stackglyembed_score"]


def get_nglyco_corroboration(
    sequence: str,
    positions: Sequence[int],
    result_dir: Path = None,
    filename_prefix: str = "",
) -> Dict[int, Dict]:
    """Corre StackGlyEmbed sobre ``sequence`` para ``positions`` y devuelve la corroboracion por posicion.

    Args:
        sequence: Secuencia COMPLETA de la proteina (Fase 1 saneada o Fase
            1.5 ATMSEQ) -- NUNCA un fragmento/peptido, ver docstring del
            runner para por que la posicion 1-based coincide exactamente sin
            conversion de offsets.
        positions: Posiciones 1-based de la Asparagina de cada secuon
            N-X-[S/T] candidato a corroborar (tipicamente las filas con
            ``tipo_ptm`` en {``n_linked_glycosylation``, ``glycosylation_n``}
            y ``pasa_umbral=True`` -- ver ``ptm_annotation.py``).
        result_dir: Carpeta donde guardar el CSV crudo. Si es ``None``, usa
            ``Settings.FASTA_OUTPUT_DIR``.
        filename_prefix: Prefijo (tipicamente ``f"{accession}_"``) para el
            archivo crudo persistido en ``result_dir``.

    Returns:
        Dict ``{posicion: {"stackglyembed_veredicto": str, "stackglyembed_score": float}}``,
        una entrada por posicion que StackGlyEmbed pudo evaluar (puede
        faltar alguna si cae fuera de ``sequence``, ver runner). Dict vacio
        si StackGlyEmbed no esta instalado, el subproceso falla, excede el
        timeout, o la lista de posiciones esta vacia -- NUNCA lanza.
    """
    if not positions:
        return {}

    if not _validate_installation():
        return {}

    base_output_dir = Path(result_dir) if result_dir else Settings.FASTA_OUTPUT_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = base_output_dir / f"{filename_prefix}stackglyembed.csv"

    cmd = [
        Settings.STACKGLYEMBED_PYTHON_BIN, str(Settings.STACKGLYEMBED_RUNNER_SCRIPT.resolve()),
        "--sequence", sequence,
        "--positions", *[str(p) for p in positions],
        "--models-dir", str(Settings.STACKGLYEMBED_MODELS_DIR),
        "--t5-model-path", str(Settings.STACKGLYEMBED_T5_MODEL_PATH),
        "--esm-model-name", Settings.STACKGLYEMBED_ESM_MODEL_NAME,
        "--out-csv", str(out_csv.resolve()),
    ]

    logger.info(
        "Ejecutando StackGlyEmbed (corroboracion de N-glicosilacion) sobre %d posicion(es).",
        len(positions),
    )
    try:
        # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
        # backend inline que no existe en el conda env aislado) -- ver
        # deepptmpred_engine.py para el caso real que disparo esto.
        env = {**os.environ, "MPLBACKEND": "Agg"}
        subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=Settings.STACKGLYEMBED_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "StackGlyEmbed termino con exit code %d -- se omite la corroboracion de N-glicosilacion "
            "(no fatal, no afecta consenso/pasa_umbral). stderr: %s",
            exc.returncode, (exc.stderr or "<vacio>")[-6000:],
        )
        return {}
    except subprocess.TimeoutExpired:
        logger.warning(
            "StackGlyEmbed excedio el tiempo limite de %ds -- se omite la corroboracion de "
            "N-glicosilacion (no fatal).", Settings.STACKGLYEMBED_TIMEOUT_SECONDS,
        )
        return {}
    except FileNotFoundError as exc:
        logger.warning(
            "No se pudo invocar STACKGLYEMBED_PYTHON_BIN='%s' -- se omite la corroboracion de "
            "N-glicosilacion (no fatal): %s", Settings.STACKGLYEMBED_PYTHON_BIN, exc,
        )
        return {}

    if not out_csv.is_file():
        logger.warning(
            "StackGlyEmbed reporto exit code 0 pero no genero '%s' -- se omite la corroboracion "
            "de N-glicosilacion (no fatal).", out_csv,
        )
        return {}

    try:
        df = pd.read_csv(out_csv)
    except Exception as exc:  # noqa: BLE001 -- degradacion no fatal, cualquier fallo de lectura cuenta
        logger.warning(
            "'%s' no se pudo leer como CSV valido -- se omite la corroboracion de N-glicosilacion "
            "(no fatal): %s", out_csv, exc,
        )
        return {}

    missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
    if missing:
        logger.warning(
            "'%s' no tiene las columnas esperadas %s (faltan: %s) -- se omite la corroboracion de "
            "N-glicosilacion (no fatal).", out_csv, OUTPUT_COLUMNS, missing,
        )
        return {}

    return {
        int(r["position"]): {
            "stackglyembed_veredicto": r["stackglyembed_veredicto"],
            "stackglyembed_score": float(r["stackglyembed_score"]),
        }
        for _, r in df.iterrows()
    }


def _validate_installation() -> bool:
    """Comprueba que el venv/script/pickles del recurso externo esten presentes.

    A diferencia de ``DeepMVPEngine``/``DeepPTMPredEngine`` (motores
    OBLIGATORIOS del consenso, su ausencia lanza una excepcion fatal
    accionable), aqui una instalacion faltante NUNCA es fatal --
    StackGlyEmbed es puramente informativo. Se registra un aviso una unica
    vez por llamada con instrucciones de instalacion.
    """
    python_bin = Path(Settings.STACKGLYEMBED_PYTHON_BIN)
    if not python_bin.is_file():
        logger.warning(
            "No se encontro el interprete Python de StackGlyEmbed en '%s' -- se omite la "
            "corroboracion de N-glicosilacion (no fatal). Requiere el venv ya instalado en el "
            "proyecto hermano 'B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed' (ver "
            "README.md) o apunta STACKGLYEMBED_PYTHON_BIN a otra instalacion.",
            python_bin,
        )
        return False
    if not Settings.STACKGLYEMBED_RUNNER_SCRIPT.is_file():
        logger.warning(
            "No se encontro el runner de StackGlyEmbed en '%s' (deberia venir con el repo, esto "
            "no deberia pasar) -- se omite la corroboracion de N-glicosilacion (no fatal).",
            Settings.STACKGLYEMBED_RUNNER_SCRIPT,
        )
        return False
    models_dir = Path(Settings.STACKGLYEMBED_MODELS_DIR)
    if not (models_dir / "base_layer_pickle_files" / "SVM_meta_layer.sav").is_file():
        logger.warning(
            "No se encontraron los pickles del clasificador de StackGlyEmbed en '%s' -- se omite "
            "la corroboracion de N-glicosilacion (no fatal). Este clon externo vive en el proyecto "
            "hermano 'B-Cell-Epitope-Prediction/StackGlyEmbed/prediction' (ver README.md), apunta "
            "STACKGLYEMBED_MODELS_DIR a su ubicacion si esta en otro sitio.",
            models_dir,
        )
        return False
    return True
