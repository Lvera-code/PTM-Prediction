"""Motor real de consenso para 'n_linked_glycosylation' (EMNGly), Camino PDB unicamente.

Reemplaza a CoNglyPred (candidato original de Decision 2, confirmado
DEFINITIVAMENTE sin pesos publicados en ningun sitio -- ver STATUS.md,
seccion Decision 2). Este modulo invoca
``src/engines/_emngly_runner.py`` (vendorizado, ver su docstring para el
detalle completo de la arquitectura, la reimplementacion local-only de ESM-1b
y el fix de alineamiento de ``structure_emb``) via subprocess, sobre el venv
dedicado ``Settings.EMNGLY_PYTHON_BIN``.

## Rol en el consenso (decision 2026-08-06, ver ``src/engines/ptm_annotation.py``)

A diferencia de MeToken/GlyGen (puramente informativos, nunca deciden
``pasa_umbral``/``consenso``), EMNGly SI es un motor de consenso real --
junto con StackGlyEmbed (promovido del mismo rol informativo a este mismo
papel, ver ``stackglyembed_engine.py``), es la segunda/tercera opinion real
que reemplaza a DeepPTMPred para este tipo especifico (confirmado modelo
muerto, AUROC~=0.51, ``CONSENSUS_EXCLUDED_TYPES`` en ``ptm_annotation.py``).

Aun asi, SIGUE degradando de forma no-fatal (mismo patron que
``stackglyembed_engine.py``/``metoken_engine.py``, NUNCA como
``DeepPTMPredEngine``/``DeepMVPEngine`` que lanzan si la instalacion falta):
si EMNGly no esta instalado en esta maquina, ``ptm_annotation.py`` cae a un
consenso de 2 motores (DeepMVP+StackGlyEmbed) en vez de 3 -- no bloquea el
resto del pipeline. Requiere ``pdb_path`` (Camino PDB unicamente, EMNGly
exige estructura real via MIF) y la tabla de mapeo de posiciones de Fase 1.5
(``position_mapping``, para alinear ``structure_emb`` correctamente -- ver
docstring del runner).
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd

from src.config.settings import Settings
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

OUTPUT_COLUMNS = ["position", "probability"]


def get_emngly_predictions(
    accession: str,
    sequence: str,
    positions: Sequence[int],
    pdb_path: Path,
    position_mapping_path: Path,
    result_dir: Path = None,
    filename_prefix: str = "",
) -> Dict[int, Dict]:
    """Corre EMNGly sobre ``pdb_path``/``sequence`` para ``positions`` y devuelve el resultado por posicion.

    Args:
        accession: Accession de la proteina (para cache/nombre de archivo).
        sequence: Secuencia ATMSEQ COMPLETA de Fase 1.5 (no un fragmento).
        positions: Posiciones 1-based de la Asparagina de cada secuon
            N-X-[S/T] candidato (tipicamente filas con ``tipo_ptm`` en
            {``n_linked_glycosylation``, ``glycosylation_n``}).
        pdb_path: PDB de UNA SOLA cadena (``record.chain_pdb_path`` de Fase
            1.5 -- nunca ``record.pdb_path`` crudo, ver docstring del runner
            para por que).
        position_mapping_path: Ruta al CSV de Fase 1.5
            (``{accession}_position_mapping.csv``, columnas
            ``fasta_position``/``pdb_seqid``) -- necesario para alinear
            ``structure_emb`` con la numeracion real del PDB (ver docstring
            del runner).
        result_dir: Carpeta donde guardar el CSV crudo. Si es ``None``, usa
            ``Settings.FASTA_OUTPUT_DIR``.
        filename_prefix: Prefijo (tipicamente ``f"{accession}_"``) para el
            archivo crudo persistido en ``result_dir``.

    Returns:
        Dict ``{posicion: {"emngly_probability": float}}``, una entrada por
        posicion que EMNGly pudo evaluar (puede faltar alguna si cae fuera
        de rango o sin mapeo PDB valido, ver runner). Dict vacio si EMNGly no
        esta instalado, el subproceso falla, excede el timeout, o la lista
        de posiciones esta vacia -- NUNCA lanza.
    """
    if not positions:
        return {}

    if not _validate_installation():
        return {}

    base_output_dir = Path(result_dir) if result_dir else Settings.FASTA_OUTPUT_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = base_output_dir / f"{filename_prefix}emngly.csv"

    cmd = [
        Settings.EMNGLY_PYTHON_BIN, str(Settings.EMNGLY_RUNNER_SCRIPT.resolve()),
        "--emngly-home", str(Settings.EMNGLY_HOME.resolve()),
        "--accession", accession,
        "--sequence", sequence,
        "--pdb-path", str(Path(pdb_path).resolve()),
        "--position-mapping-csv", str(Path(position_mapping_path).resolve()),
        "--positions", *[str(p) for p in positions],
        "--mif-weights", str(Settings.EMNGLY_MIF_WEIGHTS.resolve()),
        "--esm-checkpoint", str(Settings.EMNGLY_ESM_CHECKPOINT.resolve()),
        "--svm-checkpoint", str(Settings.EMNGLY_SVM_CHECKPOINT.resolve()),
        "--cache-dir", str(Settings.EMNGLY_CACHE_DIR.resolve()),
        "--out-csv", str(out_csv.resolve()),
    ]

    logger.info("Ejecutando EMNGly (consenso de N-glicosilacion) sobre %d posicion(es).", len(positions))
    try:
        # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
        # backend inline que no existe en el conda env aislado) -- ver
        # deepptmpred_engine.py para el caso real que disparo esto.
        env = {**os.environ, "MPLBACKEND": "Agg"}
        subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=Settings.EMNGLY_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "EMNGly termino con exit code %d -- se omite del consenso de N-glicosilacion (no "
            "fatal, degrada a los motores restantes). stderr: %s",
            exc.returncode, (exc.stderr or "<vacio>")[:2000],
        )
        return {}
    except subprocess.TimeoutExpired:
        logger.warning(
            "EMNGly excedio el tiempo limite de %ds -- se omite del consenso de N-glicosilacion "
            "(no fatal).", Settings.EMNGLY_TIMEOUT_SECONDS,
        )
        return {}
    except FileNotFoundError as exc:
        logger.warning(
            "No se pudo invocar EMNGLY_PYTHON_BIN='%s' -- se omite del consenso de N-glicosilacion "
            "(no fatal): %s", Settings.EMNGLY_PYTHON_BIN, exc,
        )
        return {}

    if not out_csv.is_file():
        logger.warning(
            "EMNGly reporto exit code 0 pero no genero '%s' -- se omite del consenso de "
            "N-glicosilacion (no fatal).", out_csv,
        )
        return {}

    try:
        df = pd.read_csv(out_csv)
    except Exception as exc:  # noqa: BLE001 -- degradacion no fatal, cualquier fallo de lectura cuenta
        logger.warning(
            "'%s' no se pudo leer como CSV valido -- se omite del consenso de N-glicosilacion "
            "(no fatal): %s", out_csv, exc,
        )
        return {}

    missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
    if missing:
        logger.warning(
            "'%s' no tiene las columnas esperadas %s (faltan: %s) -- se omite del consenso de "
            "N-glicosilacion (no fatal).", out_csv, OUTPUT_COLUMNS, missing,
        )
        return {}

    return {
        int(r["position"]): {"emngly_probability": float(r["probability"])}
        for _, r in df.iterrows()
    }


def _validate_installation() -> bool:
    """Comprueba que el venv/runner/repo/pesos esten presentes.

    NUNCA fatal (a diferencia de ``DeepPTMPredEngine``/``DeepMVPEngine``):
    EMNGly es un motor de consenso REAL pero opcional -- su ausencia degrada
    el consenso de N-glicosilacion a los motores restantes (DeepMVP +
    StackGlyEmbed), nunca tumba el pipeline principal.
    """
    python_bin = Path(Settings.EMNGLY_PYTHON_BIN)
    if not python_bin.is_file():
        logger.warning(
            "No se encontro el interprete Python de EMNGly en '%s' -- se omite del consenso de "
            "N-glicosilacion (no fatal). Ver README.md - Seccion de instalacion (venv dedicado "
            "'emngly').", python_bin,
        )
        return False
    if not Settings.EMNGLY_RUNNER_SCRIPT.is_file():
        logger.warning(
            "No se encontro el runner de EMNGly en '%s' (deberia venir con el repo, esto no "
            "deberia pasar) -- se omite del consenso de N-glicosilacion (no fatal).",
            Settings.EMNGLY_RUNNER_SCRIPT,
        )
        return False
    if not (Settings.EMNGLY_HOME / "model" / "MIF" / "__init__.py").is_file():
        logger.warning(
            "No se encontro el clon de EMNgly en '%s' -- se omite del consenso de N-glicosilacion "
            "(no fatal). Clona 'git clone https://github.com/StellaHxy/EMNgly' en la raiz del "
            "proyecto (ver README.md).", Settings.EMNGLY_HOME,
        )
        return False
    if not Settings.EMNGLY_MIF_WEIGHTS.is_file():
        logger.warning(
            "No se encontraron los pesos de MIF en '%s' (deberian venir bundled con el clon del "
            "repo, esto no deberia pasar) -- se omite del consenso de N-glicosilacion (no fatal).",
            Settings.EMNGLY_MIF_WEIGHTS,
        )
        return False
    if not Settings.EMNGLY_ESM_CHECKPOINT.is_file():
        logger.warning(
            "No se encontro el checkpoint ESM-1b en '%s' -- se omite del consenso de "
            "N-glicosilacion (no fatal). Descarga manual, ver README.md - Seccion de instalacion.",
            Settings.EMNGLY_ESM_CHECKPOINT,
        )
        return False
    if not Settings.EMNGLY_SVM_CHECKPOINT.is_file():
        logger.warning(
            "No se encontro el checkpoint del SVM en '%s' -- se omite del consenso de "
            "N-glicosilacion (no fatal). Descarga manual desde Google Drive, ver README.md - "
            "Seccion de instalacion.", Settings.EMNGLY_SVM_CHECKPOINT,
        )
        return False
    return True
