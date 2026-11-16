"""Optional corroboration of KINASE SPECIFICITY (Kinase Library), both paths (FASTA and PDB).

Biological coherence analysis 2026-08-07, point 5 (see
``src/engines/_kinase_library_runner.py`` for the complete detail of why
this source CAN give per-site guidance, unlike the polyubiquitin chain
type, a line of reasoning that lived in Phase A/3c, removed from
scope 2026-08-10). It is NEVER a consensus engine:
Kinase Library does not decide ``pasa_umbral``/``consenso``, it only corroborates
phosphorylation sites that the consensus has ALREADY accepted, with the most
probable kinase/family according to the published specificity matrices (Johnson et al.
2023 Nature, Yaron-Barir et al. 2024 Nature).

It only needs the COMPLETE SEQUENCE of the protein (``kl.Substrate`` performs its
own trimming and padding of the 15-mer window around the site), so it
applies equally to the FASTA path (``annotate_fasta_path``) and the PDB path
(``annotate_pdb_path``). It invokes ``src/engines/_kinase_library_runner.py``
via subprocess, on the DEDICATED conda environment
``Settings.KINASE_LIBRARY_PYTHON_BIN`` (the ``numpy``/``pandas`` that the
``kinase-library`` package pins are incompatible with the versions pinned in the
main venv, see requirements.txt).

The same non-deciding pattern as MeToken: any failure
(environment not installed, a subprocess that fails, a timeout, malformed output)
degrades to a warning in the log and returns ``{}``. It NEVER propagates an exception
that would bring down the main pipeline. The caller
(``src/engines/ptm_annotation.py``) decides what to do with an empty result
(it annotates no corroboration on any row, and the rest of the pipeline continues
exactly the same).
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd

from src.config.settings import Settings
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

OUTPUT_COLUMNS = [
    "position", "kinase_library_top_kinase", "kinase_library_top_family",
    "kinase_library_percentile", "kinase_library_top3_kinases",
]


def get_kinase_corroboration(
    sequence: str,
    positions: Sequence[int],
    result_dir: Path = None,
    filename_prefix: str = "",
) -> Dict[int, Dict]:
    """Corre Kinase Library sobre ``sequence`` para ``positions`` y devuelve la corroboracion por posicion.

    Args:
        sequence: Secuencia COMPLETA de la proteina (Fase 1 saneada o Fase
            1.5 ATMSEQ) -- NUNCA un fragmento/peptido: la posicion 1-based
            coincide exactamente sin conversion de offsets.
        positions: Posiciones 1-based de fosforilacion (S/T/Y) ya aceptadas
            por el consenso (``pasa_umbral=True``) a corroborar -- tipicamente
            filas con ``tipo_ptm`` en {``phosphorylation``,
            ``phosphorylation_st``, ``phosphorylation_y``}.
        result_dir: Carpeta donde guardar el CSV crudo. Si es ``None``, usa
            ``Settings.FASTA_OUTPUT_DIR``.
        filename_prefix: Prefijo (tipicamente ``f"{accession}_"``) para el
            archivo crudo persistido en ``result_dir``.

    Returns:
        Dict ``{posicion: {"kinase_library_top_kinase": str, "kinase_library_top_family": str,
        "kinase_library_percentile": float, "kinase_library_top3_kinases": str}}``, una entrada
        por posicion que Kinase Library pudo puntuar (puede faltar alguna si
        el residuo real en esa posicion no es S/T/Y, ver runner). Dict vacio
        si Kinase Library no esta instalado, el subproceso falla, excede el
        timeout, o la lista de posiciones esta vacia -- NUNCA lanza.
    """
    if not positions:
        return {}

    if not Settings.KINASE_LIBRARY_ENABLED:
        return {}

    if not _validate_installation():
        return {}

    base_output_dir = Path(result_dir) if result_dir else Settings.FASTA_OUTPUT_DIR
    base_output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = base_output_dir / f"{filename_prefix}kinase_library.csv"

    cmd = [
        Settings.KINASE_LIBRARY_PYTHON_BIN, str(Settings.KINASE_LIBRARY_RUNNER_SCRIPT.resolve()),
        "--sequence", sequence,
        "--positions", *[str(p) for p in positions],
        "--out-csv", str(out_csv.resolve()),
    ]

    logger.info(
        "Ejecutando Kinase Library (corroboracion de especificidad de quinasa) sobre %d posicion(es).",
        len(positions),
    )
    try:
        # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
        # backend inline que no existe en el conda env aislado) -- ver
        # deepptmpred_engine.py para el caso real que disparo esto.
        env = {**os.environ, "MPLBACKEND": "Agg"}
        subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=Settings.KINASE_LIBRARY_TIMEOUT_SECONDS,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "Kinase Library termino con exit code %d -- se omite la corroboracion de especificidad "
            "de quinasa (no fatal, no afecta consenso/pasa_umbral). stderr: %s",
            exc.returncode, (exc.stderr or "<vacio>")[-6000:],
        )
        return {}
    except subprocess.TimeoutExpired:
        logger.warning(
            "Kinase Library excedio el tiempo limite de %ds -- se omite la corroboracion de "
            "especificidad de quinasa (no fatal).", Settings.KINASE_LIBRARY_TIMEOUT_SECONDS,
        )
        return {}
    except FileNotFoundError as exc:
        logger.warning(
            "No se pudo invocar KINASE_LIBRARY_PYTHON_BIN='%s' -- se omite la corroboracion de "
            "especificidad de quinasa (no fatal): %s", Settings.KINASE_LIBRARY_PYTHON_BIN, exc,
        )
        return {}

    if not out_csv.is_file():
        logger.warning(
            "Kinase Library reporto exit code 0 pero no genero '%s' -- se omite la corroboracion "
            "de especificidad de quinasa (no fatal).", out_csv,
        )
        return {}

    try:
        df = pd.read_csv(out_csv)
    except Exception as exc:  # noqa: BLE001 -- degradacion no fatal, cualquier fallo de lectura cuenta
        logger.warning(
            "'%s' no se pudo leer como CSV valido -- se omite la corroboracion de especificidad de "
            "quinasa (no fatal): %s", out_csv, exc,
        )
        return {}

    missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
    if missing:
        logger.warning(
            "'%s' no tiene las columnas esperadas %s (faltan: %s) -- se omite la corroboracion de "
            "especificidad de quinasa (no fatal).", out_csv, OUTPUT_COLUMNS, missing,
        )
        return {}

    return {
        int(r["position"]): {
            "kinase_library_top_kinase": r["kinase_library_top_kinase"],
            "kinase_library_top_family": r["kinase_library_top_family"],
            "kinase_library_percentile": float(r["kinase_library_percentile"]),
            "kinase_library_top3_kinases": r["kinase_library_top3_kinases"],
        }
        for _, r in df.iterrows()
    }


def _validate_installation() -> bool:
    """Comprueba que el entorno conda dedicado y el runner esten presentes.

    A diferencia de ``DeepMVPEngine``/``DeepPTMPredEngine`` (motores
    OBLIGATORIOS del consenso, su ausencia lanza una excepcion fatal
    accionable), aqui una instalacion faltante NUNCA es fatal -- Kinase
    Library es puramente informativa. Se registra un aviso una unica vez por
    llamada con instrucciones de instalacion.
    """
    python_bin = Path(Settings.KINASE_LIBRARY_PYTHON_BIN)
    if not python_bin.is_file():
        logger.warning(
            "No se encontro el interprete Python de Kinase Library en '%s' -- se omite la "
            "corroboracion de especificidad de quinasa (no fatal). Requiere el entorno conda "
            "'kinase_library' ya instalado (ver README.md: 'conda create -n kinase_library "
            "python=3.10 -y && conda run -n kinase_library pip install kinase-library') o apunta "
            "KINASE_LIBRARY_PYTHON_BIN a otra instalacion.",
            python_bin,
        )
        return False
    if not Settings.KINASE_LIBRARY_RUNNER_SCRIPT.is_file():
        logger.warning(
            "No se encontro el runner de Kinase Library en '%s' (deberia venir con el repo, esto "
            "no deberia pasar) -- se omite la corroboracion de especificidad de quinasa (no fatal).",
            Settings.KINASE_LIBRARY_RUNNER_SCRIPT,
        )
        return False
    return True
