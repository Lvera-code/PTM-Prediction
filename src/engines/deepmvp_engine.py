"""Phase 2: PTM zone prediction via LOCAL DeepMVP (pure subprocess).

The only engine of the FASTA path; engine 1 of 2 (alongside DeepPTMPred) of
the consensus on the PDB path.

Verified by reading github.com/bzhanglab/DeepMVP directly (README.md,
DeepMVP.py, lib/PTModels.py, lib/Metrics.py) on 2026-07-27, not from a search
engine summary. A real, packaged repo (``requirements.txt`` + ``environment.yml``,
GPL-3.0), with weights NOT included (a separate manual download, see
``Settings.DEEPMVP_DOWNLOAD_URL``).

The real CLI (confirmed in ``DeepMVP.py``, matching the README example):
``python DeepMVP.py predict -m <model_dir> -d <fasta> -t 2 -o <out_dir>``.
Without ``-i``: in task mode 2 (PTM site prediction) without ``-i``, DeepMVP extracts
by itself ALL the candidate sites of each supported PTM type
directly from the FASTA (-d), one for each residue of each model's
alphabet (K for acetylation, S/T for phosphorylation_st, etc.), which is exactly
this wrapper's use case (annotating the whole input proteome, not a
single site).

Fixed output (not configurable via DeepMVP flags, confirmed in
``lib/PTModels.py::ptm_prediction_for_multiple_ptms``, with the prefix hardcoded to
``"site_prediction"``): ``<out_dir>/site_prediction.tsv`` with the columns
``protein | aa | pos | x | y_pred | fpr | ptm`` (concatenated from the 8
model subfolders per PTM type: acetylation_k, glycosylation_n,
methylation_k, methylation_r, phosphorylation_st, phosphorylation_y,
sumoylation_k, ubiquitination_k, that is, 8 residue-specific models that
cover the project's 6 biological PTM categories).

Tolerance to non-canonical residues (confirmed by reading
``lib/PeptideEncode.py::encodePeptideOneHot``): DeepMVP does NOT abort on a
residue outside the standard alphabet; 'X' is encoded as a zero vector, and
any other unrecognized character as a vector of 0.5 with a printed
warning (the corresponding ``exit(1)`` is commented out in the real code).
This is MORE PERMISSIVE than the default fatal-rejection policy of
``src.utils.fasta_parser`` (documented there as conservative, pending
verification), a decision still to be discussed, so the Phase 1 policy is not
relaxed unilaterally here.
"""

import os
import subprocess
from pathlib import Path
from typing import List, Sequence

import pandas as pd

from src.config.settings import Settings
from src.engines.base_engine import BaseEngine
from src.utils.exceptions import DeepMVPExecutionError
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

# Columnas confirmadas leyendo lib/PTModels.py (ptm_predict/dl_models_predict)
# y el ejemplo de salida real documentado en README.md.
OUTPUT_COLUMNS = ["protein", "aa", "pos", "x", "y_pred", "fpr", "ptm"]
SITE_PREDICTION_FILENAME = "site_prediction.tsv"


class DeepMVPEngine(BaseEngine[str, pd.DataFrame]):
    """Ejecuta DeepMVP LOCALMENTE (subprocess) para prediccion de sitios PTM.

    Cada llamada a :meth:`run` procesa un FASTA completo (todas las
    accessions que contenga) en una unica invocacion del CLI -- DeepMVP no
    tiene un modo de streaming por accession, procesa el ``-d`` completo de
    una vez.
    """

    def __init__(
        self,
        deepmvp_home: Path = Settings.DEEPMVP_HOME,
        script_name: str = Settings.DEEPMVP_SCRIPT_NAME,
        python_bin: str = Settings.DEEPMVP_PYTHON_BIN,
        model_dir: Path = Settings.DEEPMVP_MODEL_DIR,
    ):
        self._deepmvp_home = Path(deepmvp_home)
        self._script = self._deepmvp_home / script_name
        self._python_bin = python_bin
        self._model_dir = Path(model_dir)

    @property
    def script(self) -> Path:
        """Ruta resuelta a ``DeepMVP.py`` de la instalacion local."""
        return self._script

    def _validate_installation(self) -> None:
        """Comprueba que el repo clonado y los pesos descargados esten presentes.

        Raises:
            DeepMVPExecutionError: Con un mensaje accionable si falta el
                script del repo o la carpeta de pesos (ambos pasos de
                instalacion manual, ver README.md - Seccion de Instalacion).
        """
        if not self._script.is_file():
            raise DeepMVPExecutionError(
                f"No se encontro la instalacion local de DeepMVP en '{self._script}'. "
                "Clona el repo con 'git clone https://github.com/bzhanglab/DeepMVP' en la "
                "raiz del proyecto (o apunta DEEPMVP_HOME a su ubicacion) y vuelve a "
                "intentarlo. Ver README.md - Seccion de Instalacion."
            )
        if not self._model_dir.is_dir() or not any(self._model_dir.iterdir()):
            raise DeepMVPExecutionError(
                f"No se encontraron pesos de DeepMVP en '{self._model_dir}'. Descarga los "
                f"pesos pre-entrenados desde {Settings.DEEPMVP_DOWNLOAD_URL}, descomprime el "
                ".tar.gz y coloca su contenido en esa carpeta (o apunta DEEPMVP_MODEL_DIR a "
                "su ubicacion). Ver README.md - Seccion de Instalacion."
            )

    def run(self, items: Sequence[str], output_dir: Path = None) -> List[pd.DataFrame]:
        """Corre DeepMVP localmente (tarea 2: PTM site prediction) sobre cada FASTA de ``items``.

        Args:
            items: Rutas locales a archivos FASTA saneados (Camino FASTA:
                salida de Fase 1; Camino PDB: FASTA ATMSEQ derivado de
                Fase 1.5).
            output_dir: Carpeta base donde guardar los artefactos crudos
                (subcarpeta por ``stem`` de cada FASTA). Si es ``None``, usa
                ``Settings.FASTA_OUTPUT_DIR``.

        Returns:
            Lista de DataFrames con ``OUTPUT_COLUMNS`` (una fila por sitio
            PTM candidato, sin ningun filtrado de umbral -- eso es
            responsabilidad de Fase 3 nucleo, no de este motor), en el mismo
            orden que ``items``.

        Raises:
            DeepMVPExecutionError: Si la instalacion local (repo o pesos) no
                existe, el subproceso falla, excede el timeout, o la salida
                no tiene el formato esperado.
        """
        self._validate_installation()
        base_output_dir = Path(output_dir) if output_dir else Settings.FASTA_OUTPUT_DIR

        return [
            self._run_single(fasta_path, base_output_dir / Path(fasta_path).stem)
            for fasta_path in items
        ]

    def _run_single(self, fasta_path: str, result_dir: Path) -> pd.DataFrame:
        fasta = Path(fasta_path)
        if not fasta.is_file():
            raise FileNotFoundError(f"No se encontro el FASTA de entrada: {fasta}")

        result_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self._python_bin, str(self._script.resolve()),
            "predict",
            "-m", str(self._model_dir.resolve()),
            "-d", str(fasta.resolve()),
            "-t", Settings.DEEPMVP_TASK_PTM_SITE_PREDICTION,
            "-o", str(result_dir.resolve()),
        ]

        logger.info("Ejecutando DeepMVP local para '%s': %s", fasta.name, " ".join(cmd))
        try:
            # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
            # backend inline que no existe en el conda env aislado) -- ver
            # deepptmpred_engine.py para el caso real que disparo esto.
            env = {**os.environ, "MPLBACKEND": "Agg"}
            subprocess.run(
                cmd,
                cwd=str(self._deepmvp_home),
                check=True,
                capture_output=True,
                text=True,
                timeout=Settings.DEEPMVP_TIMEOUT_SECONDS,
                env=env,
            )
        except subprocess.CalledProcessError as exc:
            raise DeepMVPExecutionError(
                f"DeepMVP termino con exit code {exc.returncode} para '{fasta.name}'. "
                f"stderr: {(exc.stderr or '<vacio>')[:2000]}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise DeepMVPExecutionError(
                f"DeepMVP excedio el tiempo limite de {Settings.DEEPMVP_TIMEOUT_SECONDS}s "
                f"para '{fasta.name}'. Aumenta DEEPMVP_TIMEOUT_SECONDS si el FASTA tiene "
                "muchas accessions o el hardware es lento (CPU vs GPU)."
            ) from exc

        return self._load_site_predictions(result_dir, fasta.name)

    @staticmethod
    def _load_site_predictions(result_dir: Path, fasta_name: str) -> pd.DataFrame:
        """Lee y valida ``site_prediction.tsv`` (nombre fijo, ver docstring del modulo)."""
        tsv_path = result_dir / SITE_PREDICTION_FILENAME
        if not tsv_path.is_file():
            raise DeepMVPExecutionError(
                f"DeepMVP no genero '{SITE_PREDICTION_FILENAME}' para '{fasta_name}' en "
                f"'{result_dir}'. El subproceso reporto exit code 0 pero el archivo de "
                "salida esperado no existe: revisa el log completo del subproceso."
            )

        df = pd.read_csv(tsv_path, sep="\t")
        missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
        if missing:
            raise DeepMVPExecutionError(
                f"'{tsv_path}' no tiene las columnas esperadas {OUTPUT_COLUMNS} "
                f"(faltan: {missing}, columnas reales: {list(df.columns)}). "
                "Puede que una version mas nueva de DeepMVP haya cambiado su formato de salida."
            )
        return df[OUTPUT_COLUMNS]
