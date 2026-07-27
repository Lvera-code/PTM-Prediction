"""Fase 3a: prediccion de zonas PTM via DeepMVP LOCAL (subprocess puro).

Motor unico del Camino FASTA; motor 1 de 2 (junto a DeepPTMPred, no
implementado todavia) del consenso en Camino PDB.

Verificado leyendo directamente github.com/bzhanglab/DeepMVP (README.md,
DeepMVP.py, lib/PTModels.py, lib/Metrics.py) el 2026-07-27 -- no por resumen
de buscador. Repo real y empaquetado (``requirements.txt`` + ``environment.yml``,
GPL-3.0), pesos NO incluidos (descarga manual separada, ver
``Settings.DEEPMVP_DOWNLOAD_URL``).

CLI real (confirmado en ``DeepMVP.py``, coincide con el ejemplo del README):
``python DeepMVP.py predict -m <model_dir> -d <fasta> -t 2 -o <out_dir>``.
Sin ``-i``: en modo tarea 2 (PTM site prediction) sin ``-i``, DeepMVP extrae
por si mismo TODOS los sitios candidatos de cada tipo de PTM soportado
directamente del FASTA (-d), uno por cada residuo del alfabeto de cada
modelo (K para acetilacion, S/T para fosforilacion_st, etc.) -- exactamente
el caso de uso de este wrapper (anotar todo el proteoma de entrada, no un
sitio puntual).

Salida fija (no configurable via flags de DeepMVP, confirmado en
``lib/PTModels.py::ptm_prediction_for_multiple_ptms``, prefix hardcodeado a
``"site_prediction"``): ``<out_dir>/site_prediction.tsv`` con columnas
``protein | aa | pos | x | y_pred | fpr | ptm`` (concatenado de las 8
subcarpetas de modelo por tipo de PTM: acetylation_k, glycosylation_n,
methylation_k, methylation_r, phosphorylation_st, phosphorylation_y,
sumoylation_k, ubiquitination_k -- 8 modelos especificos de residuo que
cubren las 6 categorias biologicas de PTM del proyecto).

Tolerancia a residuos no canonicos (confirmado leyendo
``lib/PeptideEncode.py::encodePeptideOneHot``): DeepMVP NO aborta ante un
residuo fuera del alfabeto estandar -- 'X' se codifica como vector de ceros,
cualquier otro caracter no reconocido como vector de 0.5 con un aviso
impreso (el ``exit(1)`` correspondiente esta comentado en el codigo real).
Esto es MAS PERMISIVO que la politica de rechazo fatal por defecto de
``src.utils.fasta_parser`` (documentada ahi como conservadora, pendiente de
verificar) -- decision pendiente de discutir con el usuario, no se relaja
la politica de Fase 1 unilateralmente aqui.
"""

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
            subprocess.run(
                cmd,
                cwd=str(self._deepmvp_home),
                check=True,
                capture_output=True,
                text=True,
                timeout=Settings.DEEPMVP_TIMEOUT_SECONDS,
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
