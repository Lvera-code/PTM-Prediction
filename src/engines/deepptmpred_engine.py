"""Fase 2: prediccion de zonas PTM via DeepPTMPred LOCAL (subprocess sobre runner propio).

Motor 2 de 2 del consenso, Camino PDB unicamente (DeepPTMPred exige
``pdb_path`` obligatorio -- sin modo solo-secuencia, decision 2026-07-26).

Verificado leyendo directamente github.com/kuikui-wang/DeepPTMPred
(README.md, pred/train_PTM/predict.py, pred/train_PTM/e2_single_data.py,
pred/train_PTM/environment.yml) el 2026-07-27 -- no por resumen de buscador.
Repo real, pesos SI incluidos (.h5 por tipo, ~19MB c/u), pero SIN licencia
declarada (a diferencia de DeepMVP, GPL-3.0) -- verificar con Carlos antes
de cualquier uso mas alla de investigacion/TFG.

A diferencia de DeepMVPEngine, este wrapper NO invoca ``predict.py``
directamente por subprocess: ese script no tiene CLI real (``ptm_type``/
``pdb_path`` hardcodeados en su bloque ``__main__``, confirmado leyendo el
codigo). En su lugar invoca ``src/engines/_deepptmpred_runner.py`` (propio de
este proyecto, ver su docstring para el detalle completo de por que existe y
que bug evita de ``e2_single_data.py``).

Diferencia estructural clave con DeepMVP: DeepPTMPred predice UN tipo de PTM
por invocacion (no los 17 de una vez) -- este engine invoca el runner una
vez POR TIPO, por accession, y concatena los resultados en un unico
DataFrame por accession (mismo shape final que ``DeepMVPEngine.run()``, para
que el nucleo de Fase 3 pueda tratarlos de forma simetrica).

NO PROBADO TODAVIA contra el entorno real (PyRosetta/TensorFlow 2.15/
checkpoint ESM-2 no instalados en esta maquina, ver STATUS.md) -- los tests
mockean ``subprocess.run``.
"""

import os
import subprocess
from pathlib import Path
from typing import List, Sequence

import pandas as pd

from src.config.settings import Settings
from src.engines.base_engine import BaseEngine
from src.utils.exceptions import DeepPTMPredExecutionError
from src.utils.logger_config import setup_logger
from src.utils.structure_parser import StructureRecord

logger = setup_logger(__name__)

# Confirmadas en src/engines/_deepptmpred_runner.py::OUTPUT_COLUMNS (la
# columna 'prediction' del repo original, cutoff 0.5 hardcodeado sin
# calibrar, se descarta deliberadamente -- ver docstring del runner).
OUTPUT_COLUMNS = ["protein_id", "position", "residue", "probability", "ptm_type"]


class DeepPTMPredEngine(BaseEngine[StructureRecord, pd.DataFrame]):
    """Ejecuta DeepPTMPred LOCALMENTE (subprocess sobre el runner propio) por accession.

    Cada :class:`StructureRecord` (salida de Fase 1.5) se predice una vez
    POR TIPO de PTM (``Settings.DEEPPTMPRED_PTM_TYPES``, 17 invocaciones),
    concatenando el resultado en un unico DataFrame por accession.
    """

    def __init__(
        self,
        train_ptm_dir: Path = Settings.DEEPPTMPRED_TRAIN_PTM_DIR,
        runner_script: Path = Settings.DEEPPTMPRED_RUNNER_SCRIPT,
        python_bin: str = Settings.DEEPPTMPRED_PYTHON_BIN,
        esm_checkpoint: Path = Settings.DEEPPTMPRED_ESM_CHECKPOINT,
        custom_esm_dir: Path = Settings.DEEPPTMPRED_CUSTOM_ESM_DIR,
        ptm_types: Sequence[str] = Settings.DEEPPTMPRED_PTM_TYPES,
    ):
        self._train_ptm_dir = Path(train_ptm_dir)
        self._runner_script = Path(runner_script)
        self._python_bin = python_bin
        self._esm_checkpoint = Path(esm_checkpoint)
        self._custom_esm_dir = Path(custom_esm_dir)
        self._ptm_types = list(ptm_types)

    def _validate_installation(self) -> None:
        """Comprueba que el repo clonado y el checkpoint ESM-2 esten presentes.

        Raises:
            DeepPTMPredExecutionError: Con un mensaje accionable si falta
                ``predict.py`` del repo o el checkpoint ESM-2 (ambos pasos
                de instalacion manual, ver README.md - Seccion de
                Instalacion).
        """
        if not (self._train_ptm_dir / "predict.py").is_file():
            raise DeepPTMPredExecutionError(
                f"No se encontro la instalacion local de DeepPTMPred en "
                f"'{self._train_ptm_dir}'. Clona el repo con 'git clone "
                "https://github.com/kuikui-wang/DeepPTMPred' en la raiz del proyecto "
                "(o apunta DEEPPTMPRED_TRAIN_PTM_DIR a 'DeepPTMPred/pred/train_PTM') y "
                "vuelve a intentarlo. Ver README.md - Seccion de Instalacion."
            )
        if not self._esm_checkpoint.is_file():
            raise DeepPTMPredExecutionError(
                f"No se encontro el checkpoint ESM-2 en '{self._esm_checkpoint}'. "
                "Descargalo (ver README.md del repo, seccion 'esm model') y apunta "
                "DEEPPTMPRED_ESM_CHECKPOINT a su ubicacion."
            )

    def run(self, items: Sequence[StructureRecord], output_dir: Path = None) -> List[pd.DataFrame]:
        """Corre DeepPTMPred localmente (los 17 tipos de PTM) sobre cada estructura de ``items``.

        Args:
            items: :class:`StructureRecord` de Fase 1.5 (Camino PDB
                unicamente -- este motor no tiene modo solo-secuencia).
            output_dir: Carpeta base donde guardar los CSV crudos
                (subcarpeta por accession). Si es ``None``, usa
                ``Settings.FASTA_OUTPUT_DIR``.

        Returns:
            Lista de DataFrames con ``OUTPUT_COLUMNS`` (concatenado de los
            17 tipos de PTM, una fila por sitio candidato, sin ningun
            filtrado de umbral -- responsabilidad del nucleo de Fase 3), en
            el mismo orden que ``items``.

        Raises:
            DeepPTMPredExecutionError: Si la instalacion local (repo o
                checkpoint ESM-2) no existe, algun subproceso falla o
                excede el timeout, o alguna salida no tiene el formato
                esperado.
        """
        self._validate_installation()
        base_output_dir = Path(output_dir) if output_dir else Settings.FASTA_OUTPUT_DIR

        return [
            self._run_single(record, base_output_dir / record.accession)
            for record in items
        ]

    def _run_single(self, record: StructureRecord, result_dir: Path) -> pd.DataFrame:
        result_dir.mkdir(parents=True, exist_ok=True)

        per_type_frames = [
            self._run_single_ptm_type(record, ptm_type, result_dir)
            for ptm_type in self._ptm_types
        ]
        return pd.concat(per_type_frames, ignore_index=True) if per_type_frames else pd.DataFrame(
            columns=OUTPUT_COLUMNS
        )

    def _run_single_ptm_type(self, record: StructureRecord, ptm_type: str, result_dir: Path) -> pd.DataFrame:
        out_csv = result_dir / f"{ptm_type}.csv"

        cmd = [
            self._python_bin, str(self._runner_script.resolve()),
            "--train-ptm-dir", str(self._train_ptm_dir.resolve()),
            "--protein-id", record.accession,
            "--sequence", record.sequence,
            "--pdb-path", str(record.pdb_path.resolve()),
            "--ptm-type", ptm_type,
            "--esm-checkpoint", str(self._esm_checkpoint.resolve()),
            "--custom-esm-dir", str(self._custom_esm_dir.resolve()),
            "--out-csv", str(out_csv.resolve()),
        ]

        logger.info(
            "Ejecutando DeepPTMPred local (%s) para '%s': %s",
            ptm_type, record.accession, " ".join(cmd),
        )
        # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a
        # 'module://matplotlib_inline.backend_inline' para graficos inline) y
        # ese backend no existe en el conda env aislado de DeepPTMPred --
        # predict.py importa matplotlib.pyplot y revienta con ValueError al
        # intentar fijar ese backend invalido. 'Agg' es headless y siempre
        # valido, sin importar el proceso que invoque.
        env = {**os.environ, "MPLBACKEND": "Agg"}
        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=Settings.DEEPPTMPRED_TIMEOUT_SECONDS,
                env=env,
            )
        except subprocess.CalledProcessError as exc:
            raise DeepPTMPredExecutionError(
                f"DeepPTMPred ({ptm_type}) termino con exit code {exc.returncode} para "
                f"'{record.accession}'. stderr: {(exc.stderr or '<vacio>')[-6000:]}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise DeepPTMPredExecutionError(
                f"DeepPTMPred ({ptm_type}) excedio el tiempo limite de "
                f"{Settings.DEEPPTMPRED_TIMEOUT_SECONDS}s para '{record.accession}'. Aumenta "
                "DEEPPTMPRED_TIMEOUT_SECONDS si el hardware es lento (CPU vs GPU) o la "
                "proteina es muy larga."
            ) from exc

        return self._load_ptm_predictions(out_csv, record.accession, ptm_type)

    @staticmethod
    def _load_ptm_predictions(out_csv: Path, accession: str, ptm_type: str) -> pd.DataFrame:
        if not out_csv.is_file():
            raise DeepPTMPredExecutionError(
                f"El runner de DeepPTMPred no genero '{out_csv}' para '{accession}' "
                f"(tipo '{ptm_type}'). El subproceso reporto exit code 0 pero el archivo de "
                "salida esperado no existe: revisa el log completo del subproceso."
            )

        df = pd.read_csv(out_csv)
        missing = [c for c in OUTPUT_COLUMNS if c not in df.columns]
        if missing:
            raise DeepPTMPredExecutionError(
                f"'{out_csv}' no tiene las columnas esperadas {OUTPUT_COLUMNS} "
                f"(faltan: {missing}, columnas reales: {list(df.columns)})."
            )
        return df[OUTPUT_COLUMNS]
