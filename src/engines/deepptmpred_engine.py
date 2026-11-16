"""Phase 2: PTM zone prediction via LOCAL DeepPTMPred (subprocess over our own runner).

Engine 2 of 2 of the consensus, PDB path only (DeepPTMPred demands a
mandatory ``pdb_path``, with no sequence-only mode, decision 2026-07-26).

Verified by reading github.com/kuikui-wang/DeepPTMPred directly
(README.md, pred/train_PTM/predict.py, pred/train_PTM/e2_single_data.py,
pred/train_PTM/environment.yml) on 2026-07-27, not from a search engine summary.
A real repo, with weights included (.h5 per type, ~19MB each), but WITHOUT a
declared licence (unlike DeepMVP, GPL-3.0). The licence was subsequently
confirmed as CC BY-NC 4.0, the same terms as the paper, directly by a
corresponding author, so non-commercial research use is covered.

Unlike DeepMVPEngine, this wrapper does NOT invoke ``predict.py``
directly by subprocess: that script has no real CLI (``ptm_type``/
``pdb_path`` hardcoded in its ``__main__`` block, confirmed by reading the
code). Instead it invokes ``src/engines/_deepptmpred_runner.py`` (our own,
see its docstring for the complete detail of why it exists and
which bug of ``e2_single_data.py`` it avoids).

The key structural difference with DeepMVP: DeepPTMPred predicts ONE PTM type
per invocation (not all 17 at once), so this engine invokes the runner once
PER TYPE, per accession, and concatenates the results into a single
DataFrame per accession (the same final shape as ``DeepMVPEngine.run()``, so
that the Phase 3 core can treat them symmetrically).
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

# Confirmed in src/engines/_deepptmpred_runner.py::OUTPUT_COLUMNS (the
# 'prediction' column of the original repo, with a hardcoded uncalibrated
# cutoff of 0.5, is deliberately discarded; see the runner's docstring).
OUTPUT_COLUMNS = ["protein_id", "position", "residue", "probability", "ptm_type"]


class DeepPTMPredEngine(BaseEngine[StructureRecord, pd.DataFrame]):
    """Runs DeepPTMPred LOCALLY (subprocess over our own runner) per accession.

    Each :class:`StructureRecord` (the output of Phase 1.5) is predicted once
    PER PTM TYPE (``Settings.DEEPPTMPRED_PTM_TYPES``, 17 invocations),
    concatenating the result into a single DataFrame per accession.
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
        """Checks that the cloned repo and the ESM-2 checkpoint are present.

        Raises:
            DeepPTMPredExecutionError: With an actionable message if the repo's
                ``predict.py`` or the ESM-2 checkpoint is missing (both
                manual installation steps, see README.md, the
                Installation section).
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
        """Runs DeepPTMPred locally (all 17 PTM types) over each structure of ``items``.

        Args:
            items: :class:`StructureRecord` from Phase 1.5 (PDB path
                only; this engine has no sequence-only mode).
            output_dir: Base folder in which to save the raw CSVs
                (one subfolder per accession). If ``None``, it uses
                ``Settings.FASTA_OUTPUT_DIR``.

        Returns:
            A list of DataFrames with ``OUTPUT_COLUMNS`` (concatenated from the
            17 PTM types, one row per candidate site, without any
            threshold filtering, which is the responsibility of the Phase 3 core), in
            the same order as ``items``.

        Raises:
            DeepPTMPredExecutionError: If the local installation (repo or
                ESM-2 checkpoint) does not exist, a subprocess fails or
                exceeds the timeout, or an output does not have the
                expected format.
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
        # MPLBACKEND is inherited from the parent process (Jupyter/Colab sets it to
        # 'module://matplotlib_inline.backend_inline' for inline plots) and
        # that backend does not exist in DeepPTMPred's isolated conda env:
        # predict.py imports matplotlib.pyplot and fails with a ValueError when
        # trying to set that invalid backend. 'Agg' is headless and always
        # valid, regardless of the invoking process.
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
