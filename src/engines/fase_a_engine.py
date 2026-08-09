"""Fase A: modelado estructural real por sitio PTM (subprocess sobre `_fase_a_runner.py`).

Conecta Extension 3 (ddG, clase 1) + Fase A clase 2 (glicosilacion) + clase 3
(ubiquitinacion/sumoilacion) -- todas implementadas y verificadas con
corridas reales desde 2026-07-28/08-01 (ver STATUS.md), pero nunca conectadas
al orquestador principal hasta ahora (decision 2026-07-27: D era un filtro de
responsabilidad unica sin rutas a fases entonces inexistentes). Revertido
2026-08-03 para la demo del 2026-08-10: ver
``01-Proyectos/PTM-Prediction/Decisiones/`` para la decision completa.

NUNCA es un motor de consenso ni cambia ``pasa_umbral``: opera DESPUES del
filtro de Fase 3 (D), sobre un subconjunto ya aceptado de sitios (ver
``ptm_annotation.select_fase_a_candidates`` -- top-N por tipo, no todos los
sitios aceptados, por costo computacional real: modelar los ~572 sitios
reales que acepta un caso como Tau tardaria horas).

Requiere PyRosetta (conda env ``deepptmpred``, mismo env ya usado para
DeepPTMPred) -- NUNCA se importa ``src.structural.*`` desde este modulo ni
desde el proceso principal del pipeline. Cada sitio se modela en un
subprocess FRESCO (``src/engines/_fase_a_runner.py``), un proceso por sitio:
las 3 clases de Fase A inicializan PyRosetta con flags incompatibles entre
si (ver docstring de ``fase_a_dispatch.py``), por eso no se pueden agrupar
varios sitios de clases distintas en el mismo proceso.

Degradacion NO fatal en todos los modos de fallo (interprete/runner
ausentes, subprocess que revienta, timeout, JSON de salida ausente o
malformado) -- mismo patron que ``stackglyembed_engine.py``/
``metoken_engine.py``: nunca lanza una excepcion que tumbe el pipeline
principal, un sitio individual que falla no debe impedir que el resto del
reporte final se genere. Quien llama decide que hacer con
``estado="error"`` (se documenta en el reporte, no se oculta).
"""

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

from src.config.settings import Settings
from src.engines.base_engine import BaseEngine
from src.utils.logger_config import setup_logger

logger = setup_logger(__name__)

_ESTADO_NO_DISPONIBLE = "no_disponible"


@dataclass(frozen=True)
class FaseASiteRequest:
    """Un sitio PTM real a modelar: accession + pdb_path de la cadena + posicion + tipo."""

    accession: str
    pdb_path: Path
    position: int
    ptm_type: str
    uniprot_accession: Optional[str] = None


def _empty_result(estado: str, error: Optional[str] = None) -> dict:
    return {**Settings.FASE_A_RESULT_TEMPLATE, "estado": estado, "error": error}


class FaseAEngine(BaseEngine[FaseASiteRequest, dict]):
    """Modela cada :class:`FaseASiteRequest` individualmente via subprocess."""

    def __init__(
        self,
        python_bin: str = Settings.FASE_A_PYTHON_BIN,
        runner_script: Path = Settings.FASE_A_RUNNER_SCRIPT,
        timeout_seconds: int = Settings.FASE_A_TIMEOUT_SECONDS,
    ):
        self._python_bin = python_bin
        self._runner_script = Path(runner_script)
        self._timeout_seconds = timeout_seconds

    def _validate_installation(self) -> bool:
        python_bin = Path(self._python_bin)
        if not python_bin.is_file():
            logger.warning(
                "No se encontro el interprete Python de Fase A (con PyRosetta) en '%s' -- se "
                "omite el modelado estructural para todo el lote (no fatal). Requiere el conda "
                "env 'deepptmpred' ya instalado (ver README.md) o apunta FASE_A_PYTHON_BIN a otra "
                "instalacion.", python_bin,
            )
            return False
        if not self._runner_script.is_file():
            logger.warning(
                "No se encontro el runner de Fase A en '%s' (deberia venir con el repo) -- se "
                "omite el modelado estructural para todo el lote (no fatal).", self._runner_script,
            )
            return False
        return True

    def run(self, items: Sequence[FaseASiteRequest], output_dir: Optional[Path] = None) -> List[dict]:
        """Modela cada sitio de ``items`` en su propio subprocess.

        Args:
            items: Sitios a modelar (tipicamente la salida de
                ``ptm_annotation.select_fase_a_candidates`` -- ya acotada a
                los tipos con Fase A implementada y a un top-N por tipo).
            output_dir: Carpeta base donde persistir los PDB/JSON generados
                (subcarpeta por accession). Si es ``None``, usa
                ``Settings.FASTA_OUTPUT_DIR``.

        Returns:
            Lista de dicts (ver ``fase_a_dispatch.run_fase_a_for_site`` para
            las claves), mismo orden que ``items``. Si ``Settings.FASE_A_ENABLED``
            es ``False`` o la instalacion no esta disponible, cada entrada
            tiene ``estado="no_disponible"`` -- nunca lanza.
        """
        if not Settings.FASE_A_ENABLED:
            return [_empty_result(_ESTADO_NO_DISPONIBLE) for _ in items]
        if not items:
            return []
        if not self._validate_installation():
            return [_empty_result(_ESTADO_NO_DISPONIBLE) for _ in items]

        base_output_dir = Path(output_dir) if output_dir else Settings.FASTA_OUTPUT_DIR
        return [self._run_single(item, base_output_dir / item.accession / "fase_a") for item in items]

    def _run_single(self, item: FaseASiteRequest, result_dir: Path) -> dict:
        result_dir.mkdir(parents=True, exist_ok=True)
        out_pdb = result_dir / f"{item.ptm_type}_{item.position}.pdb"
        out_json = result_dir / f"{item.ptm_type}_{item.position}.json"

        cmd = [
            self._python_bin, str(self._runner_script.resolve()),
            "--pdb-path", str(Path(item.pdb_path).resolve()),
            "--position", str(item.position),
            "--ptm-type", item.ptm_type,
            "--out-pdb", str(out_pdb.resolve()),
            "--out-json", str(out_json.resolve()),
        ]
        if item.uniprot_accession:
            cmd += ["--uniprot-accession", item.uniprot_accession]

        logger.info(
            "Ejecutando Fase A (%s @ %d) para '%s': %s",
            item.ptm_type, item.position, item.accession, " ".join(cmd),
        )
        try:
            # MPLBACKEND se hereda del proceso padre (Jupyter/Colab lo fija a un
            # backend inline que no existe en el conda env aislado) -- ver
            # deepptmpred_engine.py para el caso real que disparo esto.
            env = {**os.environ, "MPLBACKEND": "Agg"}
            subprocess.run(
                cmd, check=True, capture_output=True, text=True, timeout=self._timeout_seconds,
                env=env,
            )
        except subprocess.CalledProcessError as exc:
            logger.warning(
                "Fase A (%s @ %d) termino con exit code %d -- se omite el modelado de este sitio "
                "(no fatal). stderr: %s",
                item.ptm_type, item.position, exc.returncode, (exc.stderr or "<vacio>")[:2000],
            )
            return _empty_result("error", f"exit code {exc.returncode}: {(exc.stderr or '')[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(
                "Fase A (%s @ %d) excedio el tiempo limite de %ds -- se omite el modelado de este "
                "sitio (no fatal).", item.ptm_type, item.position, self._timeout_seconds,
            )
            return _empty_result("error", f"timeout tras {self._timeout_seconds}s")

        if not out_json.is_file():
            logger.warning(
                "Fase A (%s @ %d) reporto exit code 0 pero no genero '%s' -- se omite (no fatal).",
                item.ptm_type, item.position, out_json,
            )
            return _empty_result("error", f"'{out_json}' no generado")

        try:
            with open(out_json) as f:
                return json.load(f)
        except Exception as exc:  # noqa: BLE001 -- degradacion no fatal, cualquier fallo de lectura cuenta
            logger.warning(
                "'%s' no se pudo leer como JSON valido -- se omite el modelado de este sitio "
                "(no fatal): %s", out_json, exc,
            )
            return _empty_result("error", f"JSON invalido: {exc}")
