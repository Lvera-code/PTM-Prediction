"""Corroboracion opcional de ESPECIFICIDAD DE QUINASA (Kinase Library), ambos caminos (FASTA y PDB).

Analisis de coherencia biologica 2026-08-07, punto 5 (ver
``src/engines/_kinase_library_runner.py`` para el detalle completo de por
que esta fuente SI puede orientar por-sitio, a diferencia del punto 4/
``ubiquitin_sumo.CHAIN_TYPE_DISCLAIMER``). NUNCA es un motor de consenso --
Kinase Library no decide ``pasa_umbral``/``consenso``, solo corrobora
sitios de fosforilacion que el consenso YA acepto, con la quinasa/familia
mas probable segun las matrices de especificidad publicadas (Johnson et al.
2023 Nature, Yaron-Barir et al. 2024 Nature).

Igual que StackGlyEmbed, solo necesita la SECUENCIA COMPLETA de la proteina
(``kl.Substrate`` hace su propio recorte+relleno de la ventana de 15-mer
alrededor del sitio) -- aplica igual a Camino FASTA (``annotate_fasta_path``)
y Camino PDB (``annotate_pdb_path``). Invoca
``src/engines/_kinase_library_runner.py`` via subprocess, sobre el entorno
conda DEDICADO ``Settings.KINASE_LIBRARY_PYTHON_BIN`` (``numpy``/``pandas``
que fija el paquete ``kinase-library`` son incompatibles con las versiones
fijadas del venv principal, ver requirements.txt).

Mismo patron no-decisorio que GlyGen/MeToken/StackGlyEmbed: cualquier fallo
(entorno no instalado, subproceso que revienta, timeout, salida malformada)
degrada a un aviso en el log y devuelve ``{}`` -- NUNCA propaga una excepcion
que tumbe el pipeline principal. Quien llama
(``src/engines/ptm_annotation.py``) decide que hacer con un resultado vacio
(no anota corroboracion en ninguna fila, el resto del pipeline sigue
exactamente igual).
"""

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
            1.5 ATMSEQ) -- NUNCA un fragmento/peptido, mismo criterio que
            ``stackglyembed_engine.get_nglyco_corroboration``: la posicion
            1-based coincide exactamente sin conversion de offsets.
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
        subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=Settings.KINASE_LIBRARY_TIMEOUT_SECONDS,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "Kinase Library termino con exit code %d -- se omite la corroboracion de especificidad "
            "de quinasa (no fatal, no afecta consenso/pasa_umbral). stderr: %s",
            exc.returncode, (exc.stderr or "<vacio>")[:2000],
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
