"""Fase A: enrutador de un sitio PTM real a su modulo de modelado estructural.

Conectado al pipeline principal 2026-08-03 (revierte la decision 2026-07-27
de que D -- ``ptm_annotation.apply_workflow_filter`` -- no rutea a Extension
3/Fase A porque esas fases no existian todavia; ya existen y estan
verificadas con corridas reales desde el 2026-08-01, ver STATUS.md).

De los 17 tipos de PTM del pipeline, solo 9 tienen un modulo de Fase A
implementado, repartidos en 3 clases mutuamente excluyentes:

- **Clase 1** (``pyrosetta_ptm_patch.SUPPORTED_PTM_TYPES``, 5 tipos:
  phosphorylation, acetylation, hydroxylation, gamma_carboxyglutamic_acid,
  lys_methylation): parche quimico nativo de Rosetta + relax de vecindario +
  estimacion de ddG (Extension 3, ``ddg_estimate.py``) comparando la pose WT
  relajada contra la parcheada relajada.
- **Clase 2** (``pyrosetta_glycan_patch.SUPPORTED_PTM_TYPES``, 2 tipos:
  n_linked_glycosylation, o_linked_glycosylation): adjunta el nucleo
  biosintetico conservado del glicano + refina su conformacion, con
  corroboracion opcional de GlyGen si se da un accession UniProt.
- **Clase 3** (``ubiquitin_sumo.SUPPORTED_PTM_TYPES``, 2 tipos:
  ubiquitination, sumoylation): conjugacion isopeptidica real via
  ``UBQ_GTPaseMover``.

Los otros 8 tipos (malonylation, arg_methylation, crotonylation,
succinylation, glutathionylation, s_nitrosylation, glutarylation,
citrullination) NO tienen ningun modulo de Fase A -- requeririan construir un
residuo no-canonico propio (params file real de cheminformatica), tarea NO
implementada, NO fabricada (ver docstring de ``pyrosetta_ptm_patch.py``).
``run_fase_a_for_site`` devuelve ``estado="sin_soporte_fase_a"`` para estos
sin intentar nada (ni siquiera inicializa PyRosetta).

Cada llamada a esta funcion asume un proceso PyRosetta fresco (invocado
siempre por subprocess, un sitio por proceso -- ver
``src/engines/_fase_a_runner.py``): las 3 clases inicializan PyRosetta con
flags DISTINTAS e incompatibles entre si (confirmado en los 3 modulos
reusados: ``pyrosetta.init()`` solo aplica las flags de su PRIMERA llamada
por proceso, llamadas repetidas con otras flags no tienen efecto real) -- por
eso esta funcion solo inicializa PyRosetta UNA vez, con las flags de la clase
que corresponde al ``ptm_type`` de ESTE sitio, nunca mezclando clases dentro
del mismo proceso.
"""

import statistics
from pathlib import Path
from typing import Optional

from src.config.settings import Settings
from src.structural import pyrosetta_glycan_patch, pyrosetta_ptm_patch, ubiquitin_sumo
from src.structural import ddg_estimate as ddg_estimate_module

CLASS1_TYPES = pyrosetta_ptm_patch.SUPPORTED_PTM_TYPES
CLASS2_TYPES = pyrosetta_glycan_patch.SUPPORTED_PTM_TYPES
CLASS3_TYPES = ubiquitin_sumo.SUPPORTED_PTM_TYPES

SUPPORTED_PTM_TYPES = frozenset(CLASS1_TYPES | CLASS2_TYPES | CLASS3_TYPES)

# Settings.FASE_A_SUPPORTED_PTM_TYPES es una copia sin pyrosetta de esta misma
# lista (para que ptm_annotation.py, sin PyRosetta, pueda seleccionar
# candidatos) -- falla alto en tiempo de import si alguna vez divergen, en vez
# de dejar que ptm_annotation.py seleccione en silencio un tipo sin soporte
# real (o que este modulo soporte un tipo que la seleccion nunca ve).
if SUPPORTED_PTM_TYPES != frozenset(Settings.FASE_A_SUPPORTED_PTM_TYPES):
    raise RuntimeError(
        "fase_a_dispatch.SUPPORTED_PTM_TYPES "
        f"({sorted(SUPPORTED_PTM_TYPES)}) no coincide con "
        f"Settings.FASE_A_SUPPORTED_PTM_TYPES ({sorted(Settings.FASE_A_SUPPORTED_PTM_TYPES)}) "
        "-- actualiza ambas listas juntas (ver docstring de Settings.FASE_A_SUPPORTED_PTM_TYPES)."
    )


def _run_class1(
    pdb_path: Path, position: int, ptm_type: str, out_pdb: Path,
    radius: float, nstruct: int,
) -> dict:
    """Clase 1: parche nativo + relax de vecindario + ddG (Extension 3)."""
    pyrosetta_ptm_patch.init_pyrosetta()

    ddg, wt_score, mut_score, wt_scores, mut_scores = ddg_estimate_module.estimate_ddg(
        pdb_path, position, ptm_type, radius=radius, nstruct=nstruct,
    )
    # Desviacion estandar entre las 'nstruct' trayectorias independientes de
    # relax por estado (WT y parcheado) -- ya se calculaba en la CLI standalone
    # de ddg_estimate.py (statistics.pstdev), pero antes de 2026-08-03 se
    # descartaba aqui sin llegar nunca al reporte final. ddg_std propaga el
    # error (suma en cuadratura, trayectorias WT/mutante independientes) para
    # que el reporte no muestre un ddG mas seguro de lo que realmente es.
    wt_score_std = statistics.pstdev(wt_scores) if len(wt_scores) > 1 else 0.0
    mut_score_std = statistics.pstdev(mut_scores) if len(mut_scores) > 1 else 0.0
    ddg_std = (wt_score_std ** 2 + mut_score_std ** 2) ** 0.5

    # Vuelca ademas una pose parcheada+relajada real (un unico relax, no
    # nstruct completo) como artefacto estructural visible del sitio -- el
    # ddG ya usa internamente el minimo de nstruct trayectorias, esta pose
    # es solo para tener un PDB de salida representativo.
    pose = pyrosetta_ptm_patch.load_pose(pdb_path)
    pyrosetta_ptm_patch.apply_ptm_patch(pose, position, ptm_type)
    pyrosetta_ptm_patch.relax_neighborhood(pose, position, radius=radius)
    pose.dump_pdb(str(out_pdb))

    return {
        **Settings.FASE_A_RESULT_TEMPLATE,
        "estado": "modelado",
        "clase": "class1_patch_ddg",
        "ddg": ddg,
        "ddg_std": ddg_std,
        "wt_score": wt_score,
        "wt_score_std": wt_score_std,
        "mut_score": mut_score,
        "mut_score_std": mut_score_std,
        "output_pdb": str(out_pdb),
    }


def _run_class2(
    pdb_path: Path, position: int, ptm_type: str, out_pdb: Path,
    refine_rounds: int, uniprot_accession: Optional[str],
) -> dict:
    """Clase 2: adjunta + refina glicano, corroboracion opcional de GlyGen."""
    glygen_evidencia = None
    if uniprot_accession:
        glygen_evidencia = pyrosetta_glycan_patch.check_glygen_evidence(
            uniprot_accession, position, ptm_type
        )

    pyrosetta_glycan_patch.init_pyrosetta()
    from pyrosetta import pose_from_pdb

    pose = pose_from_pdb(str(pdb_path))
    pyrosetta_glycan_patch.attach_glycan(pose, position, ptm_type)
    pyrosetta_glycan_patch.refine_glycan(pose, rounds=refine_rounds)
    pose.dump_pdb(str(out_pdb))

    return {
        **Settings.FASE_A_RESULT_TEMPLATE,
        "estado": "modelado",
        "clase": "class2_glycan",
        "glycan_tree": pyrosetta_glycan_patch.GLYCAN_TREE_BY_TYPE[ptm_type],
        "glygen_evidencia": glygen_evidencia,
        "output_pdb": str(out_pdb),
    }


def _run_class3(
    pdb_path: Path, position: int, ptm_type: str, out_pdb: Path,
    refine_cycles: int, refine_repack_cycles: int,
) -> dict:
    """Clase 3: conjugacion isopeptidica real (ubiquitinacion/sumoilacion)."""
    ubiquitin_sumo.init_pyrosetta(
        refine_cycles=refine_cycles, refine_repack_cycles=refine_repack_cycles
    )
    from pyrosetta import pose_from_pdb

    pose = pose_from_pdb(str(pdb_path))
    pose, metrics = ubiquitin_sumo.conjugate(pose, position, ptm_type)
    pose.dump_pdb(str(out_pdb))

    return {
        **Settings.FASE_A_RESULT_TEMPLATE,
        "estado": "modelado",
        "clase": "class3_conjugation",
        "conjugation_metrics": metrics,
        "output_pdb": str(out_pdb),
    }


def run_fase_a_for_site(
    pdb_path: Path,
    position: int,
    ptm_type: str,
    out_pdb: Path,
    *,
    uniprot_accession: Optional[str] = None,
    radius: float = 6.0,
    nstruct: int = 3,
    refine_rounds: int = 1,
    refine_cycles: int = 3,
    refine_repack_cycles: int = 3,
) -> dict:
    """Rutea un sitio (``pdb_path``, ``position``, ``ptm_type``) a su modulo de Fase A.

    Devuelve un dict con, como minimo, ``estado`` (``"modelado"``,
    ``"sin_soporte_fase_a"`` o ``"error"``) y ``clase``. Cualquier excepcion
    real de los modulos reusados (``ResidueMismatchError``,
    ``ConjugationFilterError``, etc.) se captura aqui y se traduce a
    ``estado="error"`` con el mensaje en ``error`` -- un sitio individual que
    falla nunca debe tumbar el resto del barrido (ver
    ``src/engines/_fase_a_runner.py``, que llama esta funcion una vez por
    sitio en su propio proceso).
    """
    pdb_path = Path(pdb_path)
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    if ptm_type not in SUPPORTED_PTM_TYPES:
        return {**Settings.FASE_A_RESULT_TEMPLATE, "estado": "sin_soporte_fase_a"}

    try:
        if ptm_type in CLASS1_TYPES:
            return _run_class1(pdb_path, position, ptm_type, out_pdb, radius, nstruct)
        if ptm_type in CLASS2_TYPES:
            return _run_class2(pdb_path, position, ptm_type, out_pdb, refine_rounds, uniprot_accession)
        return _run_class3(pdb_path, position, ptm_type, out_pdb, refine_cycles, refine_repack_cycles)
    except Exception as exc:  # noqa: BLE001 -- un sitio individual no debe tumbar el barrido completo
        return {**Settings.FASE_A_RESULT_TEMPLATE, "estado": "error", "error": f"{type(exc).__name__}: {exc}"}
