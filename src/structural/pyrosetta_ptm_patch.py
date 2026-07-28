#!/usr/bin/env python
"""Fase A clase 1 (parches quimicos simples): utilidad compartida de parcheo PyRosetta.

Standalone -- ver docstring de ``src/structural/__init__.py``: requiere
``pyrosetta``, solo presente en el conda env ``deepptmpred``, nunca se
importa desde el paquete ``src`` principal.

Cobertura real verificada 2026-07-28 en esta maquina (no asumida): de los 17
tipos de PTM del pipeline (``Settings.DEEPPTMPRED_PTM_TYPES``), Rosetta trae
un ``VariantType``/patch nativo listo para usar SOLO para 5 -- confirmado por
DOS vias independientes:

1. ``ls $PYROSETTA_DB/chemical/residue_type_sets/fa_standard/patches/*.txt``
   (listado real del PyRosetta instalado, no documentacion): existen
   ``ser_phosphorylated.txt``/``thr_phosphorylated.txt``/
   ``tyr_phosphorylated.txt``, ``lys_acetylated.txt``,
   ``pro_hydroxylated_case1.txt``/``case2.txt``, ``glu_carboxylated.txt``,
   ``lys_monomethylated.txt``/``dimethylated``/``trimethylated``. NO existen
   patches para malonylation, arg_methylation, crotonylation, succinylation,
   glutathionylation, s_nitrosylation, glutarylation ni citrullination
   (grep de esos terminos sobre el directorio completo: cero resultados).
2. Introspeccion en caliente de ``pyrosetta.rosetta.core.chemical.VariantType``
   (enum real del build instalado): mismo resultado -- existen
   ``PHOSPHORYLATION``, ``ACETYLATION``, ``HYDROXYLATION``/``HYDROXYLATION1``/
   ``HYDROXYLATION2``, ``CARBOXYLATION``, ``METHYLATION``/``DIMETHYLATION``/
   ``TRIMETHYLATION``; ninguna variante para los otros 8 tipos.

Esto CORRIGE la estimacion "12/17 tipos" de la decision de arquitectura
2026-07-28 (vault, ``2026-07-28-investigacion-ddg-y-modelado-estructural.md``):
esa cifra asumia mayor cobertura nativa sin haber inspeccionado el enum/
patches reales todavia. La cifra verificada real es 5/17 con patch nativo
listo para usar. Los otros 8 tipos de "clase 1" (fuera de glicosilacion,
clase 2, y ubiquitinacion/sumoilacion, clase 3) requeririan construir un
residuo no-canonico propio (params file con topologia, cargas parciales e
icoor -- tarea de cheminformatica real, no un simple parche) -- NO
implementado aqui, NO fabricado. Queda como item explicito pendiente, ver
STATUS.md.

Aplicado y verificado con una ejecucion real (2026-07-28): parche de
acetilacion sobre un residuo LYS real de ``AF-P10636-F1-model_v4.pdb``
(Tau) via ``add_variant_type_to_pose_residue`` -- el residuo pasa de
``LYS`` a ``LYS:acetylated``, confirmado leyendo ``pose.residue(i).name()``
antes/despues.

``relax_neighborhood`` implementa el mismo workaround documentado en
RosettaCommons para ddG/impacto estructural sobre proteinas grandes:
restringe FastRelax (backbone+chi moviles, packing) a los residuos dentro de
un radio del sitio parcheado en vez de relajar la proteina completa --
evita un coste computacional prohibitivo sin perder la respuesta
estructural real al parche. Reusada por ``ddg_estimate.py`` (Extension 3),
mismo criterio de reuso que fijo el orden de implementacion decidido el
2026-07-28.
"""

import argparse
import sys
from pathlib import Path

# Tipo canonico (mismos nombres que Settings.DEEPPTMPRED_PTM_TYPES) -> nombre
# del atributo real en pyrosetta.rosetta.core.chemical.VariantType (verificado
# arriba). El grado de metilacion (mono/di/tri) no lo distinguen los motores
# (DeepMVP/DeepPTMPred predicen 'lys_methylation' sin especificar) -> se usa
# mono por defecto, documentado explicitamente aqui, no una eleccion oculta.
PTM_VARIANT_MAP = {
    "phosphorylation": "PHOSPHORYLATION",
    "acetylation": "ACETYLATION",
    "hydroxylation": "HYDROXYLATION",
    "gamma_carboxyglutamic_acid": "CARBOXYLATION",
    "lys_methylation": "METHYLATION",
}

SUPPORTED_PTM_TYPES = frozenset(PTM_VARIANT_MAP)

# Residuo diana esperado por tipo, verificado leyendo
# DeepPTMPred/pred/train_PTM/predict.py::PredictConfig.ptm_aa_map (mismo
# mapeo que ya usa el pipeline principal) -- valida que la posicion pedida
# sea quimicamente coherente con el tipo de PTM antes de parchear.
PTM_TARGET_RESIDUE = {
    "phosphorylation": {"S", "T", "Y"},
    "acetylation": {"K"},
    "hydroxylation": {"P"},
    "gamma_carboxyglutamic_acid": {"E"},
    "lys_methylation": {"K"},
}


class UnsupportedPTMPatchError(Exception):
    """``ptm_type`` sin VariantType nativo de Rosetta verificado (ver docstring del modulo)."""


class ResidueMismatchError(Exception):
    """El residuo real en ``position`` no es el esperado para ``ptm_type``."""


def init_pyrosetta(extra_flags: str = "-mute all -ex1 -ex2aro") -> None:
    """Inicializa PyRosetta una sola vez por proceso (idempotente en la practica: PyRosetta
    ignora llamadas repetidas a ``init`` dentro del mismo proceso)."""
    import pyrosetta

    pyrosetta.init(extra_flags)


def load_pose(pdb_path: Path):
    from pyrosetta import pose_from_pdb

    return pose_from_pdb(str(pdb_path))


def apply_ptm_patch(pose, position: int, ptm_type: str):
    """Aplica el ``VariantType`` de Rosetta correspondiente a ``ptm_type`` en ``position`` (1-based).

    Muta ``pose`` in-place y lo devuelve por conveniencia.

    Raises:
        UnsupportedPTMPatchError: Si ``ptm_type`` no esta en
            :data:`PTM_VARIANT_MAP` (ver docstring del modulo para el
            detalle completo de que tipos SI/NO estan cubiertos).
        ResidueMismatchError: Si el residuo real en ``position`` no coincide
            con el residuo diana esperado para ``ptm_type``.
    """
    from pyrosetta.rosetta.core.chemical import VariantType
    from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue

    if ptm_type not in PTM_VARIANT_MAP:
        raise UnsupportedPTMPatchError(
            f"'{ptm_type}' no tiene un VariantType nativo de Rosetta verificado "
            f"(soportados: {sorted(SUPPORTED_PTM_TYPES)}). Ver docstring de este "
            "modulo / STATUS.md para el detalle de por que los demas tipos "
            "requieren construir un residuo no-canonico propio."
        )

    actual_residue = pose.residue(position).name1()
    expected = PTM_TARGET_RESIDUE[ptm_type]
    if actual_residue not in expected:
        raise ResidueMismatchError(
            f"Posicion {position} es '{actual_residue}', pero '{ptm_type}' espera uno de "
            f"{sorted(expected)}. Verifica que la numeracion coincide con la del PDB de entrada."
        )

    variant = getattr(VariantType, PTM_VARIANT_MAP[ptm_type])
    add_variant_type_to_pose_residue(pose, variant, position)
    return pose


def relax_neighborhood(pose, position: int, radius: float = 6.0, scorefxn=None, max_iter: int = 200):
    """FastRelax (cartesiano, 1 repeat) restringido al vecindario de ``position``.

    Backbone y chi moviles solo dentro de ``radius`` angstrom del sitio
    (medido residuo-residuo, incluye el propio sitio); el resto de la
    proteina se mantiene fijo (``PreventRepackingRLT`` + movemap sin bb/chi
    fuera del vecindario). Verificado con una ejecucion real 2026-07-28
    sobre ``AF-P10636-F1-model_v4.pdb`` (758 residuos): converge sin error
    con 1 repeat / max_iter=200.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.select.residue_selector import (
        NeighborhoodResidueSelector, ResidueIndexSelector, NotResidueSelector,
    )
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, PreventRepackingRLT, OperateOnResidueSubset,
    )
    from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
    from pyrosetta.rosetta.protocols.relax import FastRelax

    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015_cart")

    site_selector = ResidueIndexSelector(str(position))
    neighborhood = NeighborhoodResidueSelector(site_selector, radius, True)

    mmf = MoveMapFactory()
    mmf.add_bb_action(move_map_action.mm_enable, neighborhood)
    mmf.add_chi_action(move_map_action.mm_enable, neighborhood)

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(neighborhood)))

    relax = FastRelax(1)
    relax.set_scorefxn(scorefxn)
    relax.cartesian(True)
    relax.set_movemap_factory(mmf)
    relax.set_task_factory(tf)
    relax.max_iter(max_iter)
    relax.apply(pose)
    return pose


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aplica un parche PTM real (Fase A clase 1) y relaja el vecindario local."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int, help="Posicion 1-based (numeracion del PDB de entrada).")
    parser.add_argument("--ptm-type", required=True, choices=sorted(SUPPORTED_PTM_TYPES))
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument("--skip-relax", action="store_true", help="Solo aplica el parche, sin FastRelax del vecindario.")
    args = parser.parse_args()

    init_pyrosetta()
    pose = load_pose(Path(args.pdb_path))
    apply_ptm_patch(pose, args.position, args.ptm_type)
    if not args.skip_relax:
        relax_neighborhood(pose, args.position, radius=args.radius)
    pose.dump_pdb(args.out_pdb)
    print(f"Parche '{args.ptm_type}' aplicado en posicion {args.position} -> '{args.out_pdb}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
