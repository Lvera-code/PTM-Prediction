#!/usr/bin/env python
"""Extension 3: estimacion de impacto en estabilidad (ddG) de una PTM real via PyRosetta.

Standalone -- ver docstring de ``src/structural/__init__.py``. Reusa
``pyrosetta_ptm_patch.py`` (parcheo + relax de vecindario), orden de
implementacion decidido 2026-07-28 (utilidad compartida primero, Extension 3
se apoya en ella). Cubre unicamente los 5 tipos con VariantType nativo
verificado (``pyrosetta_ptm_patch.SUPPORTED_PTM_TYPES`` -- ver su docstring
para el detalle de cobertura real vs. la estimacion original "12/17").

Protocolo (workaround real documentado en RosettaCommons para ddG
cartesiano, mismo criterio citado en la decision 2026-07-28): compara el
score (``ref2015_cart``) de la pose WT relajada localmente contra la pose
parcheada relajada localmente, ambas con el MISMO protocolo de relax
(``pyrosetta_ptm_patch.relax_neighborhood``) para que la diferencia refleje
el efecto del parche y no ruido de muestreo entre protocolos distintos.
``ddG = score(parcheada) - score(WT)``: positivo = parche desestabiliza,
negativo = estabiliza (convencion estandar Rosetta ref2015).

Verificado con una ejecucion real 2026-07-28: acetilacion de una lisina real
en ``AF-P10636-F1-model_v4.pdb`` (Tau), ver STATUS.md para el valor
obtenido y el tiempo de computo real (relevante para decidir timeouts en
produccion -- el pipeline principal NO invoca este script todavia, decision
2026-07-27: D es un filtro de responsabilidad unica sin rutas a Extension 3).
"""

import argparse
import sys
from pathlib import Path

from pyrosetta_ptm_patch import (
    SUPPORTED_PTM_TYPES,
    apply_ptm_patch,
    init_pyrosetta,
    load_pose,
    relax_neighborhood,
)


def estimate_ddg(pdb_path: Path, position: int, ptm_type: str, radius: float = 6.0, max_iter: int = 200):
    """Devuelve ``(ddg, wt_score, mut_score)`` para el parche ``ptm_type`` en ``position``."""
    import pyrosetta

    scorefxn = pyrosetta.create_score_function("ref2015_cart")

    wt_pose = load_pose(pdb_path)
    relax_neighborhood(wt_pose, position, radius=radius, scorefxn=scorefxn, max_iter=max_iter)
    wt_score = scorefxn(wt_pose)

    mut_pose = load_pose(pdb_path)
    apply_ptm_patch(mut_pose, position, ptm_type)
    relax_neighborhood(mut_pose, position, radius=radius, scorefxn=scorefxn, max_iter=max_iter)
    mut_score = scorefxn(mut_pose)

    return mut_score - wt_score, wt_score, mut_score


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Estima el ddG (ref2015_cart) de un sitio PTM real (Extension 3)."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ptm-type", required=True, choices=sorted(SUPPORTED_PTM_TYPES))
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    init_pyrosetta()
    ddg, wt_score, mut_score = estimate_ddg(
        Path(args.pdb_path), args.position, args.ptm_type, radius=args.radius
    )

    with open(args.out_csv, "w") as f:
        f.write("position,ptm_type,wt_score,mut_score,ddg\n")
        f.write(f"{args.position},{args.ptm_type},{wt_score},{mut_score},{ddg}\n")

    print(f"ddG({args.ptm_type} @ {args.position}) = {ddg:.4f} (ref2015_cart) -> '{args.out_csv}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
