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
produccion). Conectado al pipeline principal 2026-08-03 via
``src/structural/fase_a_dispatch.py`` (ver su docstring) -- ya no es
solo un script suelto, aunque sigue siendo invocable de forma independiente
con su propio CLI.

## Robustez de produccion (2026-07-28, actualizado): ``nstruct`` en vez de un solo relax

``FastRelax`` es un protocolo estocastico (Monte Carlo): una sola corrida
puede caer en un minimo local subobtimo, dando un score mas alto de lo que
en realidad seria el estado relajado. Practica estandar de los protocolos
de ddG de Rosetta (``cartesian_ddg``/``ddg_monomer``, mismo criterio citado
en la decision del 28-07 sobre el workaround de FastRelax cartesiano):
correr ``nstruct`` trayectorias independientes por estado (WT y parcheado)
y quedarse con la de MENOR energia por estado -- no promediar los scores
directamente, que mezclaria trayectorias mal convergidas con las buenas.
Con ``nstruct=1`` (comportamiento anterior) el resultado depende del azar
de una sola trayectoria; con ``nstruct=3`` (default nuevo) el ddG es mucho
mas estable frente a re-ejecuciones, a costa de proporcionalmente mas
tiempo de computo (cada relax individual ya tarda ~1 min sobre una
proteina de tamano medio, ver STATUS.md).
"""

import argparse
import statistics
import sys
from pathlib import Path

from src.structural.pyrosetta_ptm_patch import (
    SUPPORTED_PTM_TYPES,
    apply_ptm_patch,
    init_pyrosetta,
    load_pose,
    relax_neighborhood,
)


def _best_of_n_relax(
    pdb_path: Path, position: int, ptm_type: str, scorefxn, nstruct: int, radius: float, max_iter: int
):
    """Corre ``nstruct`` trayectorias independientes de parcheo+relax y devuelve ``(min_score, scores)``.

    ``ptm_type=None`` corre el estado WT (sin parche). Cada trayectoria
    parte de una pose fresca cargada desde ``pdb_path`` (no reusa/perturba
    la pose de la trayectoria anterior), para que sean independientes de
    verdad y no una cadena de Markov correlacionada.
    """
    scores = []
    for _ in range(nstruct):
        pose = load_pose(pdb_path)
        if ptm_type is not None:
            apply_ptm_patch(pose, position, ptm_type)
        relax_neighborhood(pose, position, radius=radius, scorefxn=scorefxn, max_iter=max_iter)
        scores.append(scorefxn(pose))
    return min(scores), scores


def estimate_ddg(
    pdb_path: Path, position: int, ptm_type: str, radius: float = 6.0, max_iter: int = 200, nstruct: int = 3
):
    """Devuelve ``(ddg, wt_score, mut_score, wt_scores, mut_scores)`` para ``ptm_type`` en ``position``.

    ``wt_score``/``mut_score`` son el minimo de ``nstruct`` trayectorias
    independientes por estado (ver docstring del modulo); ``wt_scores``/
    ``mut_scores`` son las listas completas, utiles para reportar
    dispersion (desviacion estandar) como medida de confianza del ddG.
    """
    import pyrosetta

    scorefxn = pyrosetta.create_score_function("ref2015_cart")

    wt_score, wt_scores = _best_of_n_relax(
        pdb_path, position, None, scorefxn, nstruct, radius, max_iter
    )
    mut_score, mut_scores = _best_of_n_relax(
        pdb_path, position, ptm_type, scorefxn, nstruct, radius, max_iter
    )

    return mut_score - wt_score, wt_score, mut_score, wt_scores, mut_scores


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Estima el ddG (ref2015_cart) de un sitio PTM real (Extension 3)."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ptm-type", required=True, choices=sorted(SUPPORTED_PTM_TYPES))
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument(
        "--nstruct", type=int, default=3,
        help="Trayectorias independientes de relax por estado (WT/parcheado); se usa la de menor "
        "energia por estado, practica estandar de los protocolos ddG de Rosetta. 1 = comportamiento "
        "original (una sola relajacion, mas rapido pero mas ruidoso).",
    )
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    init_pyrosetta()
    ddg, wt_score, mut_score, wt_scores, mut_scores = estimate_ddg(
        Path(args.pdb_path), args.position, args.ptm_type,
        radius=args.radius, nstruct=args.nstruct,
    )
    wt_std = statistics.pstdev(wt_scores) if len(wt_scores) > 1 else 0.0
    mut_std = statistics.pstdev(mut_scores) if len(mut_scores) > 1 else 0.0

    with open(args.out_csv, "w") as f:
        f.write("position,ptm_type,nstruct,wt_score,mut_score,ddg,wt_score_std,mut_score_std\n")
        f.write(
            f"{args.position},{args.ptm_type},{args.nstruct},{wt_score},{mut_score},{ddg},"
            f"{wt_std},{mut_std}\n"
        )

    print(
        f"ddG({args.ptm_type} @ {args.position}) = {ddg:.4f} (ref2015_cart, nstruct={args.nstruct}, "
        f"wt_std={wt_std:.2f}, mut_std={mut_std:.2f}) -> '{args.out_csv}'"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
