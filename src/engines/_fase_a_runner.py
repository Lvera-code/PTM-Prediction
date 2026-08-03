#!/usr/bin/env python
"""Runner standalone de Fase A: modela UN sitio PTM real via PyRosetta.

Invocado siempre por subprocess desde ``src/engines/fase_a_engine.py`` con el
interprete del conda env ``deepptmpred`` (``Settings.FASE_A_PYTHON_BIN`` --
mismo env que ya tiene PyRosetta instalado y verificado para DeepPTMPred,
ver STATUS.md), un sitio por proceso -- mismo patron que
``_deepptmpred_runner.py``/``_metoken_runner.py``: nunca importado
directamente por el proceso principal del pipeline (que no tiene PyRosetta
instalado).

Delega todo el enrutamiento real a ``src.structural.fase_a_dispatch`` (ver su
docstring para el detalle de las 3 clases de Fase A soportadas). Este script
solo se encarga de la interfaz CLI + serializar el resultado a JSON, mismo
rol que cumplen los ``main()`` de los 4 modulos de ``src/structural/`` para
uso manual, pero pensado para ser invocado programaticamente en un barrido de
muchos sitios (ver ``FaseAEngine``).
"""

import argparse
import json
import sys
from pathlib import Path

# src.structural.fase_a_dispatch hace imports absolutos ('from src.structural
# import ...') que requieren la raiz del repo en sys.path -- a diferencia de
# los runners de motores externos (DeepPTMPred/MeToken/StackGlyEmbed), este
# runner importa codigo PROPIO del proyecto, no un repo vendorizado, por eso
# necesita este insert explicito (subprocess.run no fija cwd para este
# runner, ver fase_a_engine.py).
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from src.config.settings import Settings  # noqa: E402
from src.structural.fase_a_dispatch import SUPPORTED_PTM_TYPES  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Modela un sitio PTM real (Fase A: parche+ddG / glicano / conjugacion)."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ptm-type", required=True)
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--uniprot-accession", default=None)
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--nstruct", type=int, default=3)
    parser.add_argument("--refine-rounds", type=int, default=1)
    parser.add_argument("--refine-cycles", type=int, default=3)
    parser.add_argument("--refine-repack-cycles", type=int, default=3)
    args = parser.parse_args()

    # Import diferido (despues del sys.path.insert de arriba, y solo cuando
    # hace falta): 'sin_soporte_fase_a' se resuelve sin tocar PyRosetta en
    # absoluto, igual que hace run_fase_a_for_site.
    from src.structural.fase_a_dispatch import run_fase_a_for_site

    if args.ptm_type not in SUPPORTED_PTM_TYPES:
        result = {**Settings.FASE_A_RESULT_TEMPLATE, "estado": "sin_soporte_fase_a"}
    else:
        result = run_fase_a_for_site(
            Path(args.pdb_path), args.position, args.ptm_type, Path(args.out_pdb),
            uniprot_accession=args.uniprot_accession, radius=args.radius, nstruct=args.nstruct,
            refine_rounds=args.refine_rounds, refine_cycles=args.refine_cycles,
            refine_repack_cycles=args.refine_repack_cycles,
        )

    result["position"] = args.position
    result["ptm_type"] = args.ptm_type

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Fase A ({args.ptm_type} @ {args.position}): estado='{result['estado']}' -> '{args.out_json}'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
