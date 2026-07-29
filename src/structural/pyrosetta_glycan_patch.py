#!/usr/bin/env python
"""Fase A clase 2 (glicosilacion): adjunta un arbol de glicano real via PyRosetta.

Standalone -- ver docstring de ``src/structural/__init__.py``.

## Herramienta

`GlycanTreeModeler` (bloqueante de documentacion del 28-07, 403 en
`docs.rosettacommons.org`) se confirmo real e instalado por introspeccion
directa (``pyrosetta.rosetta.protocols.carbohydrates.GlycanTreeModeler``),
sin necesitar la documentacion bloqueada. Pero ese Mover REFINA la
conformacion de un glicano YA adjunto -- para adjuntarlo primero hace falta
``SimpleGlycosylateMover`` (tambien real, mismo namespace), que si
construye el arbol de glicano desde cero sobre un residuo dado.

## Que arbol de glicano se adjunta

Ni DeepMVP ni DeepPTMPred predicen la COMPOSICION del glicano (glicoforma),
solo el SITIO -- la composicion real depende de la maquinaria de
glicosilacion de la celula/tejido, no es inferible desde la secuencia sola.
Se usa el nucleo biosintetico conservado como default DOCUMENTADO, no una
prediccion:
- N-linked: ``N-glycan_core`` (Man3GlcNAc2, el nucleo pentasacarido
  presente en TODO N-glicano maduro antes de cualquier ramificacion/
  procesamiento especifico de tejido -- unico default universalmente
  defendible).
- O-linked: ``core_1_O-glycan`` (Galβ1-3GalNAc, "antigeno T", el core O-glicano
  mas comun en mamiferos, pero NO universal como el core N-glicano --
  existen 8 cores O-glicano distintos en `common_glycans/`, la eleccion real
  depende del tejido/enzima, aqui se usa el mas frecuente como default
  razonable, no una prediccion).
Ambos strings IUPAC-condensados verificados leyendo
``database/chemical/carbohydrates/common_glycans/*.iupac`` reales del
PyRosetta instalado, no inventados.

Verificado con una ejecucion real 2026-07-28: adjunto N-glycan_core en
ASN484 (sequon real NAT) y core_1_O-glycan en THR52 de
``AF-P10636-F1-model_v4.pdb`` (Tau) -- 5 y 2 residuos de azucar anadidos
respectivamente, refinamiento con `GlycanTreeModeler` (1 round) corre sin
error (~64s para el caso N-glycan_core).

## Corroboracion GlyGen (2026-07-29, opcional via --uniprot-accession)

Mejora identificada 2026-07-28, implementada aqui: antes de adjuntar el
nucleo generico, se puede consultar GlyGen (``src.structural.glygen_client``,
API real, no requiere PyRosetta) para saber si el sitio ya tiene evidencia
EXPERIMENTAL reportada (no solo predicha por otra herramienta). Esto NO
cambia el glicano que se adjunta (el nucleo biosintetico generico sigue
siendo lo unico que PyRosetta puede construir directamente aqui -- el
``glytoucan_ac`` que reporta GlyGen identifica un glicano especifico que
requeriria traducirlo a un arbol IUPAC construible, un problema de mapeo
real no resuelto por esta mejora) -- es puramente informativo, para que quien
interprete el resultado sepa si el sitio esta ya caracterizado
experimentalmente o es una posicion sin evidencia externa conocida. Fallos
de red degradan a un aviso, nunca abortan el flujo principal (el adjuntado
de glicano no depende de la disponibilidad de GlyGen).
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

GLYCAN_TREE_BY_TYPE = {
    "n_linked_glycosylation": "N-glycan_core",
    "o_linked_glycosylation": "core_1_O-glycan",
}

TARGET_RESIDUE = {
    "n_linked_glycosylation": {"N"},
    "o_linked_glycosylation": {"S", "T"},
}

SUPPORTED_PTM_TYPES = frozenset(GLYCAN_TREE_BY_TYPE)


class UnsupportedGlycanTypeError(Exception):
    """``ptm_type`` no es n_linked_glycosylation ni o_linked_glycosylation."""


class ResidueMismatchError(Exception):
    """El residuo real en ``position`` no es el esperado para ``ptm_type``."""


def init_pyrosetta() -> None:
    import pyrosetta

    pyrosetta.init(
        "-mute all -include_sugars -alternate_3_letter_codes pdb_sugar "
        "-write_pdb_link_records true"
    )


def attach_glycan(pose, position: int, ptm_type: str):
    """Adjunta el nucleo de glicano por defecto de ``ptm_type`` en ``position`` (1-based).

    Muta ``pose`` in-place y lo devuelve por conveniencia.
    """
    from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover

    if ptm_type not in GLYCAN_TREE_BY_TYPE:
        raise UnsupportedGlycanTypeError(
            f"'{ptm_type}' no es un tipo de glicosilacion soportado (soportados: "
            f"{sorted(SUPPORTED_PTM_TYPES)})."
        )

    actual_residue = pose.residue(position).name1()
    expected = TARGET_RESIDUE[ptm_type]
    if actual_residue not in expected:
        raise ResidueMismatchError(
            f"Posicion {position} es '{actual_residue}', pero '{ptm_type}' espera uno de "
            f"{sorted(expected)}."
        )

    mover = SimpleGlycosylateMover()
    mover.set_position(position)
    mover.set_glycosylation(GLYCAN_TREE_BY_TYPE[ptm_type])
    mover.apply(pose)
    return pose


def refine_glycan(pose, rounds: int = 1, scorefxn=None):
    """Refina la conformacion de TODOS los glicanos ya adjuntos en ``pose`` via ``GlycanTreeModeler``.

    No restringe a un solo arbol (``GlycanTreeModeler`` detecta automaticamente
    los residuos de carbohidrato presentes en la pose) -- si se llama con
    varios sitios ya adjuntos, refina todos a la vez.
    """
    import pyrosetta
    from pyrosetta.rosetta.protocols.carbohydrates import GlycanTreeModeler

    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")

    modeler = GlycanTreeModeler()
    modeler.set_rounds(rounds)
    modeler.set_refine(True)
    modeler.set_scorefxn(scorefxn)
    modeler.apply(pose)
    return pose


def check_glygen_evidence(uniprot_accession: str, position: int, ptm_type: str) -> Optional[str]:
    """Consulta GlyGen y devuelve un mensaje informativo, o ``None`` si no hay corroboracion.

    Puramente informativo (ver docstring del modulo) -- nunca lanza: un
    fallo de red o accession no reconocido se reporta como un mensaje de
    aviso, no como una excepcion, para que el flujo principal (adjuntar el
    nucleo generico) nunca dependa de que GlyGen este disponible.
    """
    from src.structural.glygen_client import GlyGenLookupError, lookup_site

    try:
        site = lookup_site(uniprot_accession, position, ptm_type)
    except GlyGenLookupError as exc:
        return f"GlyGen no disponible ({exc}) -- se usa el nucleo generico sin corroboracion."

    if site is None:
        return (
            f"GlyGen no reporta ningun sitio de '{ptm_type}' en la posicion {position} de "
            f"'{uniprot_accession}' -- sin evidencia externa conocida, se usa el nucleo generico."
        )

    glytoucan_ac = site.get("glytoucan_ac")
    category = site.get("site_category", "desconocida")
    if glytoucan_ac:
        return (
            f"GlyGen SI reporta evidencia experimental en este sitio (categoria='{category}', "
            f"glycan={glytoucan_ac}) -- el nucleo generico adjuntado aqui es un placeholder "
            "documentado, no ese glicano especifico (ver docstring del modulo)."
        )
    return (
        f"GlyGen reporta este sitio con categoria='{category}' (sin glicano especifico "
        "identificado) -- se usa el nucleo generico."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Adjunta (y opcionalmente refina) un glicano real (Fase A clase 2)."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ptm-type", required=True, choices=sorted(SUPPORTED_PTM_TYPES))
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument("--skip-refine", action="store_true")
    parser.add_argument("--refine-rounds", type=int, default=1)
    parser.add_argument(
        "--uniprot-accession", default=None,
        help="Accession UniProt canonico (sin sufijo de isoforma) para consultar GlyGen como "
             "corroboracion informativa antes de adjuntar el nucleo generico. Opcional: si se "
             "omite, no se consulta GlyGen (comportamiento identico al de antes de esta mejora).",
    )
    args = parser.parse_args()

    if args.uniprot_accession:
        message = check_glygen_evidence(args.uniprot_accession, args.position, args.ptm_type)
        print(f"[GlyGen] {message}")

    init_pyrosetta()
    from pyrosetta import pose_from_pdb

    pose = pose_from_pdb(args.pdb_path)
    attach_glycan(pose, args.position, args.ptm_type)
    if not args.skip_refine:
        refine_glycan(pose, rounds=args.refine_rounds)
    pose.dump_pdb(args.out_pdb)
    print(
        f"Glicano '{GLYCAN_TREE_BY_TYPE[args.ptm_type]}' adjunto en posicion {args.position} "
        f"-> '{args.out_pdb}'"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
