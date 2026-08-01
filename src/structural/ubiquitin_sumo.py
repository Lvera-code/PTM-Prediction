#!/usr/bin/env python
"""Fase A clase 3 (ubiquitinacion/sumoilacion): conjugacion real via PyRosetta.

Standalone -- ver docstring de ``src/structural/__init__.py``.

## Correccion de una conclusion previa (2026-08-01)

STATUS.md (seccion "Fase A clase 3", 2026-07-28) concluia que esto exigia
compilar Rosetta completo desde fuente, porque las apps publicadas
(``UBQ_E2_thioester``, serie ``UBQ_Gp``) son binarios C++. **Esa conclusion
era incorrecta**, verificado leyendo el .cc real de ambas piezas
(``github.com/RosettaCommons/rosetta``,
``source/src/protocols/chemically_conjugated_docking/UBQ_GTPaseMover.cc`` y
``source/src/apps/public/scenarios/chemically_conjugated_docking/UBQ_Gp_LYX-Cterm.cc``):
las apps son wrappers de ~90 lineas alrededor de
``protocols::chemically_conjugated_docking::UBQ_GTPaseMover``, un
``protocols::moves::Mover`` normal -- PyRosetta expone todos los Movers de
``protocols::`` automaticamente, sin compilar nada nuevo. Confirmado real
contra el wheel YA instalado (``deepptmpred``, PyRosetta 2026.30):
``pyrosetta.rosetta.protocols.chemically_conjugated_docking.UBQ_GTPaseMover``
existe e instancia sin error.

El patch de Rosetta que si estaba desactivado (``branching/C-terminal_conjugation.txt``,
"something is broken with this patch", ver seccion "Fase A clase 3" mas
abajo en STATUS.md) **no interviene en este camino** -- ``UBQ_GTPaseMover``
no usa ese patch. Usa en su lugar un ``ResidueType`` completo, ``LYX``
(lisina modificada para aceptar el enlace isopeptidico), que se activa con
una flag de init distinta (ver ``init_pyrosetta()`` abajo).

## Mecanismo real

``UBQ_GTPaseMover.initialize()`` (leido linea por linea del .cc real):
1. Toma el PRIMER residuo que devuelve el ``ResidueSelector`` configurado
   (aqui, ``make_index_selector(position)`` -- selecciona exactamente ese
   indice).
2. Carga el PDB de ubiquitina/SUMO (``set_UBQpdb``), borra su ULTIMO
   residuo y anade una glicina fresca en su lugar (irrelevante si el PDB de
   referencia ya termina en Gly nativa, como es nuestro caso).
3. Reemplaza el residuo objetivo en la proteina query por el ``ResidueType``
   ``LYX`` (lisina modificada) y forma el enlace isopeptidico real
   (``append_residue_by_bond`` con las conexiones ``lyx_connid=3``/
   ``ubq_connid=2`` verificadas en el .cc).
4. Corre un protocolo de refinamiento Monte Carlo completo (Small/Shear/
   SidechainMover + minimizacion + repacking) sobre la cola movil.

**Hallazgo real, no documentado por Rosetta**: la linea 207 del .cc
(``runtime_assert(GTPase.residue_type(GTPase_lys_).aa() == core::chemical::aa_lys)``)
esta COMENTADA -- Rosetta NO verifica que el residuo objetivo sea
realmente una lisina antes de forzarlo a ``LYX``. ``_validate_target()``
aqui abajo si lo hace (mismo patron que ``pyrosetta_glycan_patch.ResidueMismatchError``).

**Segundo hallazgo real**: ``initialize()`` (linea 206) hace
``GTPase.pdb_info()->pdb2pose(chain, GTPase_residue)`` sobre un indice que
YA es de numeracion de pose (el propio ``ResidueSelector`` devuelve indices
de pose, no de PDB) -- una doble conversion que solo da el resultado
correcto si la numeracion PDB coincide 1:1 con la numeracion de pose
(cierto para los PDBs de AlphaFold de una sola cadena que usa este
pipeline, verificado empiricamente: ``Tau`` residuo 24 -> ``GTPase_lys``
tras ``initialize()`` = 24). ``_validate_target()`` verifica esta
precondicion explicitamente en vez de asumirla -- si algun dia se usa un
PDB real con huecos/numeracion no secuencial, fallaria alto y claro en vez
de conjugar silenciosamente el residuo equivocado.

**Tercer hallazgo real (gotcha de integracion, no de correctud)**: fuera de
una corrida real via ``JobDistributor`` (que no usamos, invocamos el Mover
directo), ``protocols::jd2::JobDistributor::get_instance()->current_job()``
devuelve un job dummy ``EMPTY_JOB_use_jd2`` en vez de ``nullptr`` --
confirmado empiricamente que NO revienta, pero dos efectos secundarios
reales:
- ``initialize()`` escribe SIEMPRE ``starting_pose_EMPTY_JOB_use_jd2_0000.pdb``
  en el directorio de trabajo ACTUAL (verificado con una corrida real:
  aparecio en el cwd del proceso). ``conjugate()`` aqui abajo confina esto a
  un directorio temporal desechable, igual que otros modulos de
  ``src/structural/`` evitan efectos secundarios en el cwd del caller.
- ``analyze_and_filter()`` escribe las 5 metricas por-decoy (angulos de la
  cola movil: ``lysine_chi3``, ``lysine_chi4``, ``amide``, ``glycine_psi``,
  ``glycine_phi``) en ese job dummy, inaccesible desde Python. En vez de
  intentar leerlas de ahi, ``conjugate()`` las recalcula directamente sobre
  la pose de salida usando ``mover.get_atom_IDs()`` (expuesto, verificado
  real) + ``pose.atom_tree().torsion_angle(...)`` -- mismos 5 angulos, mismo
  valor, fuente independiente del job descartado.

## Verificado con corridas reales 2026-08-01

Contra ``DeepPTMPred/data/AF-P10636-F1-model_v4.pdb`` (Tau, LYS24, mismo
residuo ya usado en la Extension 3 para el ddG de acetilacion):
- ``initialize()`` solo (sin refinar): pose de salida 834 residuos
  (758 Tau + 76 ubiquitina), residuo en la posicion objetivo = ``LYX``.
- ``apply()`` completo (``refine_cycles=3``, valor de desarrollo rapido, ver
  mas abajo): ``MS_SUCCESS``, 16.4s, 834 residuos, 5 metricas de torsion
  recalculadas con valores reales no triviales.
- SUMO1 (``sumo1_reference.pdb``): mismo mecanismo, pose de salida 857
  residuos (758 + 99), residuo objetivo = ``LYX`` -- confirma que el mismo
  Mover sirve para ambos tipos sin cambios de codigo, solo cambiando el PDB
  de referencia (SUMO tambien termina en diglicina nativa expuesta tras
  maduracion, misma quimica isopeptidica que ubiquitina).

## PDBs de referencia empaquetados (``src/structural/data/``)

- ``ubiquitin_reference.pdb``: RCSB 1UBQ (Vijay-Kumar et al. 1987, 1.8A),
  cadena A, aguas eliminadas via ``gemmi``. Mismo PDB que usa el propio
  equipo de Rosetta en sus tests de integracion
  (``tests/integration/tests/UBQ_Gp_LYX-Cterm/options``:
  ``-chemically_conjugated_docking:UBQpdb ../UBQ_Gp_CYD-CYD/1UBQ.pdb.gz``).
  Termina nativamente en ``...LRGG`` (verificado).
- ``sumo1_reference.pdb``: RCSB 1A5R (SUMO-1 humano, NMR, 10 modelos) --
  se extrajo SOLO el modelo 1 (los modelos NMR son conformeros alternativos
  de la misma cadena, no aportan mas para este proposito) y se recorto la
  cadena hasta el residuo nativo de la diglicina C-terminal madura
  (``...YQEQTGG``), localizado programaticamente buscando el motivo "QTGG"
  en la secuencia en vez de asumir un numero de residuo fijo -- el
  constructo cristalografico/RMN original tiene 4 residuos extra
  (``HSTV``) despues de la diglicina real (cola de maduracion sin procesar
  por SENP), que NO deben tratarse como el C-terminal real. Incluye un tag
  N-terminal de expresion (``GS-``) irrelevante quimicamente (el Mover solo
  manipula el extremo C-terminal).
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

CONJUGATE_PDB_BY_TYPE = {
    "ubiquitination": DATA_DIR / "ubiquitin_reference.pdb",
    "sumoylation": DATA_DIR / "sumo1_reference.pdb",
}

SUPPORTED_PTM_TYPES = frozenset(CONJUGATE_PDB_BY_TYPE)

DEFAULT_N_TAIL_RES = 3
# Emparejados deliberadamente con el default rapido de refine_cycles (3, ver
# init_pyrosetta): permisivos, igual que el propio test de integracion de
# Rosetta usa (scorefilter=50000/SASAfilter=0) cuando corre pocos ciclos --
# CONFIRMADO REAL 2026-08-01 probando el modulo end-to-end: con
# refine_cycles=3 y los "valores buenos" documentados por Rosetta para
# produccion (scorefilter=0/SASAfilter=1000, los que tenia esta constante
# hasta este hallazgo), la conjugacion NUNCA converge lo suficiente en solo
# 3 ciclos y ConjugationFilterError salta siempre -- un footgun real, no
# hipotetico. Si se sube refine_cycles hacia el valor de produccion (20000,
# ver docstring de init_pyrosetta), estos filtros deben apretarse de vuelta
# a scorefilter=0.0/SASAfilter=1000.0 para que sigan siendo un filtro real.
DEFAULT_SCOREFILTER = 50000.0
DEFAULT_SASAFILTER = 0.0

METRIC_NAMES = (
    "lysine_chi3_CG-CD-CE-NZ",
    "lysine_chi4_CD-CE-NZ-C",
    "amide_CE-NZ-C-CA",
    "glycine_psi_NZ-C-CA-N",
    "glycine_phi_C-CA-N-C",
)


class UnsupportedConjugationTypeError(Exception):
    """``ptm_type`` no es ubiquitination ni sumoylation."""


class ResidueMismatchError(Exception):
    """El residuo real en ``position`` no es Lisina."""


class MultiChainPoseError(Exception):
    """``UBQ_GTPaseMover`` exige ``num_chains() == 1`` (runtime_assert real en el .cc)."""


class NonSequentialNumberingError(Exception):
    """La numeracion PDB no coincide 1:1 con la numeracion de pose (ver docstring del modulo)."""


class ConjugationFilterError(Exception):
    """La conjugacion corrio pero no paso los filtros de score/SASA de ``UBQ_GTPaseMover``."""


def _lyx_params_path() -> Path:
    """Resuelve la ruta de ``LYX.params`` relativa a la instalacion real de PyRosetta.

    No hardcodeada a un env/maquina concreta -- ``pyrosetta.__file__`` apunta
    siempre al paquete realmente importado, portable entre instalaciones.
    """
    import pyrosetta

    path = (
        Path(pyrosetta.__file__).parent
        / "database/chemical/residue_type_sets/fa_standard/residue_types/sidechain_conjugation/LYX.params"
    )
    if not path.exists():
        raise FileNotFoundError(
            f"No se encontro LYX.params en '{path}' -- instalacion de PyRosetta inesperada "
            "(verificado 2026-08-01 contra PyRosetta 2026.30, layout database/ junto al paquete)."
        )
    return path


def init_pyrosetta(refine_cycles: int = 3, refine_repack_cycles: int = 3) -> None:
    """Inicializa PyRosetta con las DOS flags reales que ``UBQ_GTPaseMover`` necesita.

    ``-chemical:exclude_patches SidechainConjugation``: sin ella, el sistema
    de patches por defecto bloquea el uso directo de ``LYX`` como
    ``ResidueType`` (misma flag que usa el propio test de integracion de
    Rosetta para esta app: ``-chemical:exclude_patches SidechainConjugation``).

    ``-extra_res_fa <LYX.params>``: confirmado real 2026-08-01 -- sin esta
    flag, ``ResidueTypeSet`` revienta con "The residue LYX could not be
    generated" incluso con la flag anterior puesta (ambas son necesarias,
    ninguna sustituye a la otra).

    ``refine_cycles``/``refine_repack_cycles``: NO son parametros del Mover
    (no hay setter en ``UBQ_GTPaseMover`` -- se leen de
    ``basic::options::option`` en el momento de ``apply()``), por eso se
    fijan aqui, una sola vez por proceso, igual que el resto de flags de
    ``pyrosetta.init()``. El valor por defecto (3/3) es el usado para
    verificar este modulo (16s/sitio en Tau) -- deliberadamente rapido para
    desarrollo/barridos por sitio, NO el valor recomendado por Rosetta para
    una estructura final publicable. Rosetta documenta como "production"
    ``refine_cycles=20000``/``refine_repack_cycles=100`` (comentario real en
    ``UBQ_Gp_LYX-Cterm/options``) -- mismo orden de magnitud de coste que
    ``ddg_estimate.py`` con ``nstruct`` alto, subir esto para un resultado
    final, no para un barrido exploratorio de muchos sitios.
    """
    import pyrosetta

    lyx_params = _lyx_params_path()
    pyrosetta.init(
        f"-mute all -chemical:exclude_patches SidechainConjugation "
        f"-extra_res_fa {lyx_params} -run:min_type dfpmin_armijo "
        f"-AnchoredDesign:refine_cycles {refine_cycles} "
        f"-AnchoredDesign:refine_repack_cycles {refine_repack_cycles}"
    )


def _validate_target(pose, position: int, ptm_type: str) -> None:
    if ptm_type not in CONJUGATE_PDB_BY_TYPE:
        raise UnsupportedConjugationTypeError(
            f"'{ptm_type}' no es un tipo de conjugacion soportado (soportados: "
            f"{sorted(SUPPORTED_PTM_TYPES)})."
        )

    if pose.conformation().num_chains() != 1:
        raise MultiChainPoseError(
            f"UBQ_GTPaseMover exige una pose de una sola cadena (runtime_assert real en "
            f"UBQ_GTPaseMover.cc); esta pose tiene {pose.conformation().num_chains()}."
        )

    pdb_info = pose.pdb_info()
    for i in range(1, pose.size() + 1):
        if pdb_info.number(i) != i:
            raise NonSequentialNumberingError(
                f"La numeracion PDB del residuo de pose {i} es {pdb_info.number(i)}, no {i} -- "
                "UBQ_GTPaseMover.initialize() hace una doble conversion pose->pdb->pose que solo "
                "es correcta si ambas numeraciones coinciden 1:1 (ver docstring del modulo). "
                "Conjugar aqui podria mutar silenciosamente el residuo equivocado."
            )

    actual_residue = pose.residue(position).name1()
    if actual_residue != "K":
        raise ResidueMismatchError(
            f"Posicion {position} es '{actual_residue}', se esperaba Lisina ('K') para "
            f"conjugar '{ptm_type}' -- Rosetta mismo NO valida esto (runtime_assert real "
            "comentada en el .cc), asi que lo hacemos aqui."
        )


def _recompute_metrics(pose, mover) -> dict:
    """Recalcula las 5 metricas por-decoy directamente desde la pose de salida.

    Ver docstring del modulo (tercer hallazgo): el job real donde
    ``analyze_and_filter()`` las escribe es un dummy inaccesible fuera de
    ``JobDistributor``. Misma formula, misma fuente de atom IDs
    (``mover.get_atom_IDs()``, expuesto y verificado real), calculo
    independiente.
    """
    from pyrosetta.rosetta.numeric.conversions import degrees

    atom_ids = mover.get_atom_IDs()
    atom_tree = pose.atom_tree()
    torsions = [
        (atom_ids[1], atom_ids[2], atom_ids[3], atom_ids[4]),
        (atom_ids[2], atom_ids[3], atom_ids[4], atom_ids[5]),
        (atom_ids[3], atom_ids[4], atom_ids[5], atom_ids[6]),
        (atom_ids[4], atom_ids[5], atom_ids[6], atom_ids[7]),
        (atom_ids[5], atom_ids[6], atom_ids[7], atom_ids[8]),
    ]
    return {
        name: degrees(atom_tree.torsion_angle(*atoms))
        for name, atoms in zip(METRIC_NAMES, torsions)
    }


def conjugate(
    pose,
    position: int,
    ptm_type: str,
    n_tail_res: int = DEFAULT_N_TAIL_RES,
    scorefilter: float = DEFAULT_SCOREFILTER,
    sasa_filter: float = DEFAULT_SASAFILTER,
):
    """Conjuga ubiquitina o SUMO en ``position`` (1-based) de ``pose`` via ``UBQ_GTPaseMover``.

    Muta ``pose`` in-place (misma convencion que ``pyrosetta_glycan_patch.attach_glycan``)
    y devuelve ``(pose, metrics)`` -- ``metrics`` son los 5 angulos de la cola
    movil recalculados (ver ``_recompute_metrics``).

    Lanza ``ConjugationFilterError`` si el resultado no pasa los filtros de
    score/SASA configurados (comportamiento real de ``UBQ_GTPaseMover``:
    ``last_move_status`` queda en ``FAIL_RETRY`` en vez de lanzar una
    excepcion nativa -- se traduce aqui a una excepcion real de Python en
    vez de dejar que el caller olvide comprobar el status).
    """
    import pyrosetta
    from pyrosetta.rosetta.protocols.chemically_conjugated_docking import UBQ_GTPaseMover
    from pyrosetta.rosetta.protocols.moves import MoverStatus

    _validate_target(pose, position, ptm_type)

    mover = UBQ_GTPaseMover()
    mover.set_UBQpdb(str(CONJUGATE_PDB_BY_TYPE[ptm_type].resolve()))
    mover.make_index_selector(position)
    mover.set_n_tail_res(n_tail_res)
    mover.set_scorefilter(scorefilter)
    mover.set_SASAfilter(sasa_filter)

    original_cwd = Path.cwd()
    scratch_dir = tempfile.mkdtemp(prefix="ubq_gtpase_scratch_")
    try:
        os.chdir(scratch_dir)
        mover.apply(pose)
    finally:
        os.chdir(original_cwd)
        shutil.rmtree(scratch_dir, ignore_errors=True)

    if mover.get_last_move_status() != MoverStatus.MS_SUCCESS:
        raise ConjugationFilterError(
            f"Conjugacion de '{ptm_type}' en posicion {position} no paso los filtros "
            f"(scorefilter={scorefilter}, SASAfilter={sasa_filter}) -- last_move_status="
            f"{mover.get_last_move_status()}."
        )

    return pose, _recompute_metrics(pose, mover)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Conjuga ubiquitina o SUMO en un residuo real via PyRosetta (Fase A clase 3)."
    )
    parser.add_argument("--pdb-path", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ptm-type", required=True, choices=sorted(SUPPORTED_PTM_TYPES))
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument("--n-tail-res", type=int, default=DEFAULT_N_TAIL_RES)
    parser.add_argument("--scorefilter", type=float, default=DEFAULT_SCOREFILTER)
    parser.add_argument("--sasa-filter", type=float, default=DEFAULT_SASAFILTER)
    parser.add_argument(
        "--refine-cycles", type=int, default=3,
        help="Ciclos de refinamiento Monte Carlo. Default rapido (desarrollo/barrido); "
             "Rosetta recomienda 20000 para una estructura final publicable (ver docstring).",
    )
    parser.add_argument("--refine-repack-cycles", type=int, default=3)
    args = parser.parse_args()

    init_pyrosetta(refine_cycles=args.refine_cycles, refine_repack_cycles=args.refine_repack_cycles)
    from pyrosetta import pose_from_pdb

    pose = pose_from_pdb(args.pdb_path)
    pose, metrics = conjugate(
        pose, args.position, args.ptm_type,
        n_tail_res=args.n_tail_res, scorefilter=args.scorefilter, sasa_filter=args.sasa_filter,
    )
    pose.dump_pdb(args.out_pdb)

    print(f"Conjugacion '{args.ptm_type}' en posicion {args.position} -> '{args.out_pdb}'")
    for name, value in metrics.items():
        print(f"  {name} = {value:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
