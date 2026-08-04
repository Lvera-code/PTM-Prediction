"""Ground truth del panel de validacion biologica (punto 8 del plan de robustez
post-demo-prep, STATUS.md/diario del vault 2026-08-04).

Hasta ahora el pipeline solo se habia corrido de punta a punta contra UNA proteina real
(Tau/MAPT, P10636). Este modulo añade 7 proteinas humanas adicionales, con sitios PTM
reales documentados en literatura -- para poder medir recall real (no solo eyeballing)
sobre datos que el pipeline nunca vio durante el desarrollo/calibracion.

**Regla dura de este proyecto, no negociable** (ver feedback_no_fabricar_datos_cientificos
en la memoria de sesion): cada triple (posicion, tipo, sitio) de aqui viene de una fuente
primaria real y verificable -- una anotacion de UniProt con evidencia ``ECO:0000269``
(experimental, no ``ECO:0000250``/"por similitud") y un PMID confirmado (via NCBI eutils)
de que el titulo del paper corresponde a la modificacion reclamada. Investigado por un
subagente Opus (criterio ya autorizado: verificar una afirmacion cientifica contra fuente
primaria antes de fijarla en produccion) el 2026-08-04 -- ver STATUS.md para el reporte
completo y las 3 trampas de numeracion encontradas.

## Por que AlphaFold, no PDBs experimentales

Verificado empiricamente (``fasta_inputs/validation_panel/*.pdb``, descargados de
``alphafold.ebi.ac.uk/files/AF-<accession>-F1-model_v6.pdb``): la numeracion de residuos
de AlphaFold coincide 1:1 con la numeracion canonica de UniProt para las 7 proteinas de
este panel (confirmado corriendo ``src.utils.structure_parser.parse_structure`` sobre cada
archivo real: longitud de cadena = longitud UniProt exacta en los 7 casos). Los PDBs
experimentales NO tienen esta garantia -- 3 casos reales confirmados via sus registros
DBREF: histonas (PDB = UniProt - 1, el Met inicial se recorta), protrombina (PDB =
UniProt - 43, numeracion de cadena madura tras quitar peptido señal + propeptido), EPO
(PDB = UniProt - 27, numeracion de cadena madura). Por eso el panel usa exclusivamente
AlphaFold v6 (la v4 ya no existe en el servidor -- confirmado 404 durante la
investigacion) y todas las posiciones de este modulo son numeracion UniProt = numeracion
del PDB descargado, offset 0, sin traduccion necesaria.

## Dos niveles de confianza (recomendacion explicita del subagente que investigo esto)

``tier="A"``: sitio confirmado por multiples PMIDs independientes, o por evidencia
estructural directa (ej. registros MODRES de un cristal real), o por una caracterizacion
quimica directa dedicada (no un screening de MS de alto rendimiento). Un fallo del
pipeline aqui es una señal fuerte de problema real.
``tier="B"``: sitio real y con PMID verificado, pero proveniente de un unico screening de
MS a gran escala -- localizacion solida, pero estequiometria real desconocida/posiblemente
baja. Un fallo aqui es evidencia mas debil que en tier A (reportar recall de A y B por
separado, nunca mezclados en un solo numero).

Excluidos deliberadamente (ver el reporte completo de investigacion, STATUS.md): toda
anotacion ``ECO:0000250`` (inferida por similitud, no experimental en humano), sitios con
mecanismo distinto al tipo que aparentan (ej. N-terminal Ser/Ala acetylation, que NO es el
mismo mecanismo que la N6-acetil-lisina que los motores modelan), y el sitio contestado
(HIF-1a K532 acetylation, ver nota en la entrada).
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple

PANEL_DIR = Path(__file__).resolve().parent.parent.parent / "fasta_inputs" / "validation_panel"


@dataclass(frozen=True)
class GroundTruthSite:
    position: int  # 1-based, numeracion UniProt == numeracion del PDB de AlphaFold usado
    residue: str  # letra 1-caracter esperada en esa posicion (verificacion de coherencia)
    ptm_type: str
    tier: str  # "A" o "B" -- ver docstring del modulo
    pmids: Tuple[int, ...]
    is_negative: bool = False  # True = sitio confirmado en literatura como NO modificado
    note: str = ""

    def __post_init__(self):
        if self.tier not in ("A", "B"):
            raise ValueError(f"tier debe ser 'A' o 'B', no {self.tier!r}")
        if not self.pmids:
            raise ValueError(f"GroundTruthSite en posicion {self.position} sin ningun PMID real")


@dataclass(frozen=True)
class PanelEntry:
    name: str
    uniprot_accession: str
    pdb_filename: str
    length: int
    sites: Tuple[GroundTruthSite, ...]

    @property
    def pdb_path(self) -> Path:
        return PANEL_DIR / self.pdb_filename

    @property
    def positives(self) -> Tuple[GroundTruthSite, ...]:
        return tuple(s for s in self.sites if not s.is_negative)

    @property
    def negatives(self) -> Tuple[GroundTruthSite, ...]:
        return tuple(s for s in self.sites if s.is_negative)


# ---------------------------------------------------------------------------
# 1. p53 / TP53 -- P04637 (isoforma canonica, 393 aa)
# ---------------------------------------------------------------------------
_P53 = PanelEntry(
    name="p53", uniprot_accession="P04637", pdb_filename="p53_P04637.pdb", length=393,
    sites=(
        # Tier A: multiples PMIDs independientes, sitios de referencia del campo.
        GroundTruthSite(15, "S", "phosphorylation", "A", (10570149, 11554766, 15866171, 17108107, 17591690, 17967874, 21317932, 28842590)),
        GroundTruthSite(20, "S", "phosphorylation", "A", (10570149, 11447225, 11551930, 12810724, 20041275)),
        GroundTruthSite(46, "S", "phosphorylation", "A", (11740489, 11780126, 16377624, 17349958, 17591690)),
        GroundTruthSite(392, "S", "phosphorylation", "A", (10884347, 11239457, 17108107, 21317932, 22214662, 35618207)),
        GroundTruthSite(382, "K", "acetylation", "A", (10656795, 15448695, 19608861, 20228809, 23431171, 29681526)),
        GroundTruthSite(372, "K", "lys_methylation", "A", (15525938, 16415881), note="SETD7 -- Chuikov et al. Nature 2004"),
        GroundTruthSite(370, "K", "lys_methylation", "A", (17108971, 22864287), note="SMYD2 -- Huang et al. Nature 2006"),
        # Tier B: 1-2 PMIDs, sitios reales pero con menos confirmacion independiente.
        GroundTruthSite(9, "S", "phosphorylation", "B", (18022393,)),
        GroundTruthSite(18, "T", "phosphorylation", "B", (10606744, 10951572, 16704422, 31527692)),
        GroundTruthSite(33, "S", "phosphorylation", "B", (9372954, 17591690)),
        GroundTruthSite(37, "S", "phosphorylation", "B", (17254968,)),
        GroundTruthSite(55, "T", "phosphorylation", "B", (15053879, 20124405)),
        GroundTruthSite(120, "K", "acetylation", "B", (17189187, 19854137, 23431171)),
        GroundTruthSite(305, "K", "acetylation", "B", (12724314,)),
        GroundTruthSite(381, "K", "acetylation", "B", (19608861, 29474172)),
        GroundTruthSite(373, "K", "lys_methylation", "B", (20118233,)),
        GroundTruthSite(333, "R", "arg_methylation", "B", (19011621,)),
        GroundTruthSite(335, "R", "arg_methylation", "B", (19011621,)),
        GroundTruthSite(337, "R", "arg_methylation", "B", (19011621,)),
        GroundTruthSite(24, "K", "ubiquitination", "B", (21597459,)),
        GroundTruthSite(351, "K", "ubiquitination", "B", (19033443,)),
        GroundTruthSite(386, "K", "sumoylation", "B", (11124955, 22214662)),
    ),
)

# ---------------------------------------------------------------------------
# 2. Histona H3.1 -- P68431 (136 aa). Nota: nombres de campo (H3K4 etc) = UniProt - 1
# (el Met inicial se recorta post-traduccionalmente, ver docstring del modulo).
# ---------------------------------------------------------------------------
_HISTONE_H3 = PanelEntry(
    name="histone_h3", uniprot_accession="P68431", pdb_filename="histone_h3_P68431.pdb", length=136,
    sites=(
        # Tier A: marcas clasicas, multiples PMIDs independientes (decadas de evidencia).
        GroundTruthSite(11, "S", "phosphorylation", "A", (10464286, 11856369, 12560483, 15681610, 16185088, 16457588), note="H3S10 (campo)"),
        GroundTruthSite(29, "S", "phosphorylation", "A", (10464286, 11856369, 15681610, 15684425, 16185088, 16457588), note="H3S28 (campo)"),
        GroundTruthSite(5, "K", "lys_methylation", "A", (16267050, 16457588, 17194708), note="H3K4 (campo)"),
        GroundTruthSite(10, "K", "lys_methylation", "A", (11242053, 16185088, 16267050, 16457588, 17194708, 37938770), note="H3K9 (campo)"),
        GroundTruthSite(28, "K", "lys_methylation", "A", (16185088, 16267050, 16627869, 17194708), note="H3K27 (campo)"),
        GroundTruthSite(37, "K", "lys_methylation", "A", (15983376, 16185088, 16267050, 16627869, 17194708), note="H3K36 (campo)"),
        GroundTruthSite(80, "K", "lys_methylation", "A", (15525939, 16267050, 16627869, 17194708), note="H3K79 (campo)"),
        GroundTruthSite(10, "K", "acetylation", "A", (16185088, 16267050, 16457588, 16497732, 16627869, 17194708), note="H3K9 (campo)"),
        GroundTruthSite(15, "K", "acetylation", "A", (16185088, 16267050, 16457588, 16497732, 16627869, 17194708), note="H3K14 (campo)"),
        GroundTruthSite(19, "K", "acetylation", "A", (16267050, 16627869, 17194708, 35939806), note="H3K18 (campo)"),
        GroundTruthSite(24, "K", "acetylation", "A", (16267050, 16457588, 16627869, 17194708), note="H3K23 (campo)"),
        # Tier B: sitios reales con PMID pero de un unico screening MS a gran escala.
        GroundTruthSite(12, "T", "phosphorylation", "B", (12560483, 18066052, 18243098, 22901803), note="H3T11 (campo)"),
        GroundTruthSite(3, "R", "arg_methylation", "B", (17898714, 18077460, 18079182), note="H3R2 (campo)"),
        GroundTruthSite(3, "R", "citrullination", "B", (16567635,), note="H3R2 (campo)"),
        GroundTruthSite(9, "R", "citrullination", "B", (15345777, 16567635), note="H3R8 (campo)"),
        GroundTruthSite(5, "K", "crotonylation", "B", (21925322, 28497810), note="H3K4 (campo)"),
        GroundTruthSite(10, "K", "crotonylation", "B", (21925322, 28497810), note="H3K9 (campo)"),
        GroundTruthSite(15, "K", "succinylation", "B", (22389435,), note="H3K14 (campo)"),
        GroundTruthSite(123, "K", "succinylation", "B", (22389435, 27436229), note="H3K122 (campo)"),
        GroundTruthSite(15, "K", "glutarylation", "B", (31542297,), note="H3K14 (campo)"),
        GroundTruthSite(19, "K", "ubiquitination", "B", (27595565, 29053958), note="H3K18 (campo)"),
    ),
)

# ---------------------------------------------------------------------------
# 3. Histona H4 -- P62805 (103 aa). Mismo caveat de nombres de campo que H3.
# ---------------------------------------------------------------------------
_HISTONE_H4 = PanelEntry(
    name="histone_h4", uniprot_accession="P62805", pdb_filename="histone_h4_P62805.pdb", length=103,
    sites=(
        # Tier A.
        GroundTruthSite(48, "S", "phosphorylation", "A", (17081983, 18669648, 19690332, 20068231, 21406692, 21724829, 23186163, 24275569), note="H4S47 (campo)"),
        GroundTruthSite(6, "K", "acetylation", "A", (2474456, 7664735, 19608861), note="H4K5 (campo); PMID 2474456 = anticuerpo especifico de sitio"),
        GroundTruthSite(9, "K", "acetylation", "A", (2474456, 7664735, 19608861), note="H4K8 (campo)"),
        GroundTruthSite(13, "K", "acetylation", "A", (2474456, 7664735, 17967882, 19608861), note="H4K12 (campo)"),
        GroundTruthSite(17, "K", "acetylation", "A", (2474456, 7664735, 17967882, 19608861), note="H4K16 (campo)"),
        GroundTruthSite(21, "K", "lys_methylation", "A", (12086618, 15964846, 17967882, 20622854, 27338793, 32209475), note="H4K20 (campo); PR-Set7, Fang et al. Mol Cell 2002"),
        # Tier B.
        GroundTruthSite(52, "Y", "phosphorylation", "B", (15592455, 20068231), note="H4Y51 (campo)"),
        GroundTruthSite(4, "R", "arg_methylation", "B", (11387442, 11448779, 15345777), note="H4R3 (campo); PRMT1"),
        GroundTruthSite(4, "R", "citrullination", "B", (15345777, 16567635), note="H4R3 (campo)"),
        GroundTruthSite(6, "K", "crotonylation", "B", (21925322, 28497810), note="H4K5 (campo)"),
        GroundTruthSite(13, "K", "succinylation", "B", (22389435,), note="H4K12 (campo)"),
        GroundTruthSite(6, "K", "glutarylation", "B", (31542297,), note="H4K5 (campo)"),
        GroundTruthSite(13, "K", "sumoylation", "B", (25218447,), note="H4K12 (campo)"),
        GroundTruthSite(92, "K", "ubiquitination", "B", (19818714,), note="H4K91 (campo)"),
        GroundTruthSite(92, "K", "acetylation", "B", (19818714,), note="H4K91 (campo)"),
    ),
)

# ---------------------------------------------------------------------------
# 4. Protrombina (F2) -- P00734 (622 aa). El sitio mas fuerte del panel: 10/10
# Gla confirmados independientemente por registros MODRES de un cristal real (6C2W),
# ademas de la secuenciacion original.
# ---------------------------------------------------------------------------
_PROTHROMBIN = PanelEntry(
    name="prothrombin", uniprot_accession="P00734", pdb_filename="prothrombin_P00734.pdb", length=622,
    sites=(
        GroundTruthSite(49, "E", "gamma_carboxyglutamic_acid", "A", (3759958, 6305407), note="Gla6 (campo); confirmado en MODRES de 6C2W"),
        GroundTruthSite(50, "E", "gamma_carboxyglutamic_acid", "A", (3759958, 6305407), note="Gla7 (campo); MODRES 6C2W"),
        GroundTruthSite(57, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla14 (campo); MODRES 6C2W"),
        GroundTruthSite(59, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla16 (campo); MODRES 6C2W"),
        GroundTruthSite(62, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla19 (campo); MODRES 6C2W"),
        GroundTruthSite(63, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla20 (campo); MODRES 6C2W"),
        GroundTruthSite(68, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla25 (campo); MODRES 6C2W"),
        GroundTruthSite(69, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla26 (campo); MODRES 6C2W"),
        GroundTruthSite(72, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla29 (campo); MODRES 6C2W"),
        GroundTruthSite(75, "E", "gamma_carboxyglutamic_acid", "A", (6305407,), note="Gla32 (campo); MODRES 6C2W"),
        GroundTruthSite(121, "N", "n_linked_glycosylation", "A", (14760718, 16335952, 22171320), note="sequon N-I-T confirmado en secuencia real"),
        GroundTruthSite(143, "N", "n_linked_glycosylation", "A", (14760718, 16335952, 22171320), note="sequon N-S-T confirmado en secuencia real"),
        GroundTruthSite(416, "N", "n_linked_glycosylation", "A", (16335952, 19139490, 19159218, 19838169, 873923), note="sequon N-F-T confirmado en secuencia real"),
    ),
)

# ---------------------------------------------------------------------------
# 5. HIF-1a -- Q16665 (826 aa). Unica fuente limpia de hydroxylation real.
# pLDDT bajo (proteina intrinsecamente desordenada fuera del dominio bHLH-PAS) --
# esperar que el motor estructural tenga menos confianza aqui, no es fallo de datos.
# ---------------------------------------------------------------------------
_HIF1A = PanelEntry(
    name="hif1a", uniprot_accession="Q16665", pdb_filename="hif1a_Q16665.pdb", length=826,
    sites=(
        GroundTruthSite(564, "P", "hydroxylation", "A", (11292861, 11566883, 12351678, 25974097), note="el sitio de hidroxilacion mejor confirmado de toda la biologia humana (PHD2/VHL)"),
        GroundTruthSite(402, "P", "hydroxylation", "A", (11566883,), note="paper landmark (EMBO J 2001), aceptado universalmente pese a 1 PMID"),
        GroundTruthSite(803, "N", "hydroxylation", "A", (12080085,), note="FIH-1, Genes Dev 2002, paper landmark"),
        GroundTruthSite(391, "K", "sumoylation", "A", (15465032, 17610843)),
        GroundTruthSite(477, "K", "sumoylation", "A", (15465032, 17610843)),
        GroundTruthSite(532, "K", "ubiquitination", "A", (16862177,)),
        GroundTruthSite(538, "K", "ubiquitination", "A", (16862177,)),
        GroundTruthSite(547, "K", "ubiquitination", "A", (16862177,)),
        # Tier B.
        GroundTruthSite(800, "C", "s_nitrosylation", "B", (12914934,)),
        GroundTruthSite(247, "S", "phosphorylation", "B", (20699359,)),
        GroundTruthSite(551, "S", "phosphorylation", "B", (20889502,)),
        GroundTruthSite(709, "K", "acetylation", "B", (24681946,)),
        # K532 acetylation (mismo residuo que el sitio de ubiquitinacion de arriba,
        # tipo distinto) DELIBERADAMENTE EXCLUIDA: contestada en la literatura, ver
        # docstring del modulo y STATUS.md.
    ),
)

# ---------------------------------------------------------------------------
# 6. Eritropoyetina (EPO) -- P01588 (193 aa, precursor con peptido señal incluido).
# Los 4 sitios vienen de una unica caracterizacion quimica DIRECTA (Lai et al. 1986,
# mapeo de peptidos de la glicoproteina purificada), no de una prediccion ni un
# screening -- por eso se tratan como Tier A pese a ser 1 solo PMID cada uno.
# ---------------------------------------------------------------------------
_EPO = PanelEntry(
    name="epo", uniprot_accession="P01588", pdb_filename="epo_P01588.pdb", length=193,
    sites=(
        GroundTruthSite(51, "N", "n_linked_glycosylation", "A", (3949763,), note="EPO N24 (numeracion madura); sequon N-I-T"),
        GroundTruthSite(65, "N", "n_linked_glycosylation", "A", (3949763,), note="EPO N38 (numeracion madura); sequon N-I-T"),
        GroundTruthSite(110, "N", "n_linked_glycosylation", "A", (3949763,), note="EPO N83 (numeracion madura); sequon N-S-S"),
        GroundTruthSite(153, "S", "o_linked_glycosylation", "A", (3949763,), note="EPO S126 (numeracion madura); O-GalNAc"),
    ),
)

# ---------------------------------------------------------------------------
# 7. KIT ligand / SCF -- P21583 (273 aa). Control NEGATIVO principal: Asn-97 es un
# sequon N-X-S perfectamente valido que la propia UniProt (citando PMID 1381905)
# reporta EXPLICITAMENTE como NO glicosilado en ninguna de las 2 isoformas conocidas
# -- distingue un motor con especificidad real de uno que solo empareja el motivo.
# ---------------------------------------------------------------------------
_KIT_LIGAND_SCF = PanelEntry(
    name="kit_ligand_scf", uniprot_accession="P21583", pdb_filename="kit_ligand_scf_P21583.pdb", length=273,
    sites=(
        GroundTruthSite(90, "N", "n_linked_glycosylation", "A", (1381905,), note="parcialmente glicosilado en LMW-SCF; sequon N-I-S"),
        GroundTruthSite(145, "N", "n_linked_glycosylation", "A", (1381905,), note="totalmente glicosilado; sequon N-R-S"),
        GroundTruthSite(167, "S", "o_linked_glycosylation", "A", (1381905,)),
        GroundTruthSite(168, "T", "o_linked_glycosylation", "A", (1381905,)),
        GroundTruthSite(180, "T", "o_linked_glycosylation", "A", (1381905,)),
        GroundTruthSite(
            97, "N", "n_linked_glycosylation", "A", (1381905,), is_negative=True,
            note="CONTROL NEGATIVO: sequon N-Y-S valido, pero UniProt/PMID 1381905 confirma "
                 "explicitamente que NO esta glicosilado en ninguna de las 2 isoformas (LMW/HMW-SCF)",
        ),
    ),
)

PANEL: Tuple[PanelEntry, ...] = (
    _P53, _HISTONE_H3, _HISTONE_H4, _PROTHROMBIN, _HIF1A, _EPO, _KIT_LIGAND_SCF,
)
