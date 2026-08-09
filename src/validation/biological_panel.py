"""Ground truth del panel de validacion biologica.

**Origen y evolucion:**
- 2026-08-04: panel inicial de 7 proteinas curadas a mano (punto 8 del plan de
  robustez post-demo-prep, ver STATUS.md), reemplazando la unica proteina real
  usada hasta entonces (Tau/MAPT, P10636).
- 2026-08-09: migrado a una seleccion SISTEMATICA de proteinas derivada de
  dbPTM (decision en el vault:
  01-Proyectos/PTM-Prediction/Decisiones/2026-08-09-panel-validacion-dbptm-reemplaza-curacion-manual.md).
  ``scripts/select_recall_subset.py`` elige ~25 proteinas humanas por greedy
  set-cover sobre los 17 tipos canonicos + evidencia tier A en
  ``data/dbptm/lookup.sqlite3`` (base grande derivada de dbPTM, ver
  ``scripts/import_dbptm_panel.py``) -- reemplaza la curacion manual de QUE
  proteinas entran, pero NO el rigor de verificacion de cada sitio (ver regla
  dura abajo): dbPTM solo aporta candidatas, cada una se verifico PMID por PMID
  contra NCBI eutils antes de entrar aqui, exactamente igual que las 7
  originales (4 subagentes Opus en paralelo, un catalogo real de basura de
  dbPTM encontrada y descartada: identificadores que no son PMIDs, PMIDs reales
  pero de campos completamente ajenos, tipo de PTM equivocado -- el agregador
  "arrastra" listas de un tipo de sitio a otro del mismo residuo -- y papers de
  predictores computacionales citados como si fueran evidencia experimental).
  p53 e histona H3 sobreviven de la seleccion original, fusionadas con sitios
  nuevos verificados el 2026-08-09 (sitios viejos intactos, nunca re-verificados
  sin motivo). ``kit_ligand_scf`` (unico control negativo real del panel) se
  preserva intacta: dbPTM no tiene anotaciones de "confirmado como NO
  modificado", asi que no puede salir de ahi.

**Regla dura de este proyecto, no negociable** (ver feedback_no_fabricar_datos_cientificos
en la memoria de sesion): cada triple (posicion, tipo, sitio) de aqui viene de una fuente
primaria real y verificable -- una anotacion de UniProt con evidencia ``ECO:0000269``
(experimental, no ``ECO:0000250``/"por similitud") y un PMID confirmado (via NCBI eutils)
de que el titulo del paper corresponde a la modificacion reclamada.

## Por que AlphaFold, no PDBs experimentales

Verificado empiricamente (``inputs/*.pdb``, descargados de la API de AlphaFold DB): la
numeracion de residuos de AlphaFold coincide 1:1 con la numeracion canonica de UniProt
para las proteinas de este panel (confirmado corriendo
``src.utils.structure_parser.parse_structure`` sobre cada archivo real: longitud de
cadena = longitud UniProt exacta). Los PDBs experimentales NO tienen esta garantia -- 3
casos reales confirmados via sus registros DBREF en la investigacion original: histonas
(PDB = UniProt - 1, el Met inicial se recorta), protrombina (PDB = UniProt - 43,
numeracion de cadena madura tras quitar peptido señal + propeptido), EPO (PDB = UniProt
- 27, numeracion de cadena madura) -- documentados aqui como advertencia para quien
agregue un PDB experimental en el futuro. Por eso el panel usa exclusivamente AlphaFold
v6 (la v4 ya no existe en el servidor) y todas las posiciones de este modulo son
numeracion UniProt = numeracion del PDB descargado, offset 0, sin traduccion necesaria.

## Dos niveles de confianza

``tier="A"``: sitio confirmado por multiples PMIDs independientes, o por evidencia
estructural directa (ej. registros MODRES de un cristal real), o por una caracterizacion
quimica directa dedicada (no un screening de MS de alto rendimiento). Un fallo del
pipeline aqui es una señal fuerte de problema real.
``tier="B"``: sitio real y con PMID verificado, pero proveniente de un unico screening de
MS a gran escala -- localizacion solida, pero estequiometria real desconocida/posiblemente
baja. Un fallo aqui es evidencia mas debil que en tier A (reportar recall de A y B por
separado, nunca mezclados en un solo numero).

Excluidos deliberadamente: toda anotacion ``ECO:0000250`` (inferida por similitud, no
experimental en humano); sitios con mecanismo distinto al tipo que aparentan (ej.
N-terminal Ser/Ala/Met acetylation, que NO es el mismo mecanismo que la N6-acetil-lisina
que los motores modelan -- encontrado repetidamente en el borrador de dbPTM); y PMIDs
que dbPTM propaga pero no sustentan la afirmacion exacta de ese sitio.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple

PANEL_DIR = Path(__file__).resolve().parent.parent.parent / "inputs"


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
# 1. p53 / TP53 -- P04637 (393 aa). ENTRADA FUSIONADA 2026-08-09: 22 sitios
# verificados en 2026-08-04 (intactos) + 12 nuevos verificados el 2026-08-09
# desde el borrador dbPTM. Descartes notables del borrador: R110/R209
# arg_methylation citaban papers de 1971 sobre permeabilidad de liposomas y
# pronasa inmovilizada (basura del agregador); los PMIDs que añadia a
# R333/R335/R337 eran articulos de fisica de particulas (Phys Rev D).
# ---------------------------------------------------------------------------
_P53 = PanelEntry(
    name="p53", uniprot_accession="P04637", pdb_filename="p53_P04637.pdb", length=393,
    sites=(
        # --- Verificados 2026-08-04 ---
        GroundTruthSite(15, "S", "phosphorylation", "A", (10570149, 11554766, 15866171, 17108107, 17591690, 17967874, 21317932, 28842590)),
        GroundTruthSite(20, "S", "phosphorylation", "A", (10570149, 11447225, 11551930, 12810724, 20041275)),
        GroundTruthSite(46, "S", "phosphorylation", "A", (11740489, 11780126, 16377624, 17349958, 17591690)),
        GroundTruthSite(392, "S", "phosphorylation", "A", (10884347, 11239457, 17108107, 21317932, 22214662, 35618207)),
        GroundTruthSite(382, "K", "acetylation", "A", (10656795, 15448695, 19608861, 20228809, 23431171, 29681526)),
        GroundTruthSite(372, "K", "lys_methylation", "A", (15525938, 16415881), note="SETD7 -- Chuikov et al. Nature 2004"),
        GroundTruthSite(370, "K", "lys_methylation", "A", (17108971, 22864287), note="SMYD2 -- Huang et al. Nature 2006"),
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
        # --- Nuevos, verificados 2026-08-09 (borrador dbPTM) ---
        GroundTruthSite(315, "S", "phosphorylation", "A", (10884347, 14702041), note="AURKA/CDK1/CDK2; UniProt ECO:0000269"),
        GroundTruthSite(320, "K", "acetylation", "A", (9891054, 22696202, 38043200),
                        note="sitio de PCAF; residuo explicito en los titulos de 22696202 ('lysine 317/320') "
                             "y 38043200 ('deacetylating p53 at lysine 320'). Cluster K319-K320-K321; "
                             "UniProt anota acetilo en 321 solo ECO:0000250, la literatura humana usa K320"),
        GroundTruthSite(370, "K", "acetylation", "A", (19155208, 21057544),
                        note="19155208 mapea por LC-MS/MS acetilacion en K305/370/372/373/381/382/386; "
                             "21057544 es dedicado ('p85alpha mediates p53 K370 acetylation by p300')"),
        GroundTruthSite(372, "K", "acetylation", "B", (19155208,), note="unico mapeo MS dedicado a p53 que nombra K372"),
        GroundTruthSite(382, "K", "lys_methylation", "A", (17707234, 20870725, 22864287),
                        note="SET8/KMT5A; residuo explicito en 2 de los 3 titulos"),
        GroundTruthSite(149, "S", "o_linked_glycosylation", "A", (16964247,),
                        note="O-GlcNAc; 'Ser 149 of p53 is O-GlcNAcylated' explicito en el abstract"),
        GroundTruthSite(124, "C", "glutathionylation", "A", (17555331,),
                        note="Velu et al. Biochemistry 2007; el abstract nombra cisteinas 124, 141 y 182 (MS directa)"),
        GroundTruthSite(141, "C", "glutathionylation", "A", (17555331,), note="idem; C141 = la mas reactiva segun el paper"),
        GroundTruthSite(182, "C", "glutathionylation", "A", (17555331,), note="idem"),
        GroundTruthSite(292, "K", "ubiquitination", "B", (19536131,), note="MKRN1; cross-link ECO:0000269 en UniProt"),
        GroundTruthSite(370, "K", "ubiquitination", "B", (11046142,),
                        note="Rodriguez et al. 2000: mutante 6KR (K370/372/373/381/382/386) resistente a "
                             "ubiquitinacion por Mdm2 -- evidencia mutacional de grupo, no mapeo por residuo"),
        GroundTruthSite(382, "K", "ubiquitination", "B", (11046142,), note="idem 6KR"),
    ),
)

# ---------------------------------------------------------------------------
# 2. Histona H3.1 -- P68431 (136 aa). Nota: nombres de campo (H3K4 etc) = UniProt - 1.
# ENTRADA FUSIONADA 2026-08-09: 21 sitios verificados en 2026-08-04 (intactos) +
# 21 nuevos verificados el 2026-08-09. Descarte notable: las 5 "sumoylation" del
# borrador eran cross-links de UBIQUITINA en UniProt, no de SUMO -- listas de
# PMIDs copiadas literalmente de acetilacion/metilacion del mismo residuo.
# ---------------------------------------------------------------------------
_HISTONE_H3 = PanelEntry(
    name="histone_h3", uniprot_accession="P68431", pdb_filename="histone_h3_P68431.pdb", length=136,
    sites=(
        # --- Verificados 2026-08-04 ---
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
        # --- Nuevos, verificados 2026-08-09 (borrador dbPTM) ---
        GroundTruthSite(37, "K", "acetylation", "A", (17189264, 17194708),
                        note="H3K36 (campo); 17189264 = 'Identification of histone H3 lysine 36 acetylation...'"),
        GroundTruthSite(80, "K", "acetylation", "B", (17194708,), note="H3K79 (campo); un unico estudio de MS de histonas"),
        GroundTruthSite(123, "K", "acetylation", "A", (19520870, 23415232),
                        note="H3K122 (campo); ambos titulos nombran el sitio del dyad/H3K122"),
        GroundTruthSite(18, "R", "arg_methylation", "A", (15345777, 15471871, 16497732),
                        note="H3R17 (campo); CARM1, residuo explicito en el titulo de 16497732"),
        GroundTruthSite(18, "R", "citrullination", "B", (15345777, 16567635), note="H3R17 (campo); PAD4"),
        GroundTruthSite(27, "R", "citrullination", "B", (16567635,), note="H3R26 (campo); PAD4"),
        GroundTruthSite(19, "K", "crotonylation", "B", (21925322,), note="H3K18 (campo)"),
        GroundTruthSite(24, "K", "crotonylation", "B", (21925322, 28497810), note="H3K23 (campo)"),
        GroundTruthSite(28, "K", "crotonylation", "B", (21925322,), note="H3K27 (campo)"),
        GroundTruthSite(19, "K", "glutarylation", "B", (31542297,), note="H3K18 (campo)"),
        GroundTruthSite(24, "K", "glutarylation", "B", (31542297,), note="H3K23 (campo)"),
        GroundTruthSite(28, "K", "glutarylation", "B", (31542297,), note="H3K27 (campo)"),
        GroundTruthSite(57, "K", "glutarylation", "B", (31542297,), note="H3K56 (campo)"),
        GroundTruthSite(123, "K", "lys_methylation", "B", (16267050, 17194708), note="H3K122 (campo)"),
        GroundTruthSite(11, "S", "o_linked_glycosylation", "B", (21896475,),
                        note="H3S10 (campo); residuo explicito en el titulo (O-GlcNAc)"),
        GroundTruthSite(33, "T", "o_linked_glycosylation", "B", (22371497,),
                        note="H3T32 (campo); 'Mass spectrometry analysis identified threonine 32' en el abstract"),
        GroundTruthSite(58, "S", "phosphorylation", "B", (20850016,), note="H3S57 (campo)"),
        GroundTruthSite(81, "T", "phosphorylation", "B", (20850016,), note="H3T80 (campo)"),
        GroundTruthSite(57, "K", "succinylation", "B", (22389435,), note="H3K56 (campo)"),
        GroundTruthSite(80, "K", "succinylation", "A", (22389435, 29211711),
                        note="H3K79 (campo); 29211711 = KAT2A como succiniltransferasa de H3 (Nature 2017)"),
        GroundTruthSite(24, "K", "ubiquitination", "B", (27595565, 29053958), note="H3K23 (campo); cross-link con ubiquitina (UHRF1/DNMT1)"),
    ),
)

# ---------------------------------------------------------------------------
# 3. Alfa-enolasa (ENO1) -- P06733 (434 aa). Verificado 2026-08-09. Del borrador
# dbPTM: 30 candidatos -> 6 aceptados. Peor señal/ruido del grupo 1: 5
# "acetilaciones tier A" eran ECO:0000250; 9 sitios de metilacion/sumoilacion/
# ubiquitinacion citaban el acetiloma de Choudhary (arrastre de tipo).
# ---------------------------------------------------------------------------
_ENOA = PanelEntry(
    name="enoa", uniprot_accession="P06733", pdb_filename="enoa_P06733.pdb", length=434,
    sites=(
        GroundTruthSite(233, "K", "malonylation", "A", (21908771, 26320211),
                        note="UniProt ECO:0000269 (21908771, primera identificacion de sustratos de malonilacion)"),
        GroundTruthSite(420, "K", "malonylation", "A", (21908771, 26320211), note="idem"),
        GroundTruthSite(202, "K", "sumoylation", "B", (28112733,),
                        note="cross-link SUMO2 en UniProt (ECO:0007744); mapeo sitio-especifico del SUMO-proteoma humano"),
        GroundTruthSite(44, "Y", "phosphorylation", "B", (15592455,), note="Rush et al. 2005, perfilado de fosfotirosina por inmunoafinidad"),
        GroundTruthSite(263, "S", "phosphorylation", "B", (18669648, 23186163, 24275569)),
        GroundTruthSite(272, "S", "phosphorylation", "B", (19690332, 23186163)),
    ),
)

# ---------------------------------------------------------------------------
# 4. Quininogeno-1 (KNG1) -- P01042 (644 aa). Verificado 2026-08-09. Del borrador:
# 17 candidatos -> 8 aceptados. P383 = Hyp3 de la bradiquinina (bradiquinina =
# residuos 381-389 de la cadena precursora).
# ---------------------------------------------------------------------------
_KNG1 = PanelEntry(
    name="kng1", uniprot_accession="P01042", pdb_filename="kng1_P01042.pdb", length=644,
    sites=(
        GroundTruthSite(383, "P", "hydroxylation", "A", (3182782, 3366244),
                        note="4-hidroxiprolina parcial; [Hyp3]-bradiquinina aislada de fluido ascitico y de orina humana"),
        GroundTruthSite(294, "N", "n_linked_glycosylation", "A", (12754519, 14760718, 16335952, 19139490, 19159218, 19838169),
                        note="6 glicoproteomicas independientes, todas con ECO:0000269 en UniProt"),
        GroundTruthSite(169, "N", "n_linked_glycosylation", "A", (3484703, 14760718, 16335952, 19159218),
                        note="3484703 = secuenciacion completa de la cadena pesada (Eur J Biochem 1986)"),
        GroundTruthSite(205, "N", "n_linked_glycosylation", "A", (3484703, 16335952, 19139490, 19159218)),
        GroundTruthSite(48, "N", "n_linked_glycosylation", "B", (16335952, 19139490), note="solo 2 screenings de MS"),
        GroundTruthSite(533, "T", "o_linked_glycosylation", "A", (4054110,),
                        note="secuenciacion de aminoacidos de la cadena ligera (caracterizacion quimica directa)"),
        GroundTruthSite(546, "T", "o_linked_glycosylation", "A", (4054110,), note="idem"),
        GroundTruthSite(332, "S", "phosphorylation", "A", (19824718, 24275569, 26091039), note="FAM20C; ECO:0000269 + 2 fosfoproteomas"),
    ),
)

# ---------------------------------------------------------------------------
# 5. Factor VII de coagulacion (F7) -- P08709 (466 aa). Verificado 2026-08-09.
# UniProt documenta 10 residuos Gla con los mismos 2 PMIDs (3264725, 3486420);
# solo se incluyen los 5 que aparecian en el borrador -- añadir E79/E85/E86/E89/E95
# seria trivial y correcto pero queda fuera del alcance de esta verificacion.
# D123 hydroxylation es tier A: caracterizacion quimica directa del FVIIa purificado.
# ---------------------------------------------------------------------------
_FA7 = PanelEntry(
    name="fa7", uniprot_accession="P08709", pdb_filename="fa7_P08709.pdb", length=466,
    sites=(
        GroundTruthSite(66, "E", "gamma_carboxyglutamic_acid", "A", (3264725, 3486420), note="Gla6 (numeracion madura)"),
        GroundTruthSite(67, "E", "gamma_carboxyglutamic_acid", "A", (3264725, 3486420), note="Gla7"),
        GroundTruthSite(74, "E", "gamma_carboxyglutamic_acid", "A", (3264725, 3486420), note="Gla14"),
        GroundTruthSite(76, "E", "gamma_carboxyglutamic_acid", "A", (3264725, 3486420), note="Gla16"),
        GroundTruthSite(80, "E", "gamma_carboxyglutamic_acid", "A", (3264725, 3486420), note="Gla20"),
        GroundTruthSite(123, "D", "hydroxylation", "A", (3264725,),
                        note="(3R)-3-hidroxiaspartato en el dominio EGF-1; caracterizacion quimica directa del FVIIa de plasma"),
        GroundTruthSite(205, "N", "n_linked_glycosylation", "A", (3264725, 19167329)),
        GroundTruthSite(382, "N", "n_linked_glycosylation", "A", (3264725, 19167329)),
        GroundTruthSite(112, "S", "o_linked_glycosylation", "A", (1904059, 2129367, 2511201, 21949356),
                        note="Ser52 (madura); O-Glc/O-Xyl del dominio EGF-1, 4 trabajos independientes"),
        GroundTruthSite(120, "S", "o_linked_glycosylation", "A", (1904059, 9023546), note="Ser60 (madura); O-fucosa"),
    ),
)

# ---------------------------------------------------------------------------
# 6. MAP4 (microtubule-associated protein 4) -- P27816 (1152 aa). Verificado
# 2026-08-09. Del borrador: 34 sitios -> 31 aceptados. S1073 es tier A por
# caracterizacion dedicada: Illenberger 1996 (JBC) mide in vitro la fosforilacion
# de los motivos KXGS de MAP4 (S1073 cae en uno) por MARK/p110mark.
# ---------------------------------------------------------------------------
_MAP4 = PanelEntry(
    name="map4", uniprot_accession="P27816", pdb_filename="map4_P27816.pdb", length=1152,
    sites=(
        GroundTruthSite(56, "K", "acetylation", "B", (23236377, 26051181)),
        GroundTruthSite(346, "K", "acetylation", "A", (23236377, 23954790, 26051181)),
        GroundTruthSite(368, "K", "acetylation", "A", (22424773, 23236377, 23749302, 23954790, 25953088, 26051181)),
        GroundTruthSite(946, "K", "acetylation", "A", (25953088, 26051181, 27452117)),
        GroundTruthSite(535, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(635, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(1011, "K", "lys_methylation", "B", (24129315,)),
        GroundTruthSite(346, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(352, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(464, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(769, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(832, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(329, "S", "o_linked_glycosylation", "A", (16408927, 20068230, 20305658, 21158410, 23301498, 27655845)),
        GroundTruthSite(349, "T", "o_linked_glycosylation", "A", (16408927, 20068230, 20305658, 21158410, 23301498, 27655845)),
        GroundTruthSite(580, "S", "o_linked_glycosylation", "A", (16408927, 20068230, 20305658, 21158410, 23301498, 27655845)),
        GroundTruthSite(707, "T", "o_linked_glycosylation", "A", (16408927, 20068230, 20305658, 21158410, 23301498, 27655845)),
        GroundTruthSite(742, "S", "o_linked_glycosylation", "A", (16408927, 20068230, 20305658, 21158410, 23301498, 27655845)),
        GroundTruthSite(280, "S", "phosphorylation", "A", (15302935, 17081983, 17287340, 18669648, 19060867, 20068231)),
        GroundTruthSite(507, "S", "phosphorylation", "A", (15302935, 17081983, 17287340, 18669648, 19060867, 20068231)),
        GroundTruthSite(521, "T", "phosphorylation", "A", (15302935, 16964243, 17081983, 17287340, 18669648, 19060867)),
        GroundTruthSite(636, "S", "phosphorylation", "A", (17081983, 17287340, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(1073, "S", "phosphorylation", "A", (8631898, 17081983, 18669648, 19060867, 20068231, 21406692, 21712546),
                        note="KVGS del dominio de union a microtubulos; Illenberger 1996 (JBC) caracteriza la fosforilacion de los motivos KXGS de MAP4 por MARK/p110mark"),
        GroundTruthSite(635, "C", "s_nitrosylation", "B", (16648260,)),
        GroundTruthSite(352, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(498, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(269, "K", "sumoylation", "B", (25114211,)),
        GroundTruthSite(838, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(312, "K", "ubiquitination", "A", (21890473, 21906983, 23000965, 29967540, 33845483)),
        GroundTruthSite(332, "K", "ubiquitination", "A", (21890473, 21906983, 29967540, 33845483)),
        GroundTruthSite(346, "K", "ubiquitination", "A", (21890473, 21906983, 29967540, 32142685, 33845483)),
        GroundTruthSite(498, "K", "ubiquitination", "A", (21890473, 21906983, 23000965, 29967540, 33845483)),
    ),
)

# ---------------------------------------------------------------------------
# 7. Nucleoporina NUP153 -- P49790 (1475 aa). Verificado 2026-08-09. Del
# borrador: 26 sitios -> 20 aceptados. Los 5 sitios O-GlcNAc son consistentes
# con la biologia conocida: NUP153 es una de las nucleoporinas O-GlcNAciladas
# clasicas, cada sitio con estudios de mapeo dedicados (no revisiones).
# ---------------------------------------------------------------------------
_NU153 = PanelEntry(
    name="nu153", uniprot_accession="P49790", pdb_filename="nu153_P49790.pdb", length=1475,
    sites=(
        GroundTruthSite(194, "K", "acetylation", "B", (23749302, 25953088)),
        GroundTruthSite(448, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(460, "K", "acetylation", "A", (23236377, 23954790, 25953088, 26051181)),
        GroundTruthSite(1120, "K", "acetylation", "A", (19608861, 25953088, 26051181)),
        GroundTruthSite(153, "S", "o_linked_glycosylation", "A", (20068230, 20305658, 21300897, 21740066, 23301498, 27114449)),
        GroundTruthSite(164, "S", "o_linked_glycosylation", "A", (20068230, 20305658, 21300897, 21740066, 23301498, 27114449)),
        GroundTruthSite(199, "S", "o_linked_glycosylation", "A", (20068230, 20305658, 21300897, 21740066, 23301498, 27114449)),
        GroundTruthSite(1179, "T", "o_linked_glycosylation", "A", (20068230, 20305658, 21300897, 21740066, 23301498, 27114449)),
        GroundTruthSite(1180, "T", "o_linked_glycosylation", "A", (20068230, 20305658, 21300897, 21740066, 23301498, 27114449)),
        GroundTruthSite(209, "S", "phosphorylation", "A", (16964243, 18669648, 19060867, 20068231, 21406692, 21712546)),
        GroundTruthSite(334, "S", "phosphorylation", "A", (16964243, 17081983, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(338, "S", "phosphorylation", "A", (16964243, 17081983, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(516, "S", "phosphorylation", "A", (17081983, 18669648, 19060867, 20068231, 21406692, 21712546)),
        GroundTruthSite(522, "S", "phosphorylation", "A", (18669648, 19060867, 19691289, 20068231, 21406692, 21712546)),
        GroundTruthSite(353, "K", "sumoylation", "A", (25114211, 25218447, 28112733)),
        GroundTruthSite(200, "K", "ubiquitination", "A", (21890473, 21906983, 23503661, 27667366, 29967540, 33845483)),
        GroundTruthSite(263, "K", "ubiquitination", "A", (21890473, 23000965, 23503661, 29967540, 33845483)),
        GroundTruthSite(968, "K", "ubiquitination", "A", (21890473, 21906983, 23503661, 29967540)),
        GroundTruthSite(1010, "K", "ubiquitination", "A", (21890473, 21906983, 23503661, 29967540)),
        GroundTruthSite(1157, "K", "ubiquitination", "A", (21890473, 21906983, 23503661, 27667366, 29967540)),
    ),
)

# ---------------------------------------------------------------------------
# 8. Treacle / TCOF1 -- Q13428 (1488 aa). Verificado 2026-08-09. Del borrador:
# 26 sitios -> 23 aceptados. Malonilaciones perdieron el PMID 33225896 (Mal-Prec,
# predictor computacional, no evidencia experimental).
# ---------------------------------------------------------------------------
_TCOF = PanelEntry(
    name="tcof", uniprot_accession="Q13428", pdb_filename="tcof_Q13428.pdb", length=1488,
    sites=(
        GroundTruthSite(146, "K", "acetylation", "A", (23749302, 23954790, 25953088, 26051181)),
        GroundTruthSite(244, "K", "acetylation", "A", (23954790, 25953088, 26051181)),
        GroundTruthSite(507, "K", "acetylation", "A", (23954790, 25953088, 26051181)),
        GroundTruthSite(746, "K", "acetylation", "A", (23749302, 23954790, 25953088, 26051181)),
        GroundTruthSite(811, "K", "acetylation", "A", (23236377, 23749302, 23954790, 25953088, 26051181)),
        GroundTruthSite(579, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(637, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(904, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(147, "T", "o_linked_glycosylation", "A", (23301498, 27655845, 28510447, 29351928, 31373491, 32119511)),
        GroundTruthSite(510, "S", "o_linked_glycosylation", "A", (23301498, 27655845, 28510447, 29351928, 31373491, 32119511)),
        GroundTruthSite(533, "T", "o_linked_glycosylation", "A", (23301498, 27655845, 28510447, 29351928, 31373491, 32119511)),
        GroundTruthSite(807, "S", "o_linked_glycosylation", "A", (23301498, 27655845, 28510447, 29351928, 31373491, 32119511)),
        GroundTruthSite(814, "T", "o_linked_glycosylation", "A", (23301498, 27655845, 28510447, 29351928, 31373491, 32119511)),
        GroundTruthSite(381, "S", "phosphorylation", "A", (17081983, 17287340, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(906, "S", "phosphorylation", "A", (17081983, 17287340, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(1228, "S", "phosphorylation", "A", (15302935, 16964243, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(1350, "S", "phosphorylation", "A", (17081983, 17287340, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(1378, "S", "phosphorylation", "A", (15302935, 16964243, 17081983, 17287340, 18669648, 19060867)),
        GroundTruthSite(126, "K", "sumoylation", "B", (25772364,)),
        GroundTruthSite(507, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(600, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(637, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(1414, "K", "sumoylation", "A", (25218447, 25772364, 28112733)),
    ),
)

# ---------------------------------------------------------------------------
# 9. Catenina delta-1 / p120-catenin (CTNND1) -- O60716 (968 aa). Verificado
# 2026-08-09. S268 es el sitio mas solido: Xia 2003 (PMID 12885254) lo localiza
# por mapeo triptico 2D + mutagenesis dirigida, caracterizacion dedicada, no
# screening.
# ---------------------------------------------------------------------------
_CTND1 = PanelEntry(
    name="ctnd1", uniprot_accession="O60716", pdb_filename="ctnd1_O60716.pdb", length=968,
    sites=(
        GroundTruthSite(25, "K", "acetylation", "B", (23236377,)),
        GroundTruthSite(355, "K", "acetylation", "B", (20167786,)),
        GroundTruthSite(433, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(574, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(533, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(692, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(749, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(230, "S", "phosphorylation", "A", (18669648, 19060867, 20068231, 21406692, 21712546, 23186163)),
        GroundTruthSite(268, "S", "phosphorylation", "A", (12885254, 17081983, 19060867, 20068231, 21406692, 21712546, 23186163),
                        note="Xia 2003 (Biochemistry) mapea S268 por triptico 2D + mutagenesis dirigida, no es un screening"),
        GroundTruthSite(269, "S", "phosphorylation", "A", (17081983, 18669648, 19060867, 20068231, 21406692, 21712546)),
        GroundTruthSite(349, "S", "phosphorylation", "A", (17081983, 19060867, 20068231, 21406692, 21712546, 23186163)),
        GroundTruthSite(352, "S", "phosphorylation", "A", (17081983, 19060867, 20068231, 21406692, 21712546, 23186163)),
        GroundTruthSite(394, "C", "s_nitrosylation", "B", (19483679,)),
        GroundTruthSite(421, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(517, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(882, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(355, "K", "ubiquitination", "A", (21890473, 21906983, 21963094, 23000965, 23503661, 27667366)),
    ),
)

# ---------------------------------------------------------------------------
# 10. LARP1 (La-related protein 1) -- Q6PKG0 (1096 aa). Verificado 2026-08-09.
# La entrada mas sucia de todo el panel: dbPTM le atribuia PMIDs de transporte
# de calcio, infeccion por Clonorchis sinensis y cirugia de diverticulo de
# Zenker -- 10/31 candidatos descartados por citas completamente ajenas.
# ---------------------------------------------------------------------------
_LARP1 = PanelEntry(
    name="larp1", uniprot_accession="Q6PKG0", pdb_filename="larp1_Q6PKG0.pdb", length=1096,
    sites=(
        GroundTruthSite(152, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(531, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(753, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(864, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(1017, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(148, "S", "o_linked_glycosylation", "A", (23301498, 28657654, 29351928, 31373491, 31637018, 32119511)),
        GroundTruthSite(584, "S", "o_linked_glycosylation", "A", (23301498, 28657654, 29351928, 31373491, 31637018, 32119511)),
        GroundTruthSite(672, "S", "o_linked_glycosylation", "A", (23301498, 28657654, 29351928, 31373491, 31637018, 32119511)),
        GroundTruthSite(770, "T", "o_linked_glycosylation", "A", (23301498, 28657654, 29351928, 31373491, 31637018, 32119511)),
        GroundTruthSite(809, "T", "o_linked_glycosylation", "A", (23301498, 28657654, 29351928, 31373491, 31637018, 32119511)),
        GroundTruthSite(90, "S", "phosphorylation", "A", (17081983, 18669648, 19060867, 21406692, 21712546, 23186163)),
        GroundTruthSite(526, "T", "phosphorylation", "A", (15302935, 16964243, 17081983, 17287340, 18669648, 19060867)),
        GroundTruthSite(548, "S", "phosphorylation", "A", (15302935, 17081983, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(627, "S", "phosphorylation", "A", (15302935, 17081983, 17287340, 18669648, 19060867, 20068231)),
        GroundTruthSite(766, "S", "phosphorylation", "A", (16964243, 17081983, 18669648, 19060867, 20068231, 21406692)),
        GroundTruthSite(260, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(311, "K", "sumoylation", "B", (25218447, 28112733)),
        GroundTruthSite(531, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(703, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(462, "K", "ubiquitination", "A", (21890473, 21906983, 29967540, 33845483)),
        GroundTruthSite(805, "K", "ubiquitination", "A", (21890473, 21906983, 23000965, 27667366, 33845483)),
    ),
)

# ---------------------------------------------------------------------------
# 11. RPRD2 -- Q5VT52 (1461 aa). Verificado 2026-08-09. O-GlcNAc T172/S173/T174
# bajados a tier B: dbPTM da la misma lista de 15 PMIDs a los 3, y son un
# tripeptido T-S-T contiguo -> localizacion ambigua dentro del peptido de MS.
# ---------------------------------------------------------------------------
_RPRD2 = PanelEntry(
    name="rprd2", uniprot_accession="Q5VT52", pdb_filename="rprd2_Q5VT52.pdb", length=1461,
    sites=(
        GroundTruthSite(374, "S", "phosphorylation", "A", (17081983, 18669648, 19367720, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(485, "S", "phosphorylation", "A", (18669648, 19690332, 23186163)),
        GroundTruthSite(593, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692)),
        GroundTruthSite(614, "S", "phosphorylation", "A", (16964243, 18669648, 19690332, 20068231, 23186163)),
        GroundTruthSite(723, "T", "phosphorylation", "A", (17081983, 18669648, 20068231, 23186163, 24275569)),
        GroundTruthSite(862, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(1032, "K", "acetylation", "B", (23954790, 26051181)),
        GroundTruthSite(89, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(228, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(1366, "R", "arg_methylation", "B", (24129315,)),
        GroundTruthSite(172, "T", "o_linked_glycosylation", "B", (20068230, 23301498, 28657654, 30379171), note="peptido T172-S173-T174: localizacion ambigua"),
        GroundTruthSite(173, "S", "o_linked_glycosylation", "B", (20068230, 23301498, 28657654, 30379171), note="peptido T172-S173-T174: localizacion ambigua"),
        GroundTruthSite(174, "T", "o_linked_glycosylation", "B", (20068230, 23301498, 28657654, 30379171), note="peptido T172-S173-T174: localizacion ambigua"),
        GroundTruthSite(396, "S", "o_linked_glycosylation", "B", (20068230, 23301498, 28657654, 30379171)),
        GroundTruthSite(401, "S", "o_linked_glycosylation", "B", (20068230, 23301498, 28657654, 30379171)),
        GroundTruthSite(293, "K", "ubiquitination", "B", (21890473, 29967540)),
        GroundTruthSite(112, "K", "ubiquitination", "B", (29967540,)),
        GroundTruthSite(118, "K", "ubiquitination", "B", (29967540,)),
    ),
)

# ---------------------------------------------------------------------------
# 12. Plakofilina-4 (PKP4) -- Q99569 (1192 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_PKP4 = PanelEntry(
    name="pkp4", uniprot_accession="Q99569", pdb_filename="pkp4_Q99569.pdb", length=1192,
    sites=(
        GroundTruthSite(314, "S", "phosphorylation", "A", (15144186, 17081983, 18220336, 18669648, 18691976, 19369195, 21406692, 23186163)),
        GroundTruthSite(281, "S", "phosphorylation", "A", (15144186, 18669648, 19413330, 23186163)),
        GroundTruthSite(337, "S", "phosphorylation", "A", (18669648, 19369195, 23186163, 24275569)),
        GroundTruthSite(273, "S", "phosphorylation", "A", (15144186, 18669648, 20363803, 21406692, 21712546)),
        GroundTruthSite(518, "K", "ubiquitination", "A", (21139048, 21890473, 21906983, 22053931, 23000965, 29967540)),
        GroundTruthSite(683, "K", "ubiquitination", "A", (21139048, 21890473, 21906983, 22053931, 23000965, 29967540)),
        GroundTruthSite(336, "S", "phosphorylation", "B", (18669648, 20068231, 21406692, 21712546), note="adyacente a S337, posible ambiguedad de localizacion"),
        GroundTruthSite(625, "K", "acetylation", "B", (19608861,)),
        GroundTruthSite(982, "K", "lys_methylation", "B", (23644510,)),
        GroundTruthSite(106, "S", "o_linked_glycosylation", "B", (28657654,), note="UniProt anota S106 como fosfoserina; crosstalk O-GlcNAc/fosfo plausible pero un unico estudio"),
    ),
)

# ---------------------------------------------------------------------------
# 13. ITPRID2 -- P28290 (1259 aa). Verificado 2026-08-09. OJO: P28290 NO es la
# inositol-tetrakisphosphate 1-kinase (esa es ITPKA/Q9Y3D2). P28290 = ITPI2_HUMAN,
# "Protein ITPRID2" (gen ITPRID2, antes SSFA2/KIAA0620). Se mantiene el slug
# "itpi2" porque coincide con el nombre de entrada de UniProt.
# ---------------------------------------------------------------------------
_ITPI2 = PanelEntry(
    name="itpi2", uniprot_accession="P28290", pdb_filename="itpi2_P28290.pdb", length=1259,
    sites=(
        GroundTruthSite(92, "S", "phosphorylation", "A", (18669648, 19369195, 19690332, 20068231, 23186163, 24275569)),
        GroundTruthSite(668, "S", "phosphorylation", "A", (18669648, 20363803, 21712546, 23186163)),
        GroundTruthSite(737, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 24275569)),
        GroundTruthSite(739, "S", "phosphorylation", "A", (17081983, 18669648, 23186163, 24275569)),
        GroundTruthSite(746, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 24275569)),
        GroundTruthSite(808, "K", "sumoylation", "A", (25218447, 28112733)),
        GroundTruthSite(360, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(1250, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(77, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(102, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(308, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(445, "K", "acetylation", "B", (26051181,)),
    ),
)

# ---------------------------------------------------------------------------
# 14. BCLAF1 -- Q9NYF8 (920 aa). Verificado 2026-08-09. La proteina mejor
# documentada del panel dbPTM: 20/22 candidatos sobrevivieron.
# ---------------------------------------------------------------------------
_BCLF1 = PanelEntry(
    name="bclf1", uniprot_accession="Q9NYF8", pdb_filename="bclf1_Q9NYF8.pdb", length=920,
    sites=(
        GroundTruthSite(177, "S", "phosphorylation", "A", (17081983, 18318008, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(385, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(397, "S", "phosphorylation", "A", (17081983, 17487921, 18318008, 18669648, 19367720, 20068231, 21406692, 23186163)),
        GroundTruthSite(512, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692, 23186163)),
        GroundTruthSite(658, "S", "phosphorylation", "A", (18669648, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(421, "K", "sumoylation", "A", (25755297, 28112733)),
        GroundTruthSite(437, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(491, "K", "sumoylation", "A", (25218447, 25772364, 28112733)),
        GroundTruthSite(501, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(676, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733)),
        GroundTruthSite(421, "K", "acetylation", "A", (19608861, 23236377, 23749302, 25953088, 26051181, 27452117)),
        GroundTruthSite(891, "K", "acetylation", "A", (21339330, 23236377, 23954790, 26051181, 26822725)),
        GroundTruthSite(335, "K", "acetylation", "A", (19608861, 25825284, 25953088, 26051181)),
        GroundTruthSite(421, "K", "ubiquitination", "A", (21890473, 21906983, 29967540, 33845483)),
        GroundTruthSite(332, "K", "ubiquitination", "A", (21890473, 21906983, 29967540, 33845483)),
        GroundTruthSite(333, "T", "o_linked_glycosylation", "A", (21740066, 23301498, 28510447, 30379171)),
        GroundTruthSite(332, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(445, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(809, "R", "arg_methylation", "B", (24129315,)),
        GroundTruthSite(593, "K", "malonylation", "B", (26320211,)),
    ),
)

# ---------------------------------------------------------------------------
# 15. SORBS1 / SRBS1 -- Q9BX66 (1292 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_SRBS1 = PanelEntry(
    name="srbs1", uniprot_accession="Q9BX66", pdb_filename="srbs1_Q9BX66.pdb", length=1292,
    sites=(
        GroundTruthSite(350, "S", "phosphorylation", "A", (18669648, 21406692, 23186163, 24275569)),
        GroundTruthSite(665, "S", "phosphorylation", "A", (18669648, 21130716, 23186163, 24275569)),
        GroundTruthSite(349, "T", "phosphorylation", "B", (18669648, 20166139, 21406692, 21712546), note="adyacente a S350, posible ambiguedad de localizacion"),
        GroundTruthSite(479, "S", "o_linked_glycosylation", "B", (28657654, 30379171)),
        GroundTruthSite(708, "T", "o_linked_glycosylation", "B", (28657654, 30379171), note="UniProt anota T708 como fosfotreonina; crosstalk O-GlcNAc/fosfo plausible"),
        GroundTruthSite(943, "S", "o_linked_glycosylation", "B", (28657654, 30379171)),
        GroundTruthSite(1048, "S", "o_linked_glycosylation", "B", (28657654, 30379171)),
        GroundTruthSite(1130, "K", "acetylation", "B", (20167786,)),
    ),
)

# ---------------------------------------------------------------------------
# 16. UBAP2L -- Q14157 (1087 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_UBP2L = PanelEntry(
    name="ubp2l", uniprot_accession="Q14157", pdb_filename="ubp2l_Q14157.pdb", length=1087,
    sites=(
        GroundTruthSite(609, "S", "phosphorylation", "A", (18669648, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(608, "S", "phosphorylation", "A", (18669648, 23186163)),
        GroundTruthSite(605, "S", "phosphorylation", "A", (18669648, 23186163)),
        GroundTruthSite(467, "S", "phosphorylation", "A", (18669648, 18691976, 19369195, 21406692)),
        GroundTruthSite(453, "S", "phosphorylation", "B", (18669648, 20068231, 23186163), note="S453/S454 ambiguo (doblete Ser); UniProt cura 454"),
        GroundTruthSite(612, "K", "acetylation", "B", (22424773, 23954790, 26051181)),
        GroundTruthSite(101, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(53, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(63, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(445, "S", "o_linked_glycosylation", "B", (20068230, 20305658, 21740066, 22661428, 23301498, 27114449)),
        GroundTruthSite(260, "T", "o_linked_glycosylation", "B", (20068230, 20305658, 21740066, 22661428, 23301498, 27114449)),
        GroundTruthSite(262, "S", "o_linked_glycosylation", "B", (20068230, 20305658, 21740066, 22661428, 23301498, 27114449)),
        GroundTruthSite(277, "T", "o_linked_glycosylation", "B", (20068230, 20305658, 21740066, 22661428, 23301498, 27114449)),
        GroundTruthSite(317, "S", "o_linked_glycosylation", "B", (20068230, 20305658, 21740066, 22661428, 23301498, 27114449)),
        GroundTruthSite(353, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(420, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(864, "K", "ubiquitination", "B", (21890473, 21906983)),
    ),
)

# ---------------------------------------------------------------------------
# 17. THRAP3 / TR150 -- Q9Y2W1 (955 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_TR150 = PanelEntry(
    name="tr150", uniprot_accession="Q9Y2W1", pdb_filename="tr150_Q9Y2W1.pdb", length=955,
    sites=(
        GroundTruthSite(682, "S", "phosphorylation", "A", (18669648, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(253, "S", "phosphorylation", "A", (17081983, 18220336, 18669648, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(248, "S", "phosphorylation", "A", (17081983, 18220336, 18669648, 20068231, 23186163, 24275569)),
        GroundTruthSite(379, "S", "phosphorylation", "A", (18669648, 20068231, 21406692, 23186163)),
        GroundTruthSite(939, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692, 23186163)),
        GroundTruthSite(486, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733)),
        GroundTruthSite(470, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(705, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(711, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(396, "K", "sumoylation", "A", (25755297, 28112733)),
        GroundTruthSite(387, "K", "acetylation", "B", (23236377, 23749302, 25953088, 26051181)),
        GroundTruthSite(401, "K", "acetylation", "B", (19608861, 22424773, 23749302, 25953088)),
        GroundTruthSite(470, "K", "acetylation", "B", (23749302, 25953088, 26051181, 26822725)),
        GroundTruthSite(420, "K", "acetylation", "B", (19608861, 23749302, 25953088, 26051181)),
        GroundTruthSite(527, "K", "acetylation", "B", (22424773, 23749302, 25953088, 26051181)),
        GroundTruthSite(66, "R", "arg_methylation", "B", (24129315,)),
        GroundTruthSite(252, "K", "lys_methylation", "B", (24129315,)),
        GroundTruthSite(519, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(390, "S", "o_linked_glycosylation", "B", (29351928, 30379171, 32119511)),
        GroundTruthSite(916, "S", "o_linked_glycosylation", "B", (29351928, 30379171, 32119511)),
        GroundTruthSite(221, "K", "ubiquitination", "B", (21890473, 23000965, 29967540)),
        GroundTruthSite(519, "K", "ubiquitination", "B", (21890473, 27667366, 29967540, 33845483)),
        GroundTruthSite(470, "K", "ubiquitination", "B", (21890473, 21906983, 27667366, 33845483)),
        GroundTruthSite(215, "K", "ubiquitination", "B", (23000965, 24816145, 29967540, 33845483)),
        GroundTruthSite(451, "K", "ubiquitination", "B", (21890473, 21906983, 29967540)),
    ),
)

# ---------------------------------------------------------------------------
# 18. SCAF11 / SCAFB -- Q99590 (1463 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_SCAFB = PanelEntry(
    name="scafb", uniprot_accession="Q99590", pdb_filename="scafb_Q99590.pdb", length=1463,
    sites=(
        GroundTruthSite(796, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(608, "S", "phosphorylation", "A", (18220336, 18669648, 19690332, 20068231, 23186163, 24275569)),
        GroundTruthSite(802, "S", "phosphorylation", "A", (18669648, 19690332, 21406692, 23186163, 24275569)),
        GroundTruthSite(776, "S", "phosphorylation", "A", (20068231, 21406692, 23186163)),
        GroundTruthSite(614, "S", "phosphorylation", "A", (18669648, 20068231, 23186163)),
        GroundTruthSite(601, "K", "sumoylation", "A", (25218447, 25772364, 28112733)),
        GroundTruthSite(1178, "K", "sumoylation", "A", (25218447, 25772364, 28112733)),
        GroundTruthSite(610, "K", "sumoylation", "A", (25218447, 28112733)),
        GroundTruthSite(676, "K", "sumoylation", "A", (25114211, 25218447, 28112733), note="UniProt aporta 2 PMIDs mas que el borrador (SUMO1 y SUMO2)"),
        GroundTruthSite(596, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(1392, "K", "acetylation", "B", (19608861, 22424773, 23954790, 26051181)),
        GroundTruthSite(1354, "K", "acetylation", "B", (23749302, 25953088, 26051181)),
        GroundTruthSite(676, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(1439, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(459, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(1428, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(344, "S", "o_linked_glycosylation", "B", (23301498, 29351928, 30379171)),
        GroundTruthSite(357, "S", "o_linked_glycosylation", "B", (23301498, 29351928, 30379171)),
        GroundTruthSite(512, "S", "o_linked_glycosylation", "B", (23301498, 29351928, 30379171)),
        GroundTruthSite(680, "S", "o_linked_glycosylation", "B", (23301498, 29351928, 30379171)),
        GroundTruthSite(882, "S", "o_linked_glycosylation", "B", (23301498, 29351928, 30379171)),
        GroundTruthSite(120, "K", "ubiquitination", "B", (21963094, 27667366, 29967540, 33845483)),
        GroundTruthSite(146, "K", "ubiquitination", "B", (21963094, 27667366, 29967540, 33845483)),
        GroundTruthSite(113, "K", "ubiquitination", "B", (21906983, 21963094, 27667366, 29967540)),
        GroundTruthSite(165, "K", "ubiquitination", "B", (21906983, 21963094, 27667366, 29967540)),
        GroundTruthSite(1126, "K", "ubiquitination", "B", (23000965, 29967540, 33845483)),
    ),
)

# ---------------------------------------------------------------------------
# 19. ZMYND8 / ZMYD8 -- Q9ULU4 (1186 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_ZMYD8 = PanelEntry(
    name="zmyd8", uniprot_accession="Q9ULU4", pdb_filename="zmyd8_Q9ULU4.pdb", length=1186,
    sites=(
        GroundTruthSite(490, "S", "phosphorylation", "A", (18669648, 19690332, 21406692, 24275569)),
        GroundTruthSite(486, "S", "phosphorylation", "A", (18669648, 19690332, 20068231)),
        GroundTruthSite(406, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 23186163, 24275569)),
        GroundTruthSite(547, "S", "phosphorylation", "A", (18669648, 21406692, 23186163, 24275569)),
        GroundTruthSite(488, "S", "phosphorylation", "B", (18669648, 20068231, 21406692), note="S486/S488/S490 ambiguos; UniProt cura 486 y 490, no 488"),
        GroundTruthSite(611, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733)),
        GroundTruthSite(645, "K", "sumoylation", "A", (25218447, 25755297, 28112733)),
        GroundTruthSite(56, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(70, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(390, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(562, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(175, "K", "acetylation", "B", (25953088,)),
        GroundTruthSite(200, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(750, "T", "o_linked_glycosylation", "B", (22661428, 23301498, 30379171, 32119511, 32574038)),
        GroundTruthSite(356, "K", "ubiquitination", "B", (21890473,)),
    ),
)

# ---------------------------------------------------------------------------
# 20. ACIN1 / ACINU -- Q9UKV3 (1341 aa). Verificado 2026-08-09. Unica de las 25
# nuevas sin un solo PMID falso o ajeno en su borrador. K654 lys_methylation es
# tier A por caracterizacion dedicada: Rathert 2008 (Nat Chem Biol) nombra
# ACINUS explicitamente como diana no-histonica de G9a/EHMT2.
# ---------------------------------------------------------------------------
_ACINU = PanelEntry(
    name="acinu", uniprot_accession="Q9UKV3", pdb_filename="acinu_Q9UKV3.pdb", length=1341,
    sites=(
        GroundTruthSite(1004, "S", "phosphorylation", "A", (17081983, 18669648, 19367720, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(490, "S", "phosphorylation", "A", (17081983, 18669648, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(710, "S", "phosphorylation", "A", (17081983, 18669648, 19367720, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(216, "S", "phosphorylation", "A", (17081983, 18318008, 18669648, 19367720, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(410, "S", "phosphorylation", "A", (18669648, 19690332, 20068231, 21406692, 23186163, 24275569)),
        GroundTruthSite(654, "K", "lys_methylation", "A", (18438403,), note="G9a/EHMT2 -- Rathert et al. Nat Chem Biol 2008, ACINUS nombrado en el abstract"),
        GroundTruthSite(532, "K", "sumoylation", "A", (25114211, 25218447, 25755297, 25772364, 28112733)),
        GroundTruthSite(315, "K", "sumoylation", "A", (25218447, 25772364, 28112733)),
        GroundTruthSite(305, "K", "sumoylation", "A", (25755297, 28112733)),
        GroundTruthSite(268, "K", "sumoylation", "B", (25114211,)),
        GroundTruthSite(375, "K", "sumoylation", "B", (28112733,)),
        GroundTruthSite(238, "K", "acetylation", "B", (23236377, 23954790, 25953088, 26051181)),
        GroundTruthSite(359, "K", "acetylation", "B", (19608861, 23749302, 25953088, 26051181)),
        GroundTruthSite(516, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(532, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(860, "K", "acetylation", "B", (25953088, 26051181)),
        GroundTruthSite(546, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(733, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(1052, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(1083, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(1223, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(81, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(1049, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(657, "S", "o_linked_glycosylation", "B", (23301498, 26853435, 28510447, 29351928, 30379171, 32119511)),
        GroundTruthSite(661, "S", "o_linked_glycosylation", "B", (23301498, 26853435, 28510447, 29351928, 30379171, 32119511)),
        GroundTruthSite(870, "S", "o_linked_glycosylation", "B", (23301498, 26853435, 28510447, 29351928, 30379171, 32119511)),
        GroundTruthSite(5, "K", "ubiquitination", "B", (17370265,)),
        GroundTruthSite(103, "K", "ubiquitination", "B", (21890473,)),
    ),
)

# ---------------------------------------------------------------------------
# 21. SORBS2 / SRBS2 -- O94875 (1100 aa). Verificado 2026-08-09. Unica de las 25
# nuevas sin ningun descarte: los 5 fosfositios estan curados por UniProt
# (ECO:0000269) en la posicion exacta.
# ---------------------------------------------------------------------------
_SRBS2 = PanelEntry(
    name="srbs2", uniprot_accession="O94875", pdb_filename="srbs2_O94875.pdb", length=1100,
    sites=(
        GroundTruthSite(259, "S", "phosphorylation", "A", (18669648, 20068231, 24275569, 20873877, 21712546, 20068230)),
        GroundTruthSite(277, "T", "phosphorylation", "A", (18669648, 20068231, 24275569, 20873877, 21712546, 20068230)),
        GroundTruthSite(298, "S", "phosphorylation", "A", (18669648, 20068231, 24275569, 20873877, 21712546, 21130716)),
        GroundTruthSite(299, "S", "phosphorylation", "A", (18669648, 20068231, 24275569, 20873877, 21712546, 17929957)),
        GroundTruthSite(304, "S", "phosphorylation", "A", (18669648, 20068231, 24275569, 20873877, 21712546, 17929957)),
        GroundTruthSite(288, "S", "o_linked_glycosylation", "B", (23271734, 28657654, 30379171), note="solo MS de O-GlcNAc, sin curacion UniProt"),
        GroundTruthSite(290, "T", "o_linked_glycosylation", "B", (23271734, 28657654, 30379171), note="solo MS de O-GlcNAc, sin curacion UniProt"),
        GroundTruthSite(68, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(875, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(880, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(919, "K", "acetylation", "B", (26051181,)),
        GroundTruthSite(925, "K", "acetylation", "B", (26051181,)),
    ),
)

# ---------------------------------------------------------------------------
# 22. MYPT1 / PPP1R12A -- O14974 (1030 aa). Verificado 2026-08-09. N67/N100/N226
# hydroxylation: UniProt anota (3S)-3-hidroxiasparagina "by HIF1AN" con
# ECO:0000269 en las 3 posiciones, Webb 2009 confirma hidroxilacion dependiente
# de FIH -- caracterizacion dedicada, no screening.
# ---------------------------------------------------------------------------
_MYPT1 = PanelEntry(
    name="mypt1", uniprot_accession="O14974", pdb_filename="mypt1_O14974.pdb", length=1030,
    sites=(
        GroundTruthSite(696, "T", "phosphorylation", "A", (11719507, 15194681, 15723050, 19701943, 10601309, 18669648),
                        note="sitio inhibitorio de ROCK/ZIPK; Feng 1999 (10601309) lo numera Thr695 (gizzard), O14974 humano = Thr696"),
        GroundTruthSite(445, "S", "phosphorylation", "A", (20354225, 18669648, 20068231, 23186163, 24275569, 21406692), note="NUAK1 (Zagorska 2010)"),
        GroundTruthSite(443, "T", "phosphorylation", "A", (18669648, 20068231, 23186163, 24275569, 21406692, 19690332)),
        GroundTruthSite(299, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 21406692, 20873877, 21712546)),
        GroundTruthSite(507, "S", "phosphorylation", "A", (18669648, 23186163, 24275569, 21406692, 19690332, 20873877)),
        GroundTruthSite(67, "N", "hydroxylation", "A", (19245366,), note="(3S)-3-hidroxiasparagina por FIH/HIF1AN, Webb 2009; UniProt ECO:0000269"),
        GroundTruthSite(100, "N", "hydroxylation", "A", (19245366,), note="idem N67"),
        GroundTruthSite(226, "N", "hydroxylation", "A", (19245366,), note="idem N67"),
        GroundTruthSite(379, "S", "o_linked_glycosylation", "A", (18840611, 20305658, 28657654, 30379171, 21740066, 22967762)),
        GroundTruthSite(381, "T", "o_linked_glycosylation", "A", (18840611, 20305658, 28657654, 30379171, 21740066, 22967762)),
        GroundTruthSite(385, "T", "o_linked_glycosylation", "A", (18840611, 20305658, 28657654, 30379171, 21740066, 22967762)),
        GroundTruthSite(585, "S", "o_linked_glycosylation", "A", (18840611, 20305658, 28657654, 30379171, 21740066, 22967762)),
        GroundTruthSite(601, "S", "o_linked_glycosylation", "A", (18840611, 20305658, 28657654, 30379171, 21740066, 22967762),
                        note="OJO: 601 tambien es fosfoserina curada (18477460); aqui el sitio es O-GlcNAc"),
        GroundTruthSite(315, "K", "acetylation", "B", (26051181, 25953088)),
        GroundTruthSite(420, "K", "acetylation", "B", (25953088,)),
        GroundTruthSite(947, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(984, "K", "malonylation", "B", (26320211,)),
    ),
)

# ---------------------------------------------------------------------------
# 23. UBE4B -- O95155 (1302 aa). Verificado 2026-08-09.
# ---------------------------------------------------------------------------
_UBE4B = PanelEntry(
    name="ube4b", uniprot_accession="O95155", pdb_filename="ube4b_O95155.pdb", length=1302,
    sites=(
        GroundTruthSite(803, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 17081983, 21406692, 20873877)),
        GroundTruthSite(88, "S", "phosphorylation", "A", (18669648, 23186163, 17081983, 24275569, 21406692, 19690332)),
        GroundTruthSite(103, "S", "phosphorylation", "A", (23186163, 21406692, 21712546, 19651622, 21130716, 21955146)),
        GroundTruthSite(105, "S", "phosphorylation", "A", (18669648, 23186163, 21406692, 15302935, 20873877, 21712546)),
        GroundTruthSite(801, "K", "ubiquitination", "A", (21890473, 21906983, 21963094, 21987572, 22505724, 23503661)),
        GroundTruthSite(794, "K", "acetylation", "B", (25953088, 23236377)),
        GroundTruthSite(555, "K", "acetylation", "B", (25953088,)),
        GroundTruthSite(794, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(1106, "K", "malonylation", "B", (26320211,)),
    ),
)

# ---------------------------------------------------------------------------
# 24. HSP90-alpha / HSP90AA1 -- P07900 (732 aa). Verificado 2026-08-09. C598
# s_nitrosylation (no C597): dbPTM asignaba el sitio a C597, pero UniProt asigna
# el mismo paper (Martinez-Ruiz 2005) a Cys598 -- par CC adyacente donde el MS
# no discrimina, se conserva el residuo que UniProt cura.
# ---------------------------------------------------------------------------
_HS90A = PanelEntry(
    name="hs90a", uniprot_accession="P07900", pdb_filename="hs90a_P07900.pdb", length=732,
    sites=(
        GroundTruthSite(231, "S", "phosphorylation", "A", (2492519, 18318008, 20068231, 23186163, 17081983, 24275569), note="sitio de CK2, Lees-Miller 1989 nombra Ser231 de la isoforma alpha"),
        GroundTruthSite(263, "S", "phosphorylation", "A", (2492519, 16807684, 18669648, 20068231, 23186163, 17081983), note="sitio de CK2, Lees-Miller 1989 nombra Ser263 de la isoforma alpha"),
        GroundTruthSite(252, "S", "phosphorylation", "A", (20068231, 23186163, 17081983, 21406692, 18088087, 19651622)),
        GroundTruthSite(315, "S", "phosphorylation", "B", (18669648, 20068231, 23186163, 20873877, 21712546, 19651622), note="solo MS; sin curacion UniProt y con riesgo de peptido compartido con HSP90AB1"),
        GroundTruthSite(317, "T", "phosphorylation", "B", (20068231, 23186163, 20873877, 21712546, 22199227, 22617229), note="idem S315"),
        GroundTruthSite(283, "K", "acetylation", "A", (19608861, 25953088, 23236377, 23749302, 22424773, 27452117)),
        GroundTruthSite(292, "K", "acetylation", "A", (19608861, 26051181, 25953088, 23236377, 23749302, 22424773), note="no confundir con K294 (Scroggins 2007, 17218278), que es el sitio caracterizado adyacente"),
        GroundTruthSite(362, "K", "acetylation", "A", (19608861, 26051181, 25953088, 23236377, 22424773, 27452117)),
        GroundTruthSite(69, "K", "acetylation", "A", (19608861, 26051181, 25953088, 23236377, 22424773, 27452117)),
        GroundTruthSite(576, "K", "acetylation", "A", (19608861, 26051181, 25953088, 23236377, 23749302, 27452117)),
        GroundTruthSite(598, "C", "s_nitrosylation", "A", (15937123,), note="Martinez-Ruiz 2005; UniProt asigna el S-nitrosocisteina a Cys598 (no Cys597)"),
        GroundTruthSite(420, "C", "s_nitrosylation", "B", (20140087,)),
        GroundTruthSite(481, "C", "s_nitrosylation", "B", (17629318,)),
        GroundTruthSite(420, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(481, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(529, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(598, "C", "glutathionylation", "B", (22555962,)),
        GroundTruthSite(615, "K", "lys_methylation", "B", (23644510, 24129315)),
        GroundTruthSite(407, "K", "malonylation", "B", (21908771, 26320211)),
        GroundTruthSite(74, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(292, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(458, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(489, "K", "malonylation", "B", (26320211,)),
        GroundTruthSite(100, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(314, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(436, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(539, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(546, "K", "succinylation", "B", (23954790,)),
        GroundTruthSite(185, "K", "ubiquitination", "A", (21890473, 21906983, 21139048)),
        GroundTruthSite(69, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(224, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(283, "K", "ubiquitination", "B", (21890473, 21906983)),
        GroundTruthSite(292, "K", "ubiquitination", "B", (21890473, 21906983)),
    ),
)

# ---------------------------------------------------------------------------
# 25. RSF1 / HBXAP -- Q96T23 (1441 aa). Verificado 2026-08-09. Las 5
# sumoilaciones son el mejor bloque: UniProt las anota como cross-link SUMO2 con
# exactamente los 4 PMIDs de proteomica SUMO sitio-especifica del borrador.
# ---------------------------------------------------------------------------
_RSF1 = PanelEntry(
    name="rsf1", uniprot_accession="Q96T23", pdb_filename="rsf1_Q96T23.pdb", length=1441,
    sites=(
        GroundTruthSite(1345, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 17081983, 24275569, 21406692)),
        GroundTruthSite(604, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 17081983, 24275569, 21406692)),
        GroundTruthSite(622, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 17081983, 24275569, 21406692)),
        GroundTruthSite(1375, "S", "phosphorylation", "A", (18669648, 20068231, 23186163, 17081983, 21406692, 19690332)),
        GroundTruthSite(1305, "T", "phosphorylation", "A", (18669648, 20068231, 23186163, 24275569, 21406692, 19690332)),
        GroundTruthSite(254, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733), note="cross-link SUMO2 curado en UniProt (ECO:0000269)"),
        GroundTruthSite(294, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733), note="cross-link SUMO2 curado en UniProt (ECO:0000269)"),
        GroundTruthSite(309, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733), note="cross-link SUMO2 curado en UniProt (ECO:0000269)"),
        GroundTruthSite(323, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733), note="cross-link SUMO2 curado en UniProt (ECO:0000269)"),
        GroundTruthSite(419, "K", "sumoylation", "A", (25218447, 25755297, 25772364, 28112733), note="cross-link SUMO2 curado en UniProt (ECO:0000269)"),
        GroundTruthSite(1386, "K", "acetylation", "A", (26051181, 25953088, 23749302, 22424773, 27452117, 21339330)),
        GroundTruthSite(1378, "K", "acetylation", "A", (26051181, 25953088, 23749302, 23954790)),
        GroundTruthSite(1390, "K", "acetylation", "A", (26051181, 25953088, 23236377, 23749302)),
        GroundTruthSite(136, "K", "acetylation", "B", (20167786, 26051181)),
        GroundTruthSite(698, "K", "ubiquitination", "B", (21890473, 21906983), note="UniProt anota este K como cross-link SUMO2 (28112733); el resto diGly no discrimina Ub/SUMO"),
        GroundTruthSite(1039, "K", "ubiquitination", "B", (21890473, 21906983), note="idem K698"),
    ),
)

# ---------------------------------------------------------------------------
# 26. KIT ligand / SCF -- P21583 (273 aa). Control NEGATIVO principal, preservado
# intacto desde 2026-08-04: Asn-97 es un sequon N-X-S perfectamente valido que la
# propia UniProt (citando PMID 1381905) reporta EXPLICITAMENTE como NO
# glicosilado en ninguna de las 2 isoformas conocidas -- distingue un motor con
# especificidad real de uno que solo empareja el motivo. dbPTM no tiene
# anotaciones de "confirmado como NO modificado", asi que este control no puede
# salir de la migracion a dbPTM y se mantiene como excepcion manual.
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
    _P53, _HISTONE_H3, _ENOA, _KNG1, _FA7, _MAP4, _NU153, _TCOF, _CTND1, _LARP1,
    _RPRD2, _PKP4, _ITPI2, _BCLF1, _SRBS1, _UBP2L, _TR150, _SCAFB, _ZMYD8, _ACINU,
    _SRBS2, _MYPT1, _UBE4B, _HS90A, _RSF1, _KIT_LIGAND_SCF,
)
