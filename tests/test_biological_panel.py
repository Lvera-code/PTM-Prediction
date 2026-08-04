"""Verifica la integridad del panel de validacion biologica (punto 8 del plan de robustez
post-demo-prep) contra las estructuras reales descargadas -- no solo consistencia interna
de los datos en Python, sino que cada ``GroundTruthSite.residue`` coincide realmente con el
residuo de esa posicion en el PDB real de AlphaFold correspondiente. Esto es exactamente el
tipo de chequeo que habria detectado cualquiera de las 3 trampas de numeracion documentadas
en STATUS.md (histonas/protrombina/EPO, PDB != UniProt) si el panel las hubiera heredado por
error -- como usa AlphaFold (numeracion == UniProt, verificado), todas deben coincidir.
"""

import pytest

from src.config.settings import Settings
from src.utils.structure_parser import parse_structure
from src.validation.biological_panel import PANEL


@pytest.mark.parametrize("entry", PANEL, ids=lambda e: e.name)
def test_pdb_existe_y_tiene_la_longitud_esperada(entry, tmp_path):
    assert entry.pdb_path.is_file(), f"falta el PDB del panel: {entry.pdb_path}"
    record = parse_structure(entry.pdb_path, tmp_path)
    assert len(record.sequence) == entry.length


@pytest.mark.parametrize("entry", PANEL, ids=lambda e: e.name)
def test_cada_sitio_coincide_con_el_residuo_real_del_pdb(entry, tmp_path):
    record = parse_structure(entry.pdb_path, tmp_path)
    for site in entry.sites:
        real_residue = record.sequence[site.position - 1]
        assert real_residue == site.residue, (
            f"{entry.name} posicion {site.position}: se esperaba '{site.residue}' "
            f"({site.ptm_type}) pero el PDB real tiene '{real_residue}' -- posible error de "
            "numeracion (ver docstring de biological_panel.py, 3 trampas de offset conocidas)"
        )


@pytest.mark.parametrize("entry", PANEL, ids=lambda e: e.name)
def test_todos_los_tipos_ptm_son_de_los_17_soportados(entry):
    for site in entry.sites:
        assert site.ptm_type in Settings.DEEPPTMPRED_PTM_TYPES, (
            f"{entry.name} posicion {site.position}: '{site.ptm_type}' no es uno de los 17 "
            "tipos que el pipeline soporta"
        )


@pytest.mark.parametrize("entry", PANEL, ids=lambda e: e.name)
def test_sin_positivos_duplicados_en_posicion_y_tipo(entry):
    keys = [(s.position, s.ptm_type) for s in entry.positives]
    assert len(keys) == len(set(keys)), f"{entry.name} tiene sitios positivos duplicados"


def test_panel_cubre_los_9_tipos_de_fase_a():
    fase_a_types = set(Settings.FASE_A_SUPPORTED_PTM_TYPES)
    covered = {site.ptm_type for entry in PANEL for site in entry.positives}
    faltantes = fase_a_types - covered
    assert not faltantes, f"el panel no cubre estos tipos de Fase A: {sorted(faltantes)}"


def test_hay_al_menos_un_control_negativo_real():
    negativos = [site for entry in PANEL for site in entry.negatives]
    assert len(negativos) >= 1
    assert all(site.tier == "A" for site in negativos), "los negativos deben ser de maxima confianza"
