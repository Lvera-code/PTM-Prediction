"""Tests dedicados de src/structural/pyrosetta_glycan_patch.py (punto 6 del plan de
robustez post-demo-prep). ``check_glygen_evidence`` no importa pyrosetta en absoluto (solo
``src.structural.glygen_client``, un cliente HTTP puro) -- totalmente testeable sin stub de
pyrosetta, mockeando ``lookup_site`` en el limite del cliente real (mismo criterio que
tests/test_glygen_client.py, que ya cubre lookup_site en si; aqui se cubre como
check_glygen_evidence traduce sus 4 resultados posibles a un mensaje informativo).
"""

from src.structural import glygen_client, pyrosetta_glycan_patch


def test_glygen_no_disponible_degrada_a_aviso(monkeypatch):
    def _boom(accession, position, ptm_type):
        raise glygen_client.GlyGenLookupError("timeout de red")

    monkeypatch.setattr(glygen_client, "lookup_site", _boom)

    message = pyrosetta_glycan_patch.check_glygen_evidence("P10636", 484, "n_linked_glycosylation")

    assert "GlyGen no disponible" in message
    assert "timeout de red" in message


def test_glygen_sin_sitio_reportado(monkeypatch):
    monkeypatch.setattr(glygen_client, "lookup_site", lambda accession, position, ptm_type: None)

    message = pyrosetta_glycan_patch.check_glygen_evidence("P10636", 484, "n_linked_glycosylation")

    assert "no reporta ningun sitio" in message


def test_glygen_con_evidencia_experimental_y_glicano_especifico(monkeypatch):
    monkeypatch.setattr(
        glygen_client, "lookup_site",
        lambda accession, position, ptm_type: {
            "site_category": "reported_with_glycan", "glytoucan_ac": "G12345AB",
        },
    )

    message = pyrosetta_glycan_patch.check_glygen_evidence("P10636", 484, "n_linked_glycosylation")

    assert "SI reporta evidencia experimental" in message
    assert "G12345AB" in message


def test_glygen_con_sitio_pero_sin_glicano_especifico(monkeypatch):
    monkeypatch.setattr(
        glygen_client, "lookup_site",
        lambda accession, position, ptm_type: {"site_category": "reported", "glytoucan_ac": None},
    )

    message = pyrosetta_glycan_patch.check_glygen_evidence("P10636", 484, "n_linked_glycosylation")

    assert "sin glicano especifico" in message
    assert "reported" in message


def test_glycan_tree_by_type_coincide_con_supported_types():
    assert pyrosetta_glycan_patch.SUPPORTED_PTM_TYPES == frozenset(pyrosetta_glycan_patch.GLYCAN_TREE_BY_TYPE)
    assert pyrosetta_glycan_patch.SUPPORTED_PTM_TYPES == frozenset(pyrosetta_glycan_patch.TARGET_RESIDUE)
