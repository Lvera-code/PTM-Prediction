"""Tests de scripts/select_recall_subset.py contra fixtures sinteticos chicos.

Nunca golpea dbPTM/UniProt/AlphaFold reales -- la corrida real (25 candidatas,
15/17 tipos cubiertos, tope de sitios por tipo evitando el caso real
encontrado 2026-08-09 de 549 sitios para una sola proteina) se hizo a mano al
implementar.
"""

import json
import sqlite3

from scripts.select_recall_subset import (
    MAX_SITES_PER_ACCESSION_TYPE,
    emit_python,
    load_candidate_sites,
    select_recall_subset,
)


def _build_fixture_db(db_path, rows):
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        """
        CREATE TABLE dbptm_sites (
            accession TEXT NOT NULL, position INTEGER NOT NULL, ptm_type TEXT NOT NULL,
            residue TEXT NOT NULL, tier TEXT NOT NULL, pmids TEXT NOT NULL,
            PRIMARY KEY (accession, position, ptm_type)
        )
        """
    )
    conn.executemany("INSERT INTO dbptm_sites VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    conn.close()


def test_select_recall_subset_prioriza_cobertura_de_tipos_nuevos():
    # A cubre {phospho, acetyl}, B cubre {ubiq}, C cubre {phospho} (redundante).
    tier_a_types = {"A": {"phosphorylation", "acetylation"}, "B": {"ubiquitination"}, "C": {"phosphorylation"}}
    tier_a_counts = {"A": 5, "B": 2, "C": 10}
    lengths = {"A": ("A_HUMAN", 100), "B": ("B_HUMAN", 100), "C": ("C_HUMAN", 100)}

    selected, covered = select_recall_subset(tier_a_types, tier_a_counts, lengths, target_size=2, max_size=2)

    assert selected[0] == "A"  # cubre 2 tipos nuevos, gana pese a menos evidencia que C
    assert "B" in selected  # unica fuente de ubiquitination
    assert covered == {"phosphorylation", "acetylation", "ubiquitination"}


def test_select_recall_subset_excluye_proteinas_sin_tier_a():
    tier_a_types = {"A": set(), "B": {"phosphorylation"}}
    tier_a_counts = {"B": 3}
    lengths = {"A": ("A_HUMAN", 100), "B": ("B_HUMAN", 100)}

    selected, covered = select_recall_subset(tier_a_types, tier_a_counts, lengths, target_size=1, max_size=5)
    assert selected == ["B"]


def test_select_recall_subset_excluye_proteinas_muy_largas():
    tier_a_types = {"A": {"phosphorylation"}, "B": {"phosphorylation"}}
    tier_a_counts = {"A": 1, "B": 1}
    lengths = {"A": ("A_HUMAN", 100), "B": ("B_HUMAN", 5000)}  # B excede MAX_PROTEIN_LENGTH

    selected, covered = select_recall_subset(tier_a_types, tier_a_counts, lengths, target_size=1, max_size=5)
    assert selected == ["A"]


def test_select_recall_subset_se_detiene_al_agotar_elegibles():
    tier_a_types = {"A": {"phosphorylation"}}
    tier_a_counts = {"A": 1}
    lengths = {"A": ("A_HUMAN", 100)}

    selected, covered = select_recall_subset(tier_a_types, tier_a_counts, lengths, target_size=10, max_size=30)
    assert selected == ["A"]  # no hay mas candidatas elegibles, no se cuelga ni inventa


def test_load_candidate_sites_topa_por_tipo_prefiriendo_tier_a_y_mas_pmids(tmp_path):
    db_path = tmp_path / "lookup.sqlite3"
    rows = []
    # 8 sitios tier A de fosforilacion para P00000 (mas que MAX_SITES_PER_ACCESSION_TYPE=5)
    for pos in range(1, 9):
        pmids = json.dumps([100 + pos] * pos)  # mas PMIDs cuanto mayor la posicion
        rows.append(("P00000", pos, "phosphorylation", "S", "A", pmids))
    _build_fixture_db(db_path, rows)

    sites = load_candidate_sites("P00000", db_path)
    assert len(sites) == MAX_SITES_PER_ACCESSION_TYPE
    # se quedan los de mas PMIDs (posiciones altas), no los primeros por posicion
    kept_positions = sorted(s[1] for s in sites)
    assert kept_positions == [4, 5, 6, 7, 8]


def test_load_candidate_sites_prefiere_tier_a_sobre_tier_b_en_el_tope(tmp_path):
    db_path = tmp_path / "lookup.sqlite3"
    rows = [
        # 1 sitio tier A con pocos PMIDs
        ("P00000", 1, "acetylation", "K", "A", json.dumps([111, 222])),
        # 5 sitios tier B con muchos PMIDs cada uno -- no deberian desplazar al tier A
        *[("P00000", 10 + i, "acetylation", "K", "B", json.dumps([900 + i] * 20)) for i in range(5)],
    ]
    _build_fixture_db(db_path, rows)

    sites = load_candidate_sites("P00000", db_path)
    assert len(sites) == MAX_SITES_PER_ACCESSION_TYPE
    assert ("acetylation", 1, "K", "A", (111, 222)) in sites


def test_emit_python_incluye_accession_longitud_y_sitios():
    sites = [("phosphorylation", 15, "S", "A", (111, 222))]
    draft = emit_python("P04637", "P53_HUMAN", 393, sites)
    assert 'uniprot_accession="P04637"' in draft
    assert "length=393" in draft
    assert 'GroundTruthSite(15, "S", "phosphorylation", "A", (111, 222))' in draft


def test_emit_python_pmid_unico_lleva_coma_final_de_tupla():
    sites = [("phosphorylation", 15, "S", "B", (111,))]
    draft = emit_python("P04637", "P53_HUMAN", 393, sites)
    assert "(111,)" in draft
