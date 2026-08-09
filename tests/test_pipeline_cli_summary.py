"""Tests de los helpers puros de pipeline.py que arman la tabla '-- Consenso real --'
del CLI (mejoras 2026-08-09/2026-08-10: fuente de literatura, PMIDs, tipo de
MeToken en desacuerdo, familia de quinasa, conteo de crosstalk, motor abreviado)."""

from pipeline import _abbrev_motor, _crosstalk_count, _fmt_pmids


def test_fmt_pmids_default_muestra_solo_el_primero():
    # Default max_shown=1 desde 2026-08-10: el primer PMID es (~92% del panel)
    # el mas viejo/original, no necesariamente "el mas importante" -- ver
    # docstring de _fmt_pmids.
    assert _fmt_pmids((111, 222)) == "111+1"


def test_fmt_pmids_un_solo_pmid_sin_sufijo():
    assert _fmt_pmids((111,)) == "111"


def test_fmt_pmids_respeta_max_shown_explicito():
    assert _fmt_pmids((111, 222, 333, 444), max_shown=2) == "111,222+2"


def test_crosstalk_count_cuenta_los_tipos():
    aviso = (
        "Compite con: acetylation, ubiquitination (mismo residuo, mutuamente "
        "excluyentes en una misma molecula/instante -- ver docstring del modulo)"
    )
    assert _crosstalk_count(aviso) == 2


def test_crosstalk_count_no_nombra_tipos_solo_cuenta():
    # Caso real encontrado 2026-08-09: K357 de p53 compite con 6 tipos. El orden
    # de dbPTM/_add_ptm_crosstalk_flag es alfabetico, no por relevancia -- por
    # eso esta funcion deliberadamente solo cuenta (decision 2026-08-10).
    aviso = (
        "Compite con: acetylation, crotonylation, glutarylation, lys_methylation, "
        "malonylation, sumoylation (mismo residuo, mutuamente excluyentes en una "
        "misma molecula/instante -- ver docstring del modulo)"
    )
    assert _crosstalk_count(aviso) == 6


def test_crosstalk_count_cero_si_no_matchea_formato():
    assert _crosstalk_count("texto sin el marcador esperado") == 0


def test_abbrev_motor_consenso():
    assert _abbrev_motor("DeepMVP+DeepPTMPred") == "DMVP+DPTMP"


def test_abbrev_motor_un_solo_motor():
    assert _abbrev_motor("DeepMVP") == "DMVP"
