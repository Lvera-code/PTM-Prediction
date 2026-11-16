#!/usr/bin/env python
"""Standalone runner for Kinase Library (informative corroboration of kinase specificity).

Biological coherence analysis 2026-08-07, point 5: no engine of this
pipeline (DeepMVP/DeepPTMPred) distinguishes WHICH kinase family phosphorylates a
site; they all predict "phosphorylatable in general", never by whom. Unlike
the polyubiquitin chain type (a later cellular event not
observable from sequence, a line of reasoning that lived in Phase A/3c,
removed from scope 2026-08-10), kinase specificity IS a
local sequence property around the site, and it was verified on 2026-08-07
that a real, published and downloadable source exists: Johnson et al. 2023
Nature ("An atlas of substrate specificities for the human serine/threonine
kinome", 303 kinases) + Yaron-Barir et al. 2024 Nature ("The intrinsic
substrate specificity of the human tyrosine kinome"), packaged together in
``kinase-library`` (PyPI, github.com/TheKinaseLibrary/kinase-library, licence
CC-BY-NC-SA-3.0, the same non-commercial family already accepted for EMNGly).

``kl.Substrate(sequence, phos_pos=position)`` accepts the COMPLETE sequence of
the protein plus the 1-based position of the phosphoacceptor; the library itself
extracts and pads (with '_', treated as a wildcard) the 15-residue window
centred on the site, so there is no need to reimplement that trimming here.
It auto-detects Ser/Thr vs Tyr from the real residue at that position (verified
2026-08-07 with a real example of each type: p53 S33 -> top hit ATM, which matches
the real literature on the DNA damage response; a synthetic Y
phosphosite -> top hits from the Tyr kinome, never from the Ser/Thr kinome).

A DEDICATED environment (the ``kinase_library`` conda env, see ``Settings.KINASE_LIBRARY_*``):
the ``numpy~=1.26.4``/``pandas~=2.2.3`` that the package itself pins are
incompatible with the versions pinned in this pipeline's main venv
(``numpy==2.2.6``/``pandas==2.3.3``, see requirements.txt), the same reason
why DeepMVP/DeepPTMPred/MeToken do not live in the main venv either.

An individual site that fails (the real residue is not S/T/Y at that position,
the position is out of range, etc.) NEVER brings down the rest of the batch; it is skipped and
execution continues, the same criterion as ``_metoken_runner.py``.
"""

import argparse
from pathlib import Path

import pandas as pd

OUTPUT_COLUMNS = [
    "position", "kinase_library_top_kinase", "kinase_library_top_family",
    "kinase_library_percentile", "kinase_library_top3_kinases",
]


def _score_position(kl_module, sequence: str, position: int) -> dict:
    """Puntua UNA posicion y devuelve el top kinasa/familia/percentil/top-3.

    ``phos_pos`` es 1-based, mismo convenio que el resto del pipeline --
    pasado directo a ``kl.Substrate``, que internamente hace su propio
    recorte+relleno de la ventana de 15-mer (ver docstring del modulo).
    """
    substrate = kl_module.Substrate(sequence, phos_pos=position)
    predictions = substrate.predict()

    ranked = predictions.sort_values("Percentile Rank")
    top_kinase = ranked.index[0]
    top_percentile = float(ranked.loc[top_kinase, "Percentile"])
    top3_kinases = ",".join(ranked.index[:3])

    top_family = kl_module.get_kinase(top_kinase).family

    return {
        "kinase_library_top_kinase": top_kinase,
        "kinase_library_top_family": top_family,
        "kinase_library_percentile": top_percentile,
        "kinase_library_top3_kinases": top3_kinases,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Runner standalone de Kinase Library (corroboracion informativa de especificidad de quinasa)."
    )
    parser.add_argument("--sequence", required=True, help="Secuencia COMPLETA de la proteina.")
    parser.add_argument(
        "--positions", required=True, type=int, nargs="+",
        help="Posiciones 1-based de fosforilacion ya aceptadas por el consenso (pasa_umbral=true) a corroborar.",
    )
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    import kinase_library as kl

    rows = []
    for position in args.positions:
        try:
            scored = _score_position(kl, args.sequence, position)
        except Exception as exc:  # noqa: BLE001 -- un sitio individual no debe tumbar el lote completo
            print(f"[kinase_library_runner] posicion {position} omitida: {type(exc).__name__}: {exc}")
            continue
        rows.append({"position": position, **scored})

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(out_csv, index=False)

    print(f"Kinase Library: {len(rows)}/{len(args.positions)} sitio(s) puntuados -> '{out_csv}'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
