#!/usr/bin/env python
"""Runner standalone de Kinase Library (corroboracion informativa de especificidad de quinasa).

Analisis de coherencia biologica 2026-08-07, punto 5: ningun motor de este
pipeline (DeepMVP/DeepPTMPred) distingue QUE familia de quinasa fosforila un
sitio -- todos predicen "fosforilable en general", nunca por quien. A
diferencia del tipo de cadena de poliubiquitina (evento celular posterior no
observable a partir de secuencia -- razonamiento que vivia en Fase A/3c,
eliminada del alcance 2026-08-10), la especificidad de quinasa SI es una
propiedad local de secuencia alrededor del sitio -- verificado 2026-08-07
que existe una fuente real, publicada y descargable: Johnson et al. 2023
Nature ("An atlas of substrate specificities for the human serine/threonine
kinome", 303 quinasas) + Yaron-Barir et al. 2024 Nature ("The intrinsic
substrate specificity of the human tyrosine kinome"), empaquetadas juntas en
``kinase-library`` (PyPI, github.com/TheKinaseLibrary/kinase-library, licencia
CC-BY-NC-SA-3.0 -- misma familia no comercial ya aceptada para EMNGly).

``kl.Substrate(sequence, phos_pos=position)`` acepta la secuencia COMPLETA de
la proteina + la posicion 1-based del fosfoaceptor -- la libreria misma
extrae y rellena (con '_', tratado como comodin) la ventana de 15 residuos
centrada en el sitio, sin necesidad de reimplementar ese recorte aqui.
Auto-detecta Ser/Thr vs Tyr por el residuo real en esa posicion (verificado
2026-08-07 con un ejemplo real de cada tipo: p53 S33 -> top hit ATM, coincide
con la literatura real de respuesta a dano en el ADN; un fosfositio Y
sintetico -> top hits del kinoma Tyr, nunca del kinoma Ser/Thr).

Entorno DEDICADO (conda env ``kinase_library``, ver ``Settings.KINASE_LIBRARY_*``):
``numpy~=1.26.4``/``pandas~=2.2.3`` que fija el propio paquete son
incompatibles con las versiones fijadas del venv principal de este pipeline
(``numpy==2.2.6``/``pandas==2.3.3``, ver requirements.txt) -- mismo motivo por
el que DeepMVP/DeepPTMPred/MeToken/StackGlyEmbed tampoco viven en el venv
principal.

Un sitio individual que falla (residuo real no es S/T/Y en esa posicion,
posicion fuera de rango, etc.) NUNCA tumba el resto del lote -- se omite y se
continua, mismo criterio que ``_metoken_runner.py``.
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
