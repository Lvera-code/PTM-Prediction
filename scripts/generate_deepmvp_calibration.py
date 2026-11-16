#!/usr/bin/env python
"""Generates the calibration `site_prediction.tsv` files missing from `DeepMVP/models/<tipo>/`.

Standalone: it requires the DeepMVP stack (TensorFlow 2.4.2, Python 3.7), so
run it with the `deepmvp` conda env (``conda run -n deepmvp python
scripts/generate_deepmvp_calibration.py``), never with the pipeline's
main interpreter.

## Why it exists

`DeepMVPEngine`/the Phase 3 core (``ptm_annotation.py``) filter by the
``fpr`` column, computed by ``DeepMVP/lib/Metrics.py::add_confidence_metrics``
reading ``DeepMVP/models/<tipo>/site_prediction.tsv`` (a VALIDATION file with
a real ``y`` column). That file is NOT included in
``models.tar.gz`` (confirmed real 2026-07-28), and is confirmed
also as a known upstream bug via issues #1/#2 of
`github.com/bzhanglab/DeepMVP` (two other independent users reported
the same thing in 2025-09 and 2026-06; the PR that closed issue #1 did not fix
this).

## How it is generated without needing an external FASTA

`https://deepmvp.ptmax.org/all_data.tar.gz` (55MB, not listed in the Shiny
app's static HTML but accessible by direct URL, confirmed real)
contains the labelled train/test dataset (columns ``protein, aa, pos,
x, y``) for DeepMVP's 8 types. The ``x`` column already comes with the
MAXIMUM width that any submodel of the ensemble needs (61 = flank 30,
confirmed by reading each ``model.json``: no submodel of the 8 types
asks for ``peptide_length`` > 61), and the submodels that ask for a narrower
window are derived by trimming that same centred window (the same criterion
as ``DeepMVP/lib/DataIO.py::getPeptideSequence``), without needing the
protein's complete sequence or a ``db`` FASTA (unlike
``DeepMVP.py predict -i <test_file>``, which DOES require it unconditionally
via ``processing_prediction_data``, another confirmed real bug). It reuses
DeepMVP's real code directly
(``lib.PeptideEncode.encodePeptides``, ``lib.Utils.combine_rts``) so that
the result is identical to what ``DeepMVP.py`` itself would produce.

Verified 2026-07-28 with all 8 types: AUROC 0.90-0.99 (matching the
figures published in the paper), which confirms that the generation is correct, not
merely that it runs without error.
"""

import argparse
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

DEEPMVP_ALL_DATA_URL = "https://deepmvp.ptmax.org/all_data.tar.gz"

TYPE_TO_PREFIX = {
    "acetylation_k": "acet_k",
    "glycosylation_n": "gly_n",
    "methylation_k": "met_k",
    "methylation_r": "met_r",
    "phosphorylation_st": "phos_st",
    "phosphorylation_y": "phos_y",
    "sumoylation_k": "sumo_k",
    "ubiquitination_k": "ubi_k",
}

MAX_WIDTH = 61  # ancho real de 'x' en all_data.tar.gz (flank=30), verificado para los 8 tipos


def download_all_data(dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    archive_path = dest_dir / "all_data.tar.gz"
    if not archive_path.is_file():
        print(f"Descargando {DEEPMVP_ALL_DATA_URL} -> {archive_path}")
        urllib.request.urlretrieve(DEEPMVP_ALL_DATA_URL, archive_path)
    extracted = dest_dir / "all_data"
    if not extracted.is_dir():
        with tarfile.open(archive_path) as tf:
            tf.extractall(dest_dir)
    return extracted


def trim_window(x: str, peptide_length: int) -> str:
    trim = (MAX_WIDTH - peptide_length) // 2
    return x if trim == 0 else x[trim: MAX_WIDTH - trim]


def generate_for_type(ptm_type: str, prefix: str, models_dir: Path, all_data_dir: Path, testing_suffix: str) -> None:
    from lib.PeptideEncode import encodePeptides
    from lib.Utils import combine_rts
    from tensorflow.keras.models import load_model

    model_dir = models_dir / ptm_type
    with open(model_dir / "model.json") as f:
        model_list = json.load(f)

    test_file = all_data_dir / f"{prefix}_testing_{testing_suffix}.tsv"
    df = pd.read_csv(test_file, sep="\t")
    if df["x"].str.len().nunique() != 1 or df["x"].str.len().iloc[0] != MAX_WIDTH:
        raise ValueError(f"{ptm_type}: ancho de 'x' inesperado en '{test_file}', revisar antes de continuar")

    per_model_preds = []
    for name, dp_model in model_list["dp_model"].items():
        peptide_length = dp_model["peptide_length"]
        windows = df["x"].apply(lambda s: trim_window(s, peptide_length))
        x_enc = encodePeptides(pd.DataFrame({"x": windows}))
        model_path = model_dir / Path(dp_model["model"]).name
        model = load_model(str(model_path))
        y_prob = model.predict(x_enc, batch_size=2048, verbose=0)
        per_model_preds.append(np.asarray(y_prob).reshape(-1))
        print(f"  [{ptm_type}] modelo {name} (peptide_length={peptide_length}) listo")

    res = np.stack(per_model_preds, axis=1)
    y_pred = np.apply_along_axis(combine_rts, 1, res, method="mean", remove_outlier=True)

    out = df.copy()
    out["y_pred"] = y_pred
    out_path = model_dir / "site_prediction.tsv"
    if out_path.is_file():
        print(
            f"  AVISO: '{out_path}' ya existe (generado con un --testing-suffix "
            f"previo, posiblemente distinto a '{testing_suffix}') -- se sobreescribe."
        )
    out.to_csv(out_path, sep="\t", index=False)
    print(f"'{ptm_type}' -> '{out_path}' ({len(out)} filas, testing_suffix='{testing_suffix}')")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deepmvp-home", default="DeepMVP", help="Raiz del repo DeepMVP clonado (contiene models/).")
    parser.add_argument("--work-dir", default="DeepMVP/.calibration_data", help="Donde descargar/extraer all_data.tar.gz.")
    parser.add_argument("--testing-suffix", default="70", choices=["70", "80", "90"],
                        help="Subconjunto de test a usar. Confirmado leyendo el Methods real del paper "
                             "(Nature Methods 2025, PMC12446062): los tres son EL MISMO 10%% de test "
                             "genuinamente separado del 90%% train+validation con el que se entrenaron "
                             "los pesos shipeados, filtrado ADEMAS por CD-HIT contra train+validation al "
                             "umbral de identidad indicado -- 70 = filtro mas estricto (excluye cualquier "
                             "peptido con >=70%% identidad a train/val), el mas conservador de los tres.")
    parser.add_argument("--types", nargs="*", default=list(TYPE_TO_PREFIX), choices=list(TYPE_TO_PREFIX))
    args = parser.parse_args()

    deepmvp_home = Path(args.deepmvp_home).resolve()
    sys.path.insert(0, str(deepmvp_home))

    all_data_dir = download_all_data(Path(args.work_dir).resolve())
    models_dir = deepmvp_home / "models"

    for ptm_type in args.types:
        print(f"=== {ptm_type} ===")
        generate_for_type(ptm_type, TYPE_TO_PREFIX[ptm_type], models_dir, all_data_dir, args.testing_suffix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
