#!/usr/bin/env python
"""Reconstruye summary.tsv a partir de los *_calibration.tsv ya generados por
generate_deepptmpred_calibration.py -- ese script solo escribe summary.tsv al
final de SU propia corrida (todos los tipos pasados en --ptm-types), asi que
tras dos corridas parciales (11 tipos + 6 tipos) el summary quedaba con solo
los ultimos 6. Este script agrega los 17 sin re-ejecutar el modelo.

## Regla de umbral: FPR-constrained TPR-max, no Youden's J

Youden's J (maximizar tpr-fpr) degenera cuando el AUROC de una muestra chica
es malo: el punto que maximiza J puede ser el "rechazar todo" sintetico que
sklearn.roc_curve antepone (threshold = max(score)+1), lo que dio
umbral=inf en o_linked_glycosylation en la primera version de este script.

En vez de inventar una regla propia, se porto la que el propio DeepPTMPred
implementa -- verificado leyendo su codigo real
(``DeepPTMPred/pred/train_PTM/trainer.py``, funciones
``optimize_threshold_high_specificity``/``_high_specificity2``/``_f1``,
lineas ~130-173) y confirmado empiricamente contra los resultados que SI
publico (``pred/train_PTM/result/results_*/``, columna FP Rate):

- El propio trainer.py tiene un bug de despacho real (L356-364, ``if`` en
  vez de ``elif``): la rama de fosforilacion (FPR<=0.18) queda pisada por el
  ``else`` (FPR<=0.20), y solo n_linked_glycosylation llega a la rama F1-max
  -- que en los resultados publicados le da FP Rate=0.94 (F1-max sobre datos
  desbalanceados es igual de fragil que Youden's J, mismo tipo de fallo).
  Ese bug esta en un script de ENTRENAMIENTO que este proyecto nunca ejecuta
  (solo se usa predict.py para inferencia), asi que no se toca el vendored
  trainer.py -- pero tampoco se replica el bug aqui.
- Se estandariza en FPR-constrained TPR-max para los 17 tipos (evita F1-max
  por completo, dado el caso confirmado de degeneracion), con el bound mas
  estricto (FPR<=0.18, la funcion ``optimize_threshold_high_specificity``
  real, nunca alcanzada por el bug pero si valida) reservado para
  phosphorylation, tal como describe el paper (Briefings in Bioinformatics
  27(3) bbag321) como el tipo que exige alta especificidad; el resto usa
  FPR<=0.20 (``optimize_threshold_high_specificity2``, la que de hecho
  produjo las metricas publicadas para 16/17 tipos por el bug de despacho).
  Fallback a 0.5 si ningun umbral candidato cumple la cota, igual que el
  codigo original.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

CALIB_DIR = Path(__file__).resolve().parent.parent / "DeepPTMPred" / "data" / "calibration"

# FPR maximo tolerado por tipo -- ver docstring del modulo.
HIGH_SPECIFICITY_MAX_FPR = 0.18
DEFAULT_MAX_FPR = 0.20
HIGH_SPECIFICITY_TYPES = {"phosphorylation"}


def _fpr_constrained_threshold(labels: np.ndarray, probs: np.ndarray, max_fpr: float) -> float:
    """Port fiel de optimize_threshold_high_specificity(2) de trainer.py: entre los
    umbrales con fpr <= max_fpr, el que maximiza tpr; fallback 0.5 si ninguno cumple.

    Extra respecto al original (necesario en muestra chica/ruidosa, no en el
    codigo de trainer.py porque nunca se corrio sobre datos tan pequenos):
    sklearn.roc_curve antepone un punto sintetico "rechazar todo"
    (fpr=0, tpr=0, threshold=max(score)+1) que satisface cualquier max_fpr
    trivialmente. thresholds[0] nunca es un umbral real (max(score)+1 > 1,
    predict.py jamas podria alcanzarlo), asi que se excluye antes de buscar
    el maximo -- si no, un desempate a tpr=0 lo elige y da threshold=inf.
    """
    fpr, tpr, thresholds = roc_curve(labels, probs)
    fpr, tpr, thresholds = fpr[1:], tpr[1:], thresholds[1:]
    valid = fpr <= max_fpr
    if not valid.any():
        return 0.5
    return float(thresholds[valid][np.argmax(tpr[valid])])


def main() -> None:
    rows = []
    for f in sorted(CALIB_DIR.glob("*_calibration.tsv")):
        ptm_type = f.stem.replace("_calibration", "")
        df = pd.read_csv(f, sep="\t")
        n_pos = int((df["label"] == 1).sum())
        n_neg = int((df["label"] == 0).sum())
        if n_pos == 0 or n_neg == 0:
            print(f"{ptm_type}: solo una clase presente ({n_pos} pos / {n_neg} neg), se omite del resumen.")
            continue
        auroc = roc_auc_score(df["label"], df["probability"])
        max_fpr = HIGH_SPECIFICITY_MAX_FPR if ptm_type in HIGH_SPECIFICITY_TYPES else DEFAULT_MAX_FPR
        threshold = _fpr_constrained_threshold(df["label"].to_numpy(), df["probability"].to_numpy(), max_fpr)
        rows.append({
            "ptm_type": ptm_type, "n_sites": len(df), "n_pos": n_pos, "n_neg": n_neg,
            "n_proteins": df["protein_id"].nunique(), "auroc": auroc,
            "max_fpr": max_fpr, "suggested_threshold": threshold,
        })

    summary = pd.DataFrame(rows).sort_values("ptm_type").reset_index(drop=True)
    out_path = CALIB_DIR / "summary.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"\nResumen escrito en {out_path} ({len(summary)} tipos):")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
