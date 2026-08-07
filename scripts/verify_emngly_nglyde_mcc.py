"""Go/no-go check 2 de EMNGly (ver STATUS.md 'Decision 2'): reproduce MCC~=0.736
corriendo el pipeline REAL (ESM-1b + MIF + SVM, exactamente los mismos componentes
que ``src/engines/_emngly_runner.py``) sobre el set independiente de N-GlyDE
(``EMNgly/data/N-GlyDE/NGLYDE_independent.txt``, 447 filas / 86 proteinas, via
``dukkakc/DeepNGlyPred`` -- el propio repo de EMNgly solo documenta la fuente, no
la empaqueta).

Requiere que ``scripts/prepare_emngly_nglyde_structures.py`` ya haya corrido (Fase
1.5 real sobre un modelo AlphaFold por proteina, ver su docstring para el porque
AlphaFold y no PDB cristalografico).

Bug real confirmado leyendo ``EMNgly/predict.py::get_scores`` (no arreglado en el
vendorizado, ver docstring de modulo alli): el script original llama
``get_scores(label_y, predict_y[0])`` -- ``predict_y[0]`` es un escalar (la
probabilidad de la PRIMERA fila), no la lista completa de probabilidades. Este
script reimplementa ``get_scores`` correctamente (todas las probabilidades, no
solo la primera) en vez de tocar el vendorizado.

Corre en ``.venv-emngly`` (fair-esm/torch/sklearn) -- NUNCA en ``cnb_pipeline``.
"""

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    matthews_corrcoef,
    recall_score,
    roc_auc_score,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src" / "engines"))

import _emngly_runner as runner_lib  # noqa: E402

DATASET_TSV = REPO_ROOT / "EMNgly" / "data" / "N-GlyDE" / "NGLYDE_independent.txt"
STRUCTURES_DIR = REPO_ROOT / "EMNgly" / "data" / "N-GlyDE" / "structures"
MIF_WEIGHTS = REPO_ROOT / "EMNgly" / "model" / "MIF" / "weights" / "mif.pt"
ESM_CHECKPOINT = REPO_ROOT / "EMNgly" / "esm" / "checkpoints" / "esm1b_t33_650M_UR50S.pt"
SVM_CHECKPOINT = REPO_ROOT / "EMNgly" / "checkpoints" / "N-GlyDE.pickle"
CACHE_DIR = REPO_ROOT / "EMNgly" / "cache" / "nglyde_verification"
PUBLISHED_MCC = 0.736


def get_scores(label_y, predict_y, th=0.5) -> dict:
    """Puerto CORREGIDO de ``EMNgly/predict.py::get_scores`` -- recibe la lista
    COMPLETA de probabilidades, no ``predict_y[0]`` (bug real confirmado en el
    vendorizado, ver docstring del modulo)."""
    predict_y = np.asarray(predict_y)
    label_y = np.asarray(label_y)

    auc = roc_auc_score(label_y, predict_y)
    predict_label = np.where(predict_y > th, 1, 0)
    mcc = matthews_corrcoef(label_y, predict_label)
    acc = accuracy_score(label_y, predict_label)
    tn, fp, fn, tp = confusion_matrix(label_y, predict_label).ravel()
    specificity = tn / (tn + fp)
    sensitivity = recall_score(label_y, predict_label)

    return {
        "auc": round(float(auc), 5),
        "mcc": round(float(mcc), 5),
        "acc": round(float(acc), 5),
        "specificity": round(float(specificity), 5),
        "sensitivity": round(float(sensitivity), 5),
        "n": int(len(label_y)),
        "n_positive": int(label_y.sum()),
    }


def main() -> int:
    import torch

    df = pd.read_csv(DATASET_TSV, sep="\t")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Mismo sys.path que _emngly_runner.py::main() -- 'MIF' se importa como
    # paquete top-level (necesita 'EMNgly/model') y sus propios submodulos
    # hacen imports BARE de 'sequence_models.xxx' (necesita
    # 'EMNgly/model/MIF' tambien, ver docstring del runner). Como este script
    # importa las funciones del runner directamente (no via su main()), este
    # paso hay que repetirlo aqui.
    emngly_home = REPO_ROOT / "EMNgly"
    sys.path.insert(0, str(emngly_home / "model"))
    sys.path.insert(0, str(emngly_home / "model" / "MIF"))

    print("[verify] cargando ESM-1b (puede tardar)...")
    extractor = runner_lib._ESMEmbeddingExtractor(ESM_CHECKPOINT, device)
    print("[verify] cargando SVM...")
    with open(SVM_CHECKPOINT, "rb") as f:
        svm = pickle.load(f)

    labels, probs, skipped = [], [], []
    accessions = sorted(df["UniProt_ID"].unique())
    for i, accession in enumerate(accessions, start=1):
        # prepare_emngly_nglyde_structures.py descarga a '{accession}_af.pdb'
        # -- parse_structure deriva el accession interno del STEM del archivo
        # de entrada, asi que los artefactos quedan con sufijo '_af', no el
        # accession UniProt puro (confirmado en el log real de esa corrida).
        struct_dir = STRUCTURES_DIR / accession
        derived_fasta = struct_dir / f"{accession}_af_derived.fasta"
        mapping_csv = struct_dir / f"{accession}_af_position_mapping.csv"
        chain_pdbs = list(struct_dir.glob(f"{accession}_af_chain_*.pdb"))
        rows = df[df["UniProt_ID"] == accession]

        if not derived_fasta.is_file() or not mapping_csv.is_file() or not chain_pdbs:
            print(f"[verify] {accession}: sin estructura preparada, se omiten {len(rows)} fila(s).")
            skipped.extend(rows.index.tolist())
            continue

        sequence = "".join(
            l.strip() for l in derived_fasta.read_text().splitlines() if not l.startswith(">")
        )
        mapping = pd.read_csv(mapping_csv)
        fasta_by_pdb_seqid = dict(zip(mapping["pdb_seqid"].astype(int), mapping["fasta_position"].astype(int)))

        try:
            structure_emb = runner_lib._load_or_compute_structure_emb(
                chain_pdbs[0], MIF_WEIGHTS, CACHE_DIR, accession
            )
            site_full = runner_lib._load_or_compute_esm_full(sequence, extractor, CACHE_DIR, accession)
        except Exception as exc:  # noqa: BLE001 -- fallo real de una proteina no debe tumbar el resto
            print(f"[verify] {accession}: fallo computando embeddings ({exc}), se omiten {len(rows)} fila(s).")
            skipped.extend(rows.index.tolist())
            continue

        for idx, row in rows.iterrows():
            asn_pos = int(row["ASN_Pos"])  # numeracion UniProt == pdb_seqid (AlphaFold, sin huecos)
            fasta_pos = fasta_by_pdb_seqid.get(asn_pos)
            structure_idx = asn_pos - 1

            if (
                fasta_pos is None
                or fasta_pos >= site_full.shape[0]
                or structure_idx < 0
                or structure_idx >= structure_emb.shape[0]
                or sequence[fasta_pos - 1] != "N"
            ):
                skipped.append(idx)
                continue

            site_vec = site_full[fasta_pos]
            left = max(0, fasta_pos - runner_lib._LOCAL_WINDOW_LEFT)
            right = min(fasta_pos + runner_lib._LOCAL_WINDOW_RIGHT, len(sequence))
            local_vec = extractor.extract([sequence[left:right]]).cpu().numpy()[0]
            struct_vec = structure_emb[structure_idx]

            feature_row = np.concatenate([site_vec, local_vec, struct_vec]).reshape(1, -1)
            probability = float(svm.predict_proba(feature_row)[0, 1])

            labels.append(int(row["Evidence"]))
            probs.append(probability)

        print(f"[verify] {i}/{len(accessions)} proteina(s) procesada(s), {len(labels)} sitio(s) acumulado(s)...")

    print(f"[verify] {len(labels)} sitio(s) evaluado(s), {len(skipped)} fila(s) omitida(s) de {len(df)} totales.")
    if not labels:
        print("[verify] sin datos evaluables, no se puede calcular MCC.")
        return 1

    scores = get_scores(labels, probs)
    print(f"[verify] scores: {scores}")
    print(f"[verify] MCC publicado (paper, N-GlyDE test set): {PUBLISHED_MCC}")
    print(f"[verify] MCC reproducido en esta maquina: {scores['mcc']}")

    out_csv = CACHE_DIR / "nglyde_predictions.csv"
    pd.DataFrame({"label": labels, "probability": probs}).to_csv(out_csv, index=False)
    print(f"[verify] predicciones crudas guardadas en '{out_csv}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
