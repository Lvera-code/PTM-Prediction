#!/usr/bin/env python
"""Runner standalone para StackGlyEmbed (corroboracion informativa de N-glicosilacion).

NUNCA se importa desde el paquete ``src`` -- requiere torch/tensorflow/
transformers/scikit-learn/ProteinBERT, dependencias SOLO presentes en el
venv dedicado del proyecto HERMANO ``B-Cell-Epitope-Prediction``
(``Settings.STACKGLYEMBED_PYTHON_BIN``, ver ``src/config/settings.py``). Se
invoca EXCLUSIVAMENTE via subprocess desde
``src/engines/stackglyembed_engine.py``, mismo patron que
``_deepptmpred_runner.py``/``_metoken_runner.py``.

## Por que existe (rol en el pipeline, decision 2026-08-01)

``n_linked_glycosylation`` en DeepPTMPred esta CONFIRMADO como modelo muerto
(AUROC ~0.51, ver STATUS.md, ya excluido del consenso en
``ptm_annotation.py::CONSENSUS_EXCLUDED_TYPES``) -- no arreglable
reentrenando, un problema real del dataset/modelo publicado (SMOTE no-op +
85.6% de positivos en el training set). StackGlyEmbed
(``github.com/GaryChan-lab/StackGlyEmbed``) es un tercer motor,
INDEPENDIENTE de arquitectura (ProteinBERT + ESM-2 650M + ProtT5 apilados ->
meta-clasificador SVM), especializado UNICAMENTE en N-glicosilacion -- ya
instalado y verificado funcionando de verdad en el proyecto hermano
``B-Cell-Epitope-Prediction`` (proyecto 1, independiente de este por
decision explicita 2026-07-26).

Este runner es una adaptacion propia (NO una copia importada -- viola la
independencia entre proyectos) de la logica real ya verificada en
``B-Cell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py``
(extraccion de features + prediccion, formato ``dataset.txt`` original del
repo). Reusa el MISMO venv/pickles ya instalados alli como recurso externo
(``Settings.STACKGLYEMBED_PYTHON_BIN``/``STACKGLYEMBED_MODELS_DIR``, mismo
criterio que cualquier otro motor externo de este proyecto) -- nunca
reinstala ESM-2/ProtT5/ProteinBERT de nuevo.

Simplificacion respecto al script del proyecto hermano: aqui SIEMPRE se
evalua una UNICA proteina por invocacion (mismo patron que
``_deepptmpred_runner.py``/``_metoken_runner.py``, un accession a la vez),
asi que no hace falta el formato multi-proteina ``dataset.txt`` -- este
runner recibe ``--sequence``/``--positions`` directamente por CLI y escribe
el CSV de salida sin pasar por un archivo intermedio.

## Convencion de posiciones (verificado en el proyecto hermano, no reverificado aqui)

El script REAL soporta proteinas completas sin truncar (ESM-2 chunkea por
bloques de 1024aa, ProtT5 por 8797aa, ver ``_get_esm2_embedding``/
``_get_prott5_embedding`` abajo) -- si se le pasa la SECUENCIA COMPLETA de
la proteina (no un fragmento/peptido), la posicion 1-indexada del secuon
coincide EXACTAMENTE con la numeracion global que ya usa este proyecto
(misma convencion que ``sequence``/``posicion`` en
``annotate_fasta_path``/``annotate_pdb_path``) -- sin necesidad de convertir
offsets de ventana. Quien llama (``stackglyembed_engine.py``) siempre pasa
``sequence`` completo (Fase 1/1.5), nunca un peptido recortado.

## Formato de salida

CSV con columnas ``position`` (1-based, la Asparagina del secuon),
``stackglyembed_veredicto`` (``'Glicosilado'``/``'No glicosilado'``) y
``stackglyembed_score`` (probabilidad cruda del meta-clasificador SVM).
"""

import argparse
import os
import pickle
import re
import sys
from pathlib import Path

os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUTPUT_COLUMNS = ["position", "stackglyembed_veredicto", "stackglyembed_score"]

# Ventana de residuos alrededor del sitio para el promedio de ESM-2 (mismo
# valor que ``B-Cell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py``,
# heredado del script original del repo -- los clasificadores de
# ``base_layer_pickle_files/`` se entrenaron contra esta ventana exacta).
_WINDOW_SIZE = 15
# Tamanos de chunk de cada embedder (identicos al script del proyecto
# hermano, ver docstring arriba) -- soportan proteinas completas sin
# truncar.
_ESM2_CHUNK_SIZE = 1024
_PROTT5_CHUNK_SIZE = 8797


def _get_model_with_global_embedding_as_outputs(model):
    """Reconstruye ProteinBERT para exponer la embedding global (identico al proyecto hermano)."""
    from tensorflow import keras

    global_layers = [
        layer.output
        for layer in model.layers
        if len(layer.output.shape) == 2 and layer.name in ["global-merge2-norm-block6"]
    ]
    concatenated = keras.layers.Concatenate(name="last-Window-layers")(global_layers)
    return keras.models.Model(inputs=model.inputs, outputs=concatenated)


def _get_proteinbert_representation(pretrained_model_generator, input_encoder, seq: str) -> np.ndarray:
    encoded_x = input_encoder.encode_X([seq], len(seq) + 2)
    model = _get_model_with_global_embedding_as_outputs(pretrained_model_generator.create_model(len(seq) + 2))
    return np.array(model.predict(encoded_x, batch_size=2))[0]


def _get_esm2_embedding(tokenizer, model, seq: str) -> np.ndarray:
    """Embedding ESM-2 por residuo (ultima capa, sin CLS/EOS), chunkeado para proteinas completas."""
    import torch

    chunks = [seq[i : i + _ESM2_CHUNK_SIZE] for i in range(0, len(seq), _ESM2_CHUNK_SIZE)]
    final = np.zeros((1, model.config.hidden_size))
    for chunk in chunks:
        tokens = tokenizer(chunk, return_tensors="pt")
        with torch.no_grad():
            out = model(**tokens)
        rep = out.last_hidden_state[0, 1:-1].numpy()
        final = np.concatenate((final, rep), axis=0)
    return np.delete(final, 0, axis=0)


def _get_prott5_embedding(tokenizer, model, seq: str) -> np.ndarray:
    """Embedding ProtT5 por residuo, chunkeado para proteinas completas."""
    import torch

    chunks = [seq[i : i + _PROTT5_CHUNK_SIZE] for i in range(0, len(seq), _PROTT5_CHUNK_SIZE)]
    final = np.zeros((1, model.config.d_model))
    for chunk in chunks:
        spaced = " ".join(list(re.sub(r"[UZOB]", "X", chunk)))
        ids = tokenizer([spaced], add_special_tokens=True, padding="longest", return_tensors="pt")
        with torch.no_grad():
            out = model(input_ids=ids["input_ids"], attention_mask=ids["attention_mask"])
        emb = out.last_hidden_state[0, : len(chunk)].numpy()
        final = np.concatenate((final, emb), axis=0)
    return np.delete(final, 0, axis=0)


def _extract_features(sequence: str, positions, t5_model_path: str, esm_model_name: str) -> np.ndarray:
    """Genera la matriz de features (ProteinBERT global + ESM-2 ventana + ProtT5 puntual), una fila por posicion."""
    from proteinbert import load_pretrained_model
    from transformers import AutoTokenizer, EsmModel, T5EncoderModel, T5Tokenizer

    print("Cargando ProteinBERT (local, ~/proteinbert_models/default.pkl)...", file=sys.stderr, flush=True)
    pretrained_model_generator, input_encoder = load_pretrained_model(download_model_dump_if_not_exists=False)

    print(f"Cargando ESM-2 650M ({esm_model_name}, offline local)...", file=sys.stderr, flush=True)
    esm_tokenizer = AutoTokenizer.from_pretrained(esm_model_name)
    esm_model = EsmModel.from_pretrained(esm_model_name).eval()

    print(f"Cargando ProtT5 ({t5_model_path}, offline local)...", file=sys.stderr, flush=True)
    t5_tokenizer = T5Tokenizer.from_pretrained(t5_model_path, do_lower_case=False)
    t5_model = T5EncoderModel.from_pretrained(t5_model_path).eval()

    pb_full = _get_proteinbert_representation(pretrained_model_generator, input_encoder, sequence)
    esm_full = _get_esm2_embedding(esm_tokenizer, esm_model, sequence)
    t5_full = _get_prott5_embedding(t5_tokenizer, t5_model, sequence)

    proteinbert_rows, esm_rows, t5_rows = [], [], []
    for pos in positions:
        proteinbert_rows.append(pb_full)
        start = max(pos - _WINDOW_SIZE - 1, 0)
        end = min(pos + _WINDOW_SIZE, len(sequence))
        esm_rows.append(np.mean(esm_full[start:end, :], axis=0))
        t5_rows.append(t5_full[pos - 1])

    return np.concatenate([np.array(proteinbert_rows), np.array(esm_rows), np.array(t5_rows)], axis=1)


def _preprocess(feature_x: np.ndarray, stage: int, models_dir: Path) -> np.ndarray:
    with open(models_dir / f"power_transformer_{stage}.sav", "rb") as f:
        pt = pickle.load(f)
    return pt.transform(feature_x)


def _base_layer_predictions(feature_x: np.ndarray, models_dir: Path) -> np.ndarray:
    test_x = _preprocess(feature_x, 2, models_dir)
    total = np.zeros((len(test_x), 1), dtype=float)
    pickle_dir = models_dir / "base_layer_pickle_files"

    for i in range(10):
        for base_classifier in ("SVM", "XGB", "KNN"):
            with open(pickle_dir / f"{base_classifier}_base_layer_{i}.sav", "rb") as f:
                model = pickle.load(f)
            y_proba = model.predict_proba(test_x)[:, 1].reshape(-1, 1)
            total = np.concatenate((total, y_proba), axis=1)

    return np.delete(total, 0, axis=1)


def _predict(feature_x: np.ndarray, models_dir: Path) -> np.ndarray:
    """Aplica el stack de clasificadores ya entrenados (base layer + meta-SVM), sin modificaciones."""
    x = _preprocess(feature_x, 1, models_dir)
    blp = _base_layer_predictions(x, models_dir)
    x = np.concatenate((x, blp), axis=1)
    x = _preprocess(x, 3, models_dir)

    with open(models_dir / "base_layer_pickle_files" / "SVM_meta_layer.sav", "rb") as f:
        clf = pickle.load(f)
    y_pred = clf.predict(x)
    y_proba = clf.predict_proba(x)[:, 1]
    return np.column_stack([y_pred, y_proba])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Runner standalone de StackGlyEmbed (corroboracion informativa de N-glicosilacion)."
    )
    parser.add_argument("--sequence", required=True, help="Secuencia COMPLETA de la proteina (no un fragmento).")
    parser.add_argument(
        "--positions", required=True, type=int, nargs="+",
        help="Posiciones 1-based de la Asparagina de cada secuon N-X-[S/T] candidato.",
    )
    parser.add_argument(
        "--models-dir", required=True,
        help="Carpeta 'prediction/' del clon de StackGlyEmbed (power_transformer_*.sav + base_layer_pickle_files/).",
    )
    parser.add_argument("--t5-model-path", required=True, help="Ruta local a los pesos de ProtT5.")
    parser.add_argument("--esm-model-name", required=True, help="ID de HF Hub del modelo ESM-2 (offline si ya esta cacheado).")
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    seq_length = len(args.sequence)
    valid_positions = [p for p in args.positions if 1 <= p <= seq_length]
    skipped = sorted(set(args.positions) - set(valid_positions))
    if skipped:
        print(
            f"[stackglyembed_runner] {len(skipped)} posicion(es) fuera de rango (secuencia de "
            f"{seq_length} residuos), omitidas: {skipped}",
            file=sys.stderr,
        )

    if not valid_positions:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
        return 0

    models_dir = Path(args.models_dir)
    feature_x = _extract_features(args.sequence, valid_positions, args.t5_model_path, args.esm_model_name)
    predictions = _predict(feature_x, models_dir)

    rows = []
    for pos, (y_pred, y_proba) in zip(valid_positions, predictions):
        rows.append({
            "position": pos,
            "stackglyembed_veredicto": "Glicosilado" if int(y_pred) == 1 else "No glicosilado",
            "stackglyembed_score": float(y_proba),
        })

    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
