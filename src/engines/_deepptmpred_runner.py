#!/usr/bin/env python
"""Runner standalone para DeepPTMPred (Fase 3a, motor 2/2 del consenso Camino PDB).

NUNCA se importa desde el paquete ``src`` -- requiere torch/tensorflow/
tensorflow-addons/pyrosetta/fair-esm, dependencias SOLO presentes en el venv
dedicado de DeepPTMPred (``Settings.DEEPPTMPRED_PYTHON_BIN``, distinto del
venv de DeepMVP). Se invoca EXCLUSIVAMENTE via subprocess desde
``src/engines/deepptmpred_engine.py``, mismo patron que
``stackglyembed_predict_local.py`` en BCell-Epitope-Prediction (script de
integracion propio que vive en ``src/engines/`` pero corre fuera del venv
principal del pipeline).

Por que existe este runner en vez de invocar los scripts del repo
directamente (a diferencia de DeepMVP, que si tiene un CLI real): verificado
leyendo el codigo fuente de github.com/kuikui-wang/DeepPTMPred el
2026-07-27, ni ``predict.py`` ni ``e2_single_data.py`` tienen CLI -- ambos
hardcodean ``ptm_type``/``pdb_path``/``protein_id``/la ruta del checkpoint
ESM dentro de su bloque ``if __name__ == "__main__":``. Este runner importa
las dos clases SI parametrizadas correctamente de ``predict.py``
(``PredictConfig``, ``PTMPredictor``, ambas reciben sus argumentos por
constructor) y REIMPLEMENTA la extraccion de features ESM-2 en vez de
llamar a ``e2_single_data.py::extract_full_sequence_esm``: esa funcion tiene
un bug real confirmado -- redefine ``custom_checkpoint_path`` como variable
LOCAL con una ruta absoluta hardcodeada de AutoDL (``/root/autodl-tmp/...``),
ignorando cualquier valor pasado como parametro o de modulo.

Tambien evita ``predict.py::extract_protein_id_from_pdb_path`` /
``extract_sequence_from_pdb`` (que exigen un nombre de archivo estilo
AlphaFold, p. ej. ``AF-P12345-F1-model_v4.pdb``, y hacen su propia
extraccion ATMSEQ redundante con ``src.utils.structure_parser``): este
runner recibe el accession y la secuencia YA saneados por Fase 1.5 como
argumentos explicitos, garantizando que DeepMVP y DeepPTMPred reporten
exactamente la misma numeracion de posicion para el mismo accession.

NO PROBADO TODAVIA contra el entorno real (PyRosetta/TensorFlow 2.15/
checkpoint ESM-2 no instalados en esta maquina) -- ver STATUS.md.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Columnas de salida: 'probability' es el score crudo de DeepPTMPred (se
# conserva siempre, decision de arquitectura 2026-07-27); la columna
# 'prediction' propia del repo (cutoff 0.5 hardcodeado, no calibrado contra
# ningun validation set) se descarta deliberadamente -- el filtro real lo
# aplica el nucleo de Fase 3, no este runner.
OUTPUT_COLUMNS = ["protein_id", "position", "residue", "probability", "ptm_type"]


def _extract_esm_features(sequence: str, checkpoint_path: Path, esm_dim: int = 1280) -> np.ndarray:
    """Reimplementacion propia de e2_single_data.py::extract_full_sequence_esm.

    Misma logica (chunking a 1022 tokens, capa de representacion 33,
    descarte de CLS/SEP), pero ``checkpoint_path`` SI se respeta -- el
    original la ignora (ver docstring del modulo).
    """
    import torch
    from esm import pretrained

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    with torch.serialization.safe_globals([argparse.Namespace]):
        model, alphabet = pretrained.load_model_and_alphabet(str(checkpoint_path))
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    max_len = 1022  # dejar espacio para tokens CLS/SEP (limite real de ESM-2)
    chunks = [sequence[i : i + max_len] for i in range(0, len(sequence), max_len)]
    full_features = np.zeros((len(sequence), esm_dim))

    for i, chunk in enumerate(chunks):
        data = [(f"chunk{i}", chunk)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[33])
            features = results["representations"][33][0, 1:-1].cpu().numpy()
        start = i * max_len
        end = start + len(chunk)
        full_features[start:end] = features[: len(chunk)]

    return full_features


def _load_predict_module(train_ptm_dir: Path):
    """Inserta ``train_ptm_dir`` en sys.path e importa el modulo ``predict`` del repo.

    Import diferido (dentro de la funcion, no a nivel de modulo): las
    dependencias pesadas de ``predict.py`` (tensorflow, pyrosetta,
    tensorflow_addons) solo deben cargarse cuando el runner realmente se
    ejecuta, nunca al parsear este archivo.
    """
    sys.path.insert(0, str(train_ptm_dir))
    import predict

    return predict


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Runner standalone de DeepPTMPred (un tipo de PTM por invocacion)."
    )
    parser.add_argument("--train-ptm-dir", required=True, help="Ruta a DeepPTMPred/pred/train_PTM")
    parser.add_argument("--protein-id", required=True)
    parser.add_argument("--sequence", required=True, help="Secuencia ATMSEQ ya saneada por Fase 1.5")
    parser.add_argument("--pdb-path", required=True, help="PDB de una sola cadena (Fase 1.5)")
    parser.add_argument("--ptm-type", required=True, choices=[
        "phosphorylation", "acetylation", "ubiquitination", "hydroxylation",
        "gamma_carboxyglutamic_acid", "lys_methylation", "malonylation",
        "arg_methylation", "crotonylation", "succinylation", "glutathionylation",
        "sumoylation", "s_nitrosylation", "glutarylation", "citrullination",
        "o_linked_glycosylation", "n_linked_glycosylation",
    ])
    parser.add_argument("--esm-checkpoint", required=True, help="Ruta a esm2_t33_650M_UR50D.pt")
    parser.add_argument("--custom-esm-dir", required=True, help="Cache de features ESM (.npz por accession)")
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    train_ptm_dir = Path(args.train_ptm_dir)
    predict = _load_predict_module(train_ptm_dir)

    custom_esm_dir = Path(args.custom_esm_dir)
    custom_esm_dir.mkdir(parents=True, exist_ok=True)
    esm_path = custom_esm_dir / f"{args.protein_id}_full_esm.npz"
    if not esm_path.is_file():
        features = _extract_esm_features(args.sequence, Path(args.esm_checkpoint))
        np.savez_compressed(
            esm_path, features=features, protein_id=args.protein_id,
            sequence=args.sequence, length=len(args.sequence),
        )

    # project_root = DeepPTMPred/ (train_ptm_dir = DeepPTMPred/pred/train_PTM)
    project_root = train_ptm_dir.parent.parent
    config = predict.PredictConfig(ptm_type=args.ptm_type, project_root=str(project_root))
    # Se sobreescribe DESPUES de construir el config (que fija su propio
    # default relativo al project_root) y ANTES de construir el predictor
    # (que pasa 'config' por referencia a su data loader interno): el
    # atributo se lee en el momento de la llamada, no se copia antes.
    config.custom_esm_dir = str(custom_esm_dir)

    predictor = predict.PTMPredictor(config)

    target_aa = config.target_aa
    positions = [i + 1 for i, aa in enumerate(args.sequence) if aa in target_aa]

    if not positions:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
        return 0

    results_df = predictor.predict_ptm_sites(
        args.protein_id, args.sequence, positions, pdb_path=args.pdb_path
    )
    results_df["protein_id"] = args.protein_id
    results_df["ptm_type"] = args.ptm_type
    results_df[OUTPUT_COLUMNS].to_csv(args.out_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
