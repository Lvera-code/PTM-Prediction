#!/usr/bin/env python
"""Standalone runner for MeToken (informative corroboration of TYPE, not a consensus engine).

It is NEVER imported from the ``src`` package: it requires torch/torch_scatter/
transformers/biopython/omegaconf, dependencies present ONLY in MeToken's
dedicated venv (``Settings.METOKEN_PYTHON_BIN``, distinct from
``DEEPMVP_PYTHON_BIN``/``DEEPPTMPRED_PYTHON_BIN``). It is invoked EXCLUSIVELY
via subprocess from ``src/engines/metoken_engine.py``, the same pattern as
``_deepptmpred_runner.py`` (Phase 2).

## Why it exists (its role in the pipeline, decision 2026-08-01)

It was investigated in depth (Opus subagent) that ``github.com/A4Bio/MeToken``
(ICLR 2025) is the most powerful structural engine evaluated so far for
PTM: it consumes real backbone coordinates (N/CA/C/O) via a 3D-kNN graph +
local quaternion frames, far richer than the 4 scalars
(SASA/phi/psi/plDDT) that DeepPTMPred uses. However, the PUBLISHED checkpoint is
a TYPE CLASSIFIER on ALREADY KNOWN sites, not a site detector:
confirmed in the repo's ``model_interface.py:40`` (``valid_idx = batch['Q'] >
0 if self.hparams.with_null_ptm == 0 else ...``, where the "Not a PTM type"
class is excluded from evaluation/training when ``with_null_ptm=0``, which is
how the published checkpoint comes). Verified with a real run
against ``AF-P10636-F1-model_v4.pdb`` (Tau): at positions WITHOUT a real PTM
(prolines, glycines) it predicts types with high confidence all the same, so it is
NOT usable for deciding whether a site is a PTM or not.

That is why its role here is the same non-deciding pattern as the secretory
pathway corroboration (``src/structural/uniprot_localization_client.py``):
purely informative corroboration of the TYPE at sites that the consensus has
ALREADY accepted (``pasa_umbral=true`` in ``ptm_annotation.py``), NEVER changing
``pasa_umbral``/``consenso``. See ``src/engines/metoken_engine.py`` for the
wiring (subprocess with non-fatal error handling) and
``src/engines/ptm_annotation.py::annotate_pdb_path`` for the optional
attachment point.

## Two real bugs confirmed by running the repo (not assumed)

1. **``inference.py:61`` calls ``PDB.Polypeptide.three_to_one``**, removed
   from Biopython in version >=1.80 (confirmed: ``hasattr(PDB.Polypeptide,
   'three_to_one')`` -> ``False`` in Biopython 1.87, the version that
   ``pip install biopython`` installs today), so it fails with an ``AttributeError`` in
   any modern environment. It is replaced here by
   ``PDB.Polypeptide.protein_letters_3to1``/``protein_letters_3to1_extended``
   (the dictionaries that DO exist in modern Biopython), monkeypatched
   over the already imported module. ``inference.py`` is NOT edited (vendored).

2. **``src/metoken_model.py:213`` has ``device='cuda'`` hardcoded**
   (``codebook_mask = torch.ones(len(codebook), dtype=torch.int32,
   device='cuda')``, inside ``MeToken_Model.__init__``), which makes it impossible
   to build the model on CPU. Confirmed real: without a patch,
   ``MeToken_Model(params)`` fails with ``AssertionError: Torch not
   compiled with CUDA enabled`` on this machine (no GPU, no
   ``nvidia-smi``). It is the ONLY line in all of ``src/`` with ``device='cuda'``
   hardcoded (verified by grep; the rest of the tensors use
   ``device=x.device``/``device=index.device``, following the device of the
   input tensor). Since it is inside an ``__init__`` (not a standalone
   replaceable function), it is patched by wrapping ``torch.ones`` in a
   context manager active ONLY during the construction of the model: if
   ``device='cuda'`` is requested and CUDA is not available, it redirects to ``'cpu'``;
   any other call to ``torch.ones`` (inside or outside that block) is
   unaffected. ``src/metoken_model.py`` is NOT edited (vendored).

Both patches verified for real: running WITHOUT patch 2 fails
(``AssertionError``); running WITH both patches over
``examples/Q16613.pdb`` (the repo's own example) reproduces EXACTLY the
result documented in its ``quick_inference.ipynb``
(``PTM type at the position 31 is Phosphorylation``), which confirms that the
patches do not alter the model's numerical behaviour, they only make it
run in this environment.

## Detection of the "no-PTM"/"rare" class (24 real classes, not 26)

``src/constant.py::PTMtype_list`` has 26 entries: index 0 = "Not a PTM
type" (the null class, masked during training, see above), indices
1-24 = the 24 real PTM types that the model distinguishes, index 25 = "in
Rare PTM Types" (an aggregate class of grouped rare PTMs, not a specific
interpretable type). This runner excludes BOTH indices (0 and 25) when searching for the
type with the highest probability, that is, "out of the 24 classes" as the task
requires, and reports the raw probability (softmax over the complete 26
classes, without renormalizing) of the winning index among those 24, not a
renormalized probability.

## Offline (the same caveat as this project's other ESM-based engines)

``AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")`` is called
TWICE in the repo (``src/metoken_model.py`` and
``src/datasets/featurizer.py``, a new tokenizer each time; it is not the
complete ESM-2 encoder, since MeToken uses its own ``nn.Embedding`` trained
from scratch, ``wo_esm``, not real ESM representations despite the attribute's
name), which downloads from HF Hub the first time and is cacheable locally afterwards
(see ``HF_HUB_OFFLINE``/``TRANSFORMERS_OFFLINE`` below).

## torch_scatter without a prebuilt wheel

``src/metoken_module.py`` imports ``torch_scatter`` (``scatter_sum``,
``scatter_softmax``, ``scatter_mean``), which has no prebuilt wheel for the
torch/CPU/Python combination of this environment on ``data.pyg.org``
(confirmed: the wheel index only lists variants up to
``torch-2.1.0+cpu``, while this machine has ``torch==2.13.0+cpu``), so
``pip install --no-build-isolation torch_scatter`` compiles from source
(a C++/CPU extension, a few minutes on this real machine, not the ~15 min
estimated beforehand).
"""

import argparse
import contextlib
import os
import sys
from pathlib import Path

# Offline SOLO si "facebook/esm2_t33_650M_UR50D" ya esta en la cache local de
# HF Hub -- en una maquina que ya lo descargo una vez (dev local) esto evita
# tocar red en cada corrida, pero forzarlo siempre e incondicionalmente (como
# hacia esto antes) rompe la primera corrida en una maquina nueva sin cache
# (p. ej. un runtime de Colab recien creado): LocalEntryNotFoundError, sin
# forma de descargar nunca. Layout de cache verificado contra la
# documentacion de huggingface_hub: "<HF_HOME>/hub/models--<org>--<name>".
_HF_HOME = Path(os.environ.get("HF_HOME", Path.home() / ".cache" / "huggingface"))
_ESM2_CACHE_DIR = _HF_HOME / "hub" / "models--facebook--esm2_t33_650M_UR50D"
if _ESM2_CACHE_DIR.is_dir():
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

import pandas as pd  # noqa: E402

OUTPUT_COLUMNS = ["position", "metoken_type", "metoken_probability"]

# Indices de PTMtype_list a excluir siempre al buscar el tipo con mayor
# probabilidad (ver docstring del modulo): 0 = "Not a PTM type" (clase null
# enmascarada en entrenamiento, nunca aprendio a dispararse de forma
# fiable), 25 = "in Rare PTM Types" (cubo de tipos raros agrupados, no
# interpretable como un tipo especifico).
_NULL_CLASS_INDEX = 0
_RARE_CLASS_INDEX = 25


@contextlib.contextmanager
def _force_cpu_ones():
    """Parche real 2 (ver docstring del modulo): ``torch.ones(..., device='cuda')``
    hardcodeado en ``src/metoken_model.py::MeToken_Model.__init__`` (linea 213).

    Activo SOLO durante la construccion del modelo -- envuelve ``torch.ones``
    para redirigir ``device='cuda'`` a ``'cpu'`` unicamente cuando CUDA no
    esta disponible, restaurando la funcion original al salir del bloque
    (nunca deja el monkeypatch activo mas alla de lo necesario). No se edita
    ``src/metoken_model.py`` (vendored).
    """
    import torch

    original_ones = torch.ones

    def _patched_ones(*args, **kwargs):
        if kwargs.get("device") == "cuda" and not torch.cuda.is_available():
            kwargs["device"] = "cpu"
        return original_ones(*args, **kwargs)

    torch.ones = _patched_ones
    try:
        yield
    finally:
        torch.ones = original_ones


def _patch_three_to_one() -> None:
    """Parche real 1 (ver docstring del modulo): ``PDB.Polypeptide.three_to_one``
    eliminado de Biopython >=1.80, usado por ``inference.py::get_seq_str``.

    Monkeypatchea el atributo sobre el modulo ``Bio.PDB.Polypeptide`` ya
    importado (no se edita ``inference.py``, vendored) -- ``inference.py``
    lo referencia como ``PDB.Polypeptide.three_to_one(...)`` en el cuerpo de
    una funcion, resuelto dinamicamente en cada llamada, asi que el parche
    aplica sin importar el orden de imports.
    """
    from Bio import PDB

    def _three_to_one(resname: str) -> str:
        from Bio.PDB.Polypeptide import protein_letters_3to1, protein_letters_3to1_extended

        if resname in protein_letters_3to1:
            return protein_letters_3to1[resname]
        if resname in protein_letters_3to1_extended:
            return protein_letters_3to1_extended[resname]
        raise KeyError(resname)

    PDB.Polypeptide.three_to_one = _three_to_one


def _load_metoken_modules(repo_dir: Path):
    """Inserta ``repo_dir`` en ``sys.path`` e importa los modulos del repo vendorizado.

    ``repo_dir`` se inserta en la posicion 0 de ``sys.path`` -- como este
    script SIEMPRE se ejecuta como archivo (``python _metoken_runner.py
    ...``, nunca ``python -c``/``-m`` desde la raiz de este proyecto),
    ``sys.path[0]`` ya es el directorio de este propio script
    (``PTM-Prediction/src/engines/``), que no tiene ningun subdirectorio
    ``src/`` propio -- verificado real que NO hay colision con el paquete
    ``src`` de este proyecto (que si tiene ``src/__init__.py`` real, a
    diferencia del ``src/`` de MeToken, que es un namespace package sin
    ``__init__.py`` en su raiz): ``import src.metoken_model`` dentro de este
    proceso aislado resuelve siempre al ``src/`` de MeToken, confirmado
    ejecutando este runner como script real (no en modo interactivo, donde
    el directorio de trabajo SI se agrega a ``sys.path`` y produciria una
    colision real -- cuidado si se prueba manualmente con ``python -c``
    desde la raiz de este proyecto).
    """
    sys.path.insert(0, str(repo_dir))
    _patch_three_to_one()

    from omegaconf import OmegaConf
    from src.constant import PTMtype_list
    from src.datasets.featurizer import featurize
    from src.metoken_model import MeToken_Model

    from inference import apply_ptm_indices, get_seq_str

    return {
        "OmegaConf": OmegaConf,
        "PTMtype_list": PTMtype_list,
        "featurize": featurize,
        "MeToken_Model": MeToken_Model,
        "apply_ptm_indices": apply_ptm_indices,
        "get_seq_str": get_seq_str,
    }


def _top_type_excluding_masked(probs_row, ptm_type_list):
    """Indice/etiqueta/probabilidad del tipo con mayor probabilidad, excluyendo
    la clase null (indice 0) y la clase "rare" (indice 25) -- ver docstring
    del modulo. ``probs_row`` es el vector softmax de 26 clases para UNA
    posicion.
    """
    best_idx = None
    best_prob = -1.0
    for idx, prob in enumerate(probs_row):
        if idx in (_NULL_CLASS_INDEX, _RARE_CLASS_INDEX):
            continue
        if prob > best_prob:
            best_prob = prob
            best_idx = idx
    return ptm_type_list[best_idx], float(best_prob)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Runner standalone de MeToken (corroboracion informativa de tipo)."
    )
    parser.add_argument("--repo-dir", required=True, help="Ruta al clon de github.com/A4Bio/MeToken")
    parser.add_argument("--checkpoint-path", required=True, help="Ruta a pretrained_model/checkpoint.ckpt")
    parser.add_argument("--pdb-path", required=True, help="PDB de una sola cadena (Fase 1.5)")
    parser.add_argument("--chain-id", default="A", help="Cadena a leer del PDB (default 'A')")
    parser.add_argument(
        "--positions", required=True, type=int, nargs="+",
        help="Posiciones 1-based ya aceptadas por el consenso (pasa_umbral=true) a corroborar.",
    )
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    repo_dir = Path(args.repo_dir)
    mods = _load_metoken_modules(repo_dir)

    protein_data = mods["get_seq_str"](args.pdb_path, chain_id=args.chain_id)
    seq_length = len(protein_data["seq"])

    # Posiciones pedidas que caen dentro de la secuencia realmente leida del
    # PDB (get_seq_str puede devolver menos residuos que el pdb_path si hay
    # huecos/residuos no estandar sin coordenadas N/CA/C/O completas) -- las
    # que no caen dentro del rango se omiten con un aviso, nunca fatal.
    valid_positions = [p for p in args.positions if 1 <= p <= seq_length]
    skipped = sorted(set(args.positions) - set(valid_positions))
    if skipped:
        print(
            f"[metoken_runner] {len(skipped)} posicion(es) fuera de rango (secuencia leida "
            f"tiene {seq_length} residuos), omitidas: {skipped}",
            file=sys.stderr,
        )

    if not valid_positions:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
        return 0

    zero_based = [p - 1 for p in valid_positions]
    protein_data = mods["apply_ptm_indices"](protein_data, zero_based)
    data = mods["featurize"]([protein_data])

    import torch

    checkpoint = torch.load(args.checkpoint_path, map_location="cpu", weights_only=False)
    params = mods["OmegaConf"].load(str(repo_dir / "configs" / "MeToken.yaml"))

    with _force_cpu_ones():
        model = mods["MeToken_Model"](params)
        model.load_state_dict(checkpoint)
    model.eval()

    with torch.no_grad():
        result = model(data)
    probs = result["log_probs"].softmax(dim=-1).cpu().tolist()

    rows = []
    for position in valid_positions:
        metoken_type, probability = _top_type_excluding_masked(probs[position - 1], mods["PTMtype_list"])
        rows.append({"position": position, "metoken_type": metoken_type, "metoken_probability": probability})

    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
