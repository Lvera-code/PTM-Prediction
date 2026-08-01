#!/usr/bin/env python
"""Runner standalone para MeToken (corroboracion informativa de TIPO, no un motor de consenso).

NUNCA se importa desde el paquete ``src`` -- requiere torch/torch_scatter/
transformers/biopython/omegaconf, dependencias SOLO presentes en el venv
dedicado de MeToken (``Settings.METOKEN_PYTHON_BIN``, distinto de
``DEEPMVP_PYTHON_BIN``/``DEEPPTMPRED_PYTHON_BIN``). Se invoca EXCLUSIVAMENTE
via subprocess desde ``src/engines/metoken_engine.py``, mismo patron que
``_deepptmpred_runner.py`` (Fase 2) y ``stackglyembed_predict_local.py`` en
``BCell-Epitope-Prediction``.

## Por que existe (rol en el pipeline, decision 2026-08-01)

Investigado a fondo (subagente Opus) que ``github.com/A4Bio/MeToken``
(ICLR 2025) es el motor estructural mas potente evaluado hasta ahora para
PTM -- consume coordenadas backbone reales (N/CA/C/O) via grafo 3D-kNN +
marcos locales por cuaternion, mucho mas rico que los 4 escalares
(SASA/phi/psi/plDDT) que usa DeepPTMPred -- pero el checkpoint PUBLICADO es
un CLASIFICADOR DE TIPO en sitios YA CONOCIDOS, no un detector de sitio:
confirmado en ``model_interface.py:40`` del repo (``valid_idx = batch['Q'] >
0 if self.hparams.with_null_ptm == 0 else ...`` -- la clase "Not a PTM type"
queda excluida de la evaluacion/entrenamiento cuando ``with_null_ptm=0``, que
es como viene el checkpoint publicado). Verificado con una corrida real
contra ``AF-P10636-F1-model_v4.pdb`` (Tau): en posiciones SIN PTM real
(prolinas, glicinas) predice tipos con alta confianza igualmente -- NO sirve
para decidir si un sitio es o no PTM.

Por eso el rol aqui es EXACTAMENTE el mismo patron no-decisorio que GlyGen
(``src/structural/glygen_client.py`` + ``pyrosetta_glycan_patch.py::check_glygen_evidence``):
corroboracion puramente informativa del TIPO en sitios que el consenso YA
acepto (``pasa_umbral=true`` en ``ptm_annotation.py``), NUNCA cambia
``pasa_umbral``/``consenso``. Ver ``src/engines/metoken_engine.py`` para el
wiring (subprocess con manejo de error no fatal) y
``src/engines/ptm_annotation.py::annotate_pdb_path`` para el punto de
enganche opcional.

## Dos bugs reales confirmados corriendo el repo (no asumidos, ver STATUS.md)

1. **``inference.py:61`` llama a ``PDB.Polypeptide.three_to_one``**, eliminado
   de Biopython en la version >=1.80 (confirmado: ``hasattr(PDB.Polypeptide,
   'three_to_one')`` -> ``False`` en Biopython 1.87, la version que instala
   ``pip install biopython`` hoy) -- revienta con ``AttributeError`` en
   cualquier entorno moderno. Reemplazado aqui por
   ``PDB.Polypeptide.protein_letters_3to1``/``protein_letters_3to1_extended``
   (los diccionarios que SI existen en Biopython moderno), monkeypatcheado
   sobre el modulo ya importado -- NO se edita ``inference.py`` (vendored).

2. **``src/metoken_model.py:213`` tiene ``device='cuda'`` hardcodeado**
   (``codebook_mask = torch.ones(len(codebook), dtype=torch.int32,
   device='cuda')``, dentro de ``MeToken_Model.__init__``) -- imposibilita
   construir el modelo en CPU. Confirmado real: sin parche,
   ``MeToken_Model(params)`` revienta con ``AssertionError: Torch not
   compiled with CUDA enabled`` en esta maquina (sin GPU, sin
   ``nvidia-smi``). Es la UNICA linea de todo ``src/`` con ``device='cuda'``
   hardcodeado (verificado por grep -- el resto de tensores usan
   ``device=x.device``/``device=index.device``, siguiendo el device del
   tensor de entrada). Como esta dentro de un ``__init__`` (no una funcion
   standalone reemplazable), se parchea envolviendo ``torch.ones`` en un
   context manager activo SOLO durante la construccion del modelo: si se
   pide ``device='cuda'`` y CUDA no esta disponible, redirige a ``'cpu'``;
   cualquier otra llamada a ``torch.ones`` (dentro o fuera de ese bloque) no
   se ve afectada. No se edita ``src/metoken_model.py`` (vendored).

Ambos parches verificados real: corriendo SIN parche 2 revienta
(``AssertionError``); corriendo CON ambos parches sobre
``examples/Q16613.pdb`` (el ejemplo del propio repo) reproduce EXACTO el
resultado documentado en su ``quick_inference.ipynb``
(``PTM type at the position 31 is Phosphorylation``) -- confirma que los
parches no alteran el comportamiento numerico del modelo, solo lo hacen
correr en este entorno.

## Deteccion de la clase "no-PTM"/"rare" (24 clases reales, no 26)

``src/constant.py::PTMtype_list`` tiene 26 entradas: indice 0 = "Not a PTM
type" (la clase null, enmascarada en entrenamiento -- ver arriba), indices
1-24 = los 24 tipos de PTM reales que el modelo distingue, indice 25 = "in
Rare PTM Types" (cubo de PTMs raros agrupados, no un tipo especifico
interpretable). Este runner excluye AMBOS indices (0 y 25) al buscar el tipo
con mayor probabilidad -- "de las 24 clases", como pide la tarea -- y reporta
la probabilidad cruda (softmax sobre las 26 clases completas, sin
renormalizar) del indice ganador entre esas 24, no una probabilidad
renormalizada.

## Offline (mismo caveat que otros motores basados en ESM de este proyecto)

``AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")`` se llama
DOS veces en el repo (``src/metoken_model.py`` y
``src/datasets/featurizer.py``, un tokenizer nuevo cada vez -- no es el
encoder ESM-2 completo, MeToken usa su propio ``nn.Embedding`` entrenado
desde cero, ``wo_esm``, no representaciones ESM reales pese al nombre del
atributo) -- descarga de HF Hub la primera vez, cacheable localmente despues
(mismo cache que ya usa ``StackGlyEmbed`` del proyecto hermano para el mismo
checkpoint, ver ``HF_HUB_OFFLINE``/``TRANSFORMERS_OFFLINE`` abajo, patron
identico a
``BCell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py``).

## torch_scatter sin wheel prebuilt

``src/metoken_module.py`` importa ``torch_scatter`` (``scatter_sum``,
``scatter_softmax``, ``scatter_mean``) -- sin wheel prebuilt para la
combinacion torch/CPU/Python de este entorno en ``data.pyg.org``
(confirmado: el indice de wheels solo lista variantes hasta
``torch-2.1.0+cpu``, esta maquina tiene ``torch==2.13.0+cpu``), asi que
``pip install --no-build-isolation torch_scatter`` compila desde fuente
(extension C++/CPU, ~pocos minutos en esta maquina real, no los ~15 min
estimados de antemano -- ver STATUS.md).
"""

import argparse
import contextlib
import os
import sys
from pathlib import Path

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
