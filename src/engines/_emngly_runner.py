#!/usr/bin/env python
"""Runner standalone para EMNGly (motor real de consenso para 'n_linked_glycosylation', Camino PDB).

NUNCA se importa desde el paquete ``src`` -- requiere fair-esm/torch/
scikit-learn, dependencias del venv dedicado ``Settings.EMNGLY_PYTHON_BIN``
(nunca compartido con DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN/otros). Se
invoca EXCLUSIVAMENTE via subprocess desde ``src/engines/emngly_engine.py``,
mismo patron que ``_deepptmpred_runner.py``/``_stackglyembed_runner.py``.

## Por que existe (rol en el pipeline, decision 2026-08-06)

Reemplaza a CoNglyPred (candidato original de Decision 2, confirmado
DEFINITIVAMENTE sin pesos publicados en ningun sitio -- ver STATUS.md).
EMNGly (``github.com/StellaHxy/EMNgly``, Hou et al., Bioinformatics
39(11):btad650, 2023) SI tiene pesos reales y verificados a nivel de bytes
(ver ``src/config/settings.py``, bloque EMNGLY_*): ESM-1b (secuencia,
``site_emb``+``local_emb``, 1280+1280) + MIF de Microsoft (estructura real
sobre backbone N/CA/C del PDB, ``structure_emb``, 256) -> SVM (2816
features). Preserva la propiedad de diseno "el segundo motor de este tipo
usa estructura 3D real", ya decidida al descartar MTPrompt-PTM como
reemplazo de DeepPTMPred por ser solo-secuencia.

Este runner importa directamente el paquete ``MIF`` vendorizado dentro del
clon de EMNgly (``EMNgly/model/MIF/``, Microsoft ``protein-sequence-models``,
licencia BSD-2 permisiva verificada leyendo el original en
``github.com/microsoft/protein-sequence-models`` -- la copia de EMNgly perdio
el archivo LICENSE al vendorizarlo, por eso este proyecto documenta la
licencia real aqui en vez de asumirla de la copia) para el calculo de
``structure_emb``, y REIMPLEMENTA ``model/get_esm_embedding.py::ESMEmbeddingExtractor``
para ``site_emb``/``local_emb`` -- misma logica de chunking EXACTA (necesaria
para reproducir bit a bit lo que el SVM aprendio), pero cargando el
checkpoint ESM-1b desde una ruta ``.pt`` LOCAL en vez de
``torch.hub.load("facebookresearch/esm:main", ...)`` (el script original
pega red en cada corrida, viola la politica local-only de este proyecto,
mismo patron ya resuelto en ``_deepptmpred_runner.py``).

## Convencion de indices de ``site_emb``/``local_emb`` (verificado leyendo
   ``model/get_esm_embedding.py`` linea a linea, no asumido)

``ESMEmbeddingExtractor.extract()`` NUNCA descarta el token BOS/inicio de la
representacion cruda de ESM para el primer chunk (``delta=0`` cuando
``i==0``, la slice ``[0:j-i+1]`` incluye el indice 0) -- el array resultante
queda con el token BOS en el indice 0 y el residuo 1-based ``k`` en el indice
``k`` (NO en ``k-1``). Esto hace que ``get_site_features``
(``emb = extract([seq])[pos]``, ``pos`` 1-based) sea CORRECTO tal cual esta
escrito -- no hay off-by-one real pese a la apariencia inicial (investigacion
2026-08-06 lo habia flagged como riesgo, este runner lo confirma NO-bug
leyendo el codigo con cuidado). ``get_local_features``
(``emb = extract([local_seq])[0]``) es un diseno deliberado distinto: toma el
embedding del token BOS de la ventana local (no un residuo), como resumen
pooled del contexto -- mismo patron que usar el token [CLS] en BERT como
representacion de todo el segmento. Este runner reproduce AMBAS convenciones
tal cual, sin "arreglar" ninguna: el SVM se entreno sobre exactamente estos
features, cambiarlos crearia un desajuste train/inferencia nuevo (misma
clase de bug que ya costo 3 iteraciones reales en ``_deepptmpred_runner.py``,
ver su docstring).

## Convencion de indices de ``structure_emb`` (riesgo real distinto, ver abajo)

``model/MIF/sequence_models/pdb_utils.py::parse_PDB`` construye el array de
salida indexado por NUMERO DE RESIDUO DEL PDB menos 1 (``resn = int(resid) - 1``,
con ``min_resn`` forzado a ``min(min_resn, 0)`` y huecos rellenados con 'X'/NaN)
-- NUNCA por orden secuencial de residuos observados. El propio
``predict.py`` de EMNgly indexa ``structure_emb[pos-1]`` asumiendo
implicitamente que ``pos`` (posicion en su CSV de entrenamiento) coincide con
la numeracion real del PDB -- cierto para los PDBs de AlphaFold2 de sus
datasets (numeracion continua 1..N sin huecos), pero FALSO en general para
estructuras cristalograficas con huecos o numeracion que no arranca en 1.
Este proyecto usa la tabla de mapeo YA construida en Fase 1.5
(``src/utils/structure_parser.py::StructureRecord.position_mapping``,
columnas ``fasta_position``/``pdb_seqid``) para traducir la posicion 1-based
de ``sequence`` (ATMSEQ, la que usa el resto del pipeline) al numero de
residuo REAL del PDB antes de indexar ``structure_emb`` -- generaliza
correctamente el supuesto implicito del script original en vez de asumirlo
ciegamente (fix proactivo del riesgo #1 identificado en la investigacion de
reemplazo de CoNglyPred, 2026-08-06, antes de tener un caso real que lo
rompiera).

``parse_PDB`` con ``chain=None`` (default, igual que ``model/get_mif_embedding.py::run``)
lee TODAS las cadenas del archivo sin filtrar -- solo es seguro si
``--pdb-path`` es un PDB de UNA SOLA cadena. Este runner exige
``record.chain_pdb_path`` (Fase 1.5, nunca ``record.pdb_path`` crudo), mismo
criterio ya establecido para MeToken/DeepPTMPred.

NO PROBADO TODAVIA contra el entorno real (checkpoint ESM-1b/SVM no
descargados en esta maquina, ver STATUS.md) -- los tests mockean
``subprocess.run``.
"""

import argparse
import hashlib
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_COLUMNS = ["position", "probability"]

# --- Reimplementacion de model/get_esm_embedding.py::ESMEmbeddingExtractor ---
# Constantes identicas al script original (ver docstring del modulo).
_ESM_EMBED_LAYER = 33
_ESM_EMBED_DIM = 1280
_LOCAL_WINDOW_LEFT = 21
_LOCAL_WINDOW_RIGHT = 20


class _ESMEmbeddingExtractor:
    """Puerto local-only de ``ESMEmbeddingExtractor`` (misma logica de chunking, checkpoint LOCAL)."""

    def __init__(self, checkpoint_path: Path, device):
        import torch
        from esm import pretrained

        # ``pretrained.load_model_and_alphabet`` entra por la rama local
        # (``load_model_and_alphabet_local``, solo ``torch.load()`` sobre
        # disco) en cuanto el nombre termina en '.pt' -- mismo mecanismo ya
        # verificado en ``_deepptmpred_runner.py`` para ESM-2. Requiere el
        # archivo companero ``<checkpoint>-contact-regression.pt`` en el
        # mismo directorio (heuristica interna de fair-esm, no excluye
        # esm1b_t33_650M_UR50S), ver README.md - Seccion de instalacion.
        #
        # Bug real confirmado 2026-08-07 corriendo la carga real (no solo
        # leyendo codigo): ``esm/pretrained.py::load_model_and_alphabet_local``
        # llama ``torch.load(path, map_location="cpu")`` SIN
        # ``weights_only=False`` -- desde PyTorch 2.6 el default de
        # ``weights_only`` cambio a True, y el checkpoint de fair-esm
        # (``argparse.Namespace`` en su estado serializado) no pasa el
        # allowlist estricto. fair-esm (paquete pip, no vendorizado en este
        # repo) nunca fue actualizado para el nuevo default. Fix acotado al
        # choke point de esta clase (unico lugar del runner que carga el
        # checkpoint), no un parche global de ``torch.load`` para todo el
        # proceso -- misma filosofia que el monkeypatch de
        # ``_deepptmpred_runner.py::_load_predict_module``. Seguro porque el
        # checkpoint viene de la URL oficial de Meta
        # (dl.fbaipublicfiles.com), no de una fuente no confiable.
        _original_torch_load = torch.load

        def _torch_load_weights_only_false(*args, **kwargs):
            kwargs.setdefault("weights_only", False)
            return _original_torch_load(*args, **kwargs)

        torch.load = _torch_load_weights_only_false
        try:
            self.model, alphabet = pretrained.load_model_and_alphabet(str(checkpoint_path))
        finally:
            torch.load = _original_torch_load
        self.model.eval()
        self.model = self.model.to(device)
        self.device = device
        self.batch_converter = alphabet.get_batch_converter()
        self.max_input_len = 1022
        self.max_step_len = 511

    def extract(self, seqs):
        """Identico a ``ESMEmbeddingExtractor.extract`` del script original (ver docstring del modulo)."""
        import torch

        max_seq_len = len(seqs[0])
        representations = []
        with torch.no_grad():
            for i in range(0, max_seq_len, self.max_step_len):
                j = min(i + self.max_input_len, max_seq_len)
                delta = 0 if i == 0 else self.max_input_len - self.max_step_len
                if i > 0 and j < i + self.max_input_len:
                    delta += i + self.max_input_len - max_seq_len
                    i = max_seq_len - self.max_input_len

                batch_seqs = [("", seq[i:j]) for seq in seqs]
                _, _, batch_tokens = self.batch_converter(batch_seqs)
                batch_tokens = batch_tokens.to(self.device)
                results = self.model(batch_tokens, repr_layers=[_ESM_EMBED_LAYER])
                representations.append(results["representations"][_ESM_EMBED_LAYER][..., delta : j - i + 1, :])

                if j >= max_seq_len:
                    break

        return torch.cat(representations, dim=1)[0]


def _sequence_cache_key(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:12]


def _load_or_compute_esm_full(sequence: str, extractor: "_ESMEmbeddingExtractor", cache_dir: Path, accession: str):
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{accession}_{_sequence_cache_key(sequence)}_esm1b_full.npy"
    if cache_path.is_file():
        return np.load(cache_path)
    full = extractor.extract([sequence]).cpu().numpy()
    np.save(cache_path, full)
    return full


def _load_or_compute_structure_emb(pdb_path: Path, mif_weights: Path, cache_dir: Path, accession: str) -> np.ndarray:
    """``structure_emb`` completo (256-dim por residuo), indexado por NUMERO DE RESIDUO DEL PDB - 1.

    Cacheado por accession + hash del contenido del PDB (no solo el
    accession): un PDB distinto bajo el mismo accession nunca reusa en
    silencio el embedding viejo -- mismo motivo que el cache de
    ``_deepptmpred_runner.py`` (auditoria 2026-07-28, item 2).
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_bytes = pdb_path.read_bytes()
    pdb_hash = hashlib.sha256(pdb_bytes).hexdigest()[:12]
    cache_path = cache_dir / f"{accession}_{pdb_hash}_mif_structure.npy"
    if cache_path.is_file():
        return np.load(cache_path)

    import torch
    from MIF.sequence_models.pdb_utils import parse_PDB, process_coords
    from MIF.sequence_models.pretrained import load_model_and_alphabet

    model, collater = load_model_and_alphabet(str(mif_weights))
    model.eval()

    coords, wt, _valid_resn = parse_PDB(str(pdb_path))
    coords_dict = {"N": coords[:, 0], "CA": coords[:, 1], "C": coords[:, 2]}
    dist, omega, theta, phi = process_coords(coords_dict)

    batch = [[
        wt,
        torch.tensor(dist, dtype=torch.float),
        torch.tensor(omega, dtype=torch.float),
        torch.tensor(theta, dtype=torch.float),
        torch.tensor(phi, dtype=torch.float),
    ]]
    src, nodes, edges, connections, edge_mask = collater(batch)
    with torch.no_grad():
        outputs = model(src, nodes, edges, connections, edge_mask)

    structure_emb = outputs[0][: len(wt)].detach().cpu().numpy()
    np.save(cache_path, structure_emb)
    return structure_emb


def _fasta_position_to_pdb_seqid(position_mapping_csv: Path) -> dict:
    df = pd.read_csv(position_mapping_csv)
    return dict(zip(df["fasta_position"].astype(int), df["pdb_seqid"].astype(int)))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Runner standalone de EMNGly (consenso de n_linked_glycosylation, Camino PDB)."
    )
    parser.add_argument("--emngly-home", required=True, help="Ruta a la raiz del clon de EMNgly")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sequence", required=True, help="Secuencia ATMSEQ completa (Fase 1.5)")
    parser.add_argument("--pdb-path", required=True, help="PDB de una sola cadena (record.chain_pdb_path)")
    parser.add_argument(
        "--position-mapping-csv", required=True,
        help="CSV de Fase 1.5 (fasta_position <-> pdb_seqid) para alinear structure_emb correctamente.",
    )
    parser.add_argument("--positions", required=True, type=int, nargs="+", help="Posiciones 1-based (Asn del secuon).")
    parser.add_argument("--mif-weights", required=True)
    parser.add_argument("--esm-checkpoint", required=True)
    parser.add_argument("--svm-checkpoint", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    # model/MIF/__init__.py hace de 'MIF' un paquete top-level valido al
    # anadir 'EMNgly/model' a sys.path (usado por el 'from MIF.sequence_models...'
    # de abajo). Ademas -- verificado 2026-08-07 corriendo el import real, no
    # solo leyendo el codigo -- 'MIF/sequence_models/*.py' (pdb_utils.py,
    # pretrained.py, etc.) hacen imports BARE de 'sequence_models.xxx' (no
    # relativos ni 'MIF.sequence_models.xxx'), igual que el script original
    # 'model/get_mif_embedding.py' que ademas de 'sys.path.append("./MIF")'
    # depende de que su propio directorio ('EMNgly/model') ya este en
    # sys.path (comportamiento automatico del interprete al correr un
    # script, no presente al importar este runner via subprocess) -- por
    # eso hace falta 'EMNgly/model/MIF' en sys.path TAMBIEN, no solo
    # 'EMNgly/model'.
    sys.path.insert(0, str(Path(args.emngly_home) / "model"))
    sys.path.insert(0, str(Path(args.emngly_home) / "model" / "MIF"))

    import torch

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cache_dir = Path(args.cache_dir)

    seq_length = len(args.sequence)
    valid_positions = [p for p in args.positions if 1 <= p <= seq_length]
    skipped = sorted(set(args.positions) - set(valid_positions))
    if skipped:
        print(
            f"[emngly_runner] {len(skipped)} posicion(es) fuera de rango (secuencia de "
            f"{seq_length} residuos), omitidas: {skipped}",
            file=sys.stderr,
        )
    if not valid_positions:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
        return 0

    fasta_to_pdb_seqid = _fasta_position_to_pdb_seqid(Path(args.position_mapping_csv))
    structure_emb = _load_or_compute_structure_emb(
        Path(args.pdb_path), Path(args.mif_weights), cache_dir, args.accession
    )

    extractor = _ESMEmbeddingExtractor(Path(args.esm_checkpoint), device)
    site_full = _load_or_compute_esm_full(args.sequence, extractor, cache_dir, args.accession)

    with open(args.svm_checkpoint, "rb") as f:
        svm = pickle.load(f)

    rows = []
    skipped_alignment = []
    for pos in valid_positions:
        pdb_seqid = fasta_to_pdb_seqid.get(pos)
        structure_idx = pdb_seqid - 1 if pdb_seqid is not None else None
        if (
            pdb_seqid is None
            or structure_idx < 0
            or structure_idx >= structure_emb.shape[0]
            or pos >= site_full.shape[0]
        ):
            skipped_alignment.append(pos)
            continue

        site_vec = site_full[pos]
        left = max(0, pos - _LOCAL_WINDOW_LEFT)
        right = min(pos + _LOCAL_WINDOW_RIGHT, seq_length)
        local_window = args.sequence[left:right]
        local_vec = extractor.extract([local_window]).cpu().numpy()[0]
        struct_vec = structure_emb[structure_idx]

        feature_row = np.concatenate([site_vec, local_vec, struct_vec]).reshape(1, -1)
        probability = float(svm.predict_proba(feature_row)[0, 1])
        rows.append({"position": pos, "probability": probability})

    if skipped_alignment:
        print(
            f"[emngly_runner] {len(skipped_alignment)} posicion(es) sin mapeo PDB valido o fuera "
            f"de rango de structure_emb/site_emb, omitidas: {skipped_alignment}",
            file=sys.stderr,
        )

    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(args.out_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
