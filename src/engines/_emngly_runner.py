#!/usr/bin/env python
"""Standalone runner for EMNGly (the real consensus engine for 'n_linked_glycosylation', PDB path).

It is NEVER imported from the ``src`` package: it requires fair-esm/torch/
scikit-learn, dependencies of the dedicated venv ``Settings.EMNGLY_PYTHON_BIN``
(never shared with DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN/others). It is
invoked EXCLUSIVELY via subprocess from ``src/engines/emngly_engine.py``,
the same pattern as ``_deepptmpred_runner.py``.

## Why it exists (its role in the pipeline, decision 2026-08-06)

It replaces CoNglyPred (the original Decision 2 candidate, confirmed
DEFINITIVELY to have no published weights anywhere). EMNGly
(``github.com/StellaHxy/EMNgly``, Hou et al., Bioinformatics
39(11):btad650, 2023) DOES have real weights, verified at byte level
(see ``src/config/settings.py``, the EMNGLY_* block): ESM-1b (sequence,
``site_emb``+``local_emb``, 1280+1280) + Microsoft's MIF (real structure
over the PDB backbone N/CA/C, ``structure_emb``, 256) -> SVM (2816
features). It preserves the design property that "the second engine of this type
uses real 3D structure", already decided when discarding MTPrompt-PTM as a
replacement for DeepPTMPred for being sequence-only.

This runner imports directly the ``MIF`` package vendored inside the
EMNgly clone (``EMNgly/model/MIF/``, Microsoft ``protein-sequence-models``,
a permissive BSD-2 licence verified by reading the original at
``github.com/microsoft/protein-sequence-models``; EMNgly's copy lost
the LICENSE file when vendoring it, which is why this project documents the real
licence here instead of assuming it from the copy) for computing
``structure_emb``, and REIMPLEMENTS ``model/get_esm_embedding.py::ESMEmbeddingExtractor``
for ``site_emb``/``local_emb``, with the EXACT same chunking logic (necessary
in order to reproduce bit for bit what the SVM learned), but loading the
ESM-1b checkpoint from a LOCAL ``.pt`` path instead of
``torch.hub.load("facebookresearch/esm:main", ...)`` (the original script
accesses the network on every run, which violates this project's local-only policy,
the same pattern already resolved in ``_deepptmpred_runner.py``).

## Index convention of ``site_emb``/``local_emb`` (verified by reading
   ``model/get_esm_embedding.py`` line by line, not assumed)

``ESMEmbeddingExtractor.extract()`` NEVER discards the BOS/start token from the
raw ESM representation for the first chunk (``delta=0`` when
``i==0``, and the slice ``[0:j-i+1]`` includes index 0), so the resulting array
is left with the BOS token at index 0 and the 1-based residue ``k`` at index
``k`` (NOT at ``k-1``). This makes ``get_site_features``
(``emb = extract([seq])[pos]``, ``pos`` 1-based) CORRECT exactly as it is
written: there is no real off-by-one despite the initial appearance (the
2026-08-06 investigation had flagged it as a risk, and this runner confirms it is
NOT a bug after reading the code carefully). ``get_local_features``
(``emb = extract([local_seq])[0]``) is a deliberately different design: it takes the
embedding of the local window's BOS token (not a residue), as a pooled
summary of the context, the same pattern as using the [CLS] token in BERT as a
representation of the whole segment. This runner reproduces BOTH conventions
as they are, without "fixing" either: the SVM was trained over exactly these
features, and changing them would create a new train/inference mismatch
(the same class of bug that already cost 3 real iterations in ``_deepptmpred_runner.py``,
see its docstring).

## Index convention of ``structure_emb`` (a different real risk, see below)

``model/MIF/sequence_models/pdb_utils.py::parse_PDB`` builds the output array
indexed by PDB RESIDUE NUMBER minus 1 (``resn = int(resid) - 1``,
with ``min_resn`` forced to ``min(min_resn, 0)`` and gaps filled with 'X'/NaN),
NEVER by sequential order of observed residues. EMNgly's own
``predict.py`` indexes ``structure_emb[pos-1]`` implicitly assuming
that ``pos`` (the position in its training CSV) matches
the real PDB numbering, which is true for the AlphaFold2 PDBs of its
datasets (continuous 1..N numbering without gaps) but FALSE in general for
crystallographic structures with gaps or with numbering that does not start at 1.
This project uses the mapping table ALREADY built in Phase 1.5
(``src/utils/structure_parser.py::StructureRecord.position_mapping``,
columns ``fasta_position``/``pdb_seqid``) to translate the 1-based position
of ``sequence`` (ATMSEQ, the one the rest of the pipeline uses) to the REAL PDB
residue number before indexing ``structure_emb``, which correctly
generalizes the original script's implicit assumption instead of assuming it
blindly (a proactive fix of risk #1 identified in the CoNglyPred
replacement investigation, 2026-08-06, before having a real case that
broke it).

``parse_PDB`` with ``chain=None`` (the default, as in ``model/get_mif_embedding.py::run``)
reads ALL the chains of the file without filtering, which is only safe if
``--pdb-path`` is a SINGLE-chain PDB. This runner therefore requires
``record.chain_pdb_path`` (Phase 1.5, never the raw ``record.pdb_path``), the same
criterion already established for MeToken/DeepPTMPred.

Both go/no-go checks were subsequently run against the real environment
(2026-08-07, with the ESM-1b checkpoint and the SVM downloaded): the
``structure_emb`` alignment was verified against real GlyGen sites in a PDB
with gaps, and the published MCC was reproduced and exceeded over the
independent N-GlyDE set.
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
