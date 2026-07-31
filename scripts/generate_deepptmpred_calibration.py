#!/usr/bin/env python
"""Calibra el umbral de probabilidad de DeepPTMPred con datos reales, por tipo de PTM.

Standalone -- requiere el stack de DeepPTMPred (TensorFlow 2.15, PyTorch,
PyRosetta, fair-esm), correr con el conda env `deepptmpred`
(``conda run -n deepptmpred python scripts/generate_deepptmpred_calibration.py``),
nunca con el interprete principal del pipeline.

## Por que existe

``Settings.DEEPPTMPRED_MIN_PROBABILITY`` (0.5) es un cutoff heredado del propio
repo, hardcodeado en ``predict.py`` sin documentar de donde sale ni contra que
validation set se calibro (a diferencia del ``fpr`` de DeepMVP, que si viene de
un archivo de validacion real). Este script corre el modelo REAL (los ``.h5``
ya instalados, via ``predict.PTMPredictor``, el mismo codigo que usa
``src/engines/_deepptmpred_runner.py`` en produccion) sobre una MUESTRA
representativa de sitios reales con etiqueta conocida, y mide AUROC/umbral
optimo (Youden's J) por cada uno de los 17 tipos soportados.

## De donde sale el dataset etiquetado

``DeepPTMPred/data/ptm_data.csv.gz`` -- mismo nombre/esquema de columnas que
``ptm_data.csv`` referenciado en el propio README de DeepPTMPred (enlace de
Google Drive, no descargable via HTTP simple), obtenido en cambio desde
``github.com/meilerlab/PTMPrediction`` (verificado real: mismas 18 columnas
que consume ``DeepPTMPred/pred/train_PTM/data_loader.py::load_dataset``
-- ``label, entry, pos, trunc_seq, sasa, phi, psi, secstruct, local_plDDT,
avg_plDDT, ptm, mod_aa, ...`` -- practicamente con total seguridad la misma
fuente). 188278 filas reales, 26796 entries unicos (verificado con pandas
2026-07-30 -- la cifra de 376557 filas en STATUS.md, previa a esta
verificacion, era incorrecta).

## Por que NO se usan las columnas estructurales del propio CSV

``sasa/phi/psi/secstruct/local_plDDT`` del CSV son snapshots precomputados
(por los autores, con metodo no documentado) usados solo en su pipeline de
ENTRENAMIENTO (``data_loader.py::load_dataset``), no en el de INFERENCIA real
que corre en produccion (``predict.py::PTMPredictor``, el que llama nuestro
runner). Ademas se confirmo un bug real en ``data_loader.py`` (linea 139-141):
``phi_center``/``psi_center`` se extraen con
``half_window = (max(window_sizes)-1)//2 = 25``, pero el array ``phi``/``psi``
del CSV solo tiene 11 elementos -- la condicion ``len(x) > half_window`` es
SIEMPRE falsa, asi que ``phi_center``/``psi_center`` fueron 0.0 constante
durante el entrenamiento real del modelo shippeado. Para calibrar el
comportamiento REAL de ese modelo en produccion (lo que le importa a este
pipeline, no una version idealizada que nunca se entreno), hay que replicar
la ruta de INFERENCIA (estructura PyRosetta real + ventanas desde la
secuencia completa + ESM-2 real), no la de entrenamiento -- de ahi que este
script descargue estructura/secuencia real por proteina en vez de reusar las
columnas del CSV directamente.

## Fuentes externas (secuencia + estructura real por proteina)

Solo se muestrean entries tipo UniProt (>=97% del dataset, ver
``_looks_like_pdb_entry``; los ~1918 sitios con entry tipo PDB p.ej.
``7V6D_A`` quedan fuera de esta muestra, decision explicita para mantener
una unica ruta de fetch). Por accession:
``https://www.alphafold.ebi.ac.uk/api/prediction/{accession}`` (verificado
real 2026-07-30) devuelve tanto la secuencia EXACTA usada por el modelo como
``pdbUrl`` -- confirmado que AlphaFold DB ya sirve ``v6`` (no ``v4``, la
version que aparece en los nombres de archivo ya presentes en
``DeepPTMPred/data/*.pdb`` de sesiones anteriores -- esos archivos siguen
siendo validos, solo que si se necesita descargar una proteina NUEVA hay que
resolver la URL real via este API en vez de asumir el patron
``AF-{acc}-F1-model_v4.pdb``, que ahora devuelve 404). Se usa la ``sequence``
que devuelve el API (no un fetch aparte a UniProt) para garantizar que
secuencia y estructura sean exactamente consistentes.

## Bug real de `predict.py` encontrado y evitado (sin tocar el vendored)

``PTMPredictor.data_loader.rosetta_calc`` se crea una sola vez (`if
self.rosetta_calc is None`) y nunca se reinicializa para un ``pdb_path``
nuevo -- si se reusa el mismo ``PTMPredictor`` (mismo modelo .h5 cargado)
para varias proteinas distintas, TODAS menos la primera reciben
silenciosamente las features estructurales de la PRIMERA proteina. Este
script reinicia ``predictor.data_loader.rosetta_calc = None`` antes de cada
proteina nueva para forzar una ``PyRosettaCalculator`` fresca -- necesario
porque, a diferencia de produccion (que crea un ``PTMPredictor`` nuevo por
invocacion, un subproceso = una proteina), aqui se reusa el mismo predictor
para muchas proteinas por tipo, evitando recargar el modelo cientos de veces.
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from sklearn.metrics import roc_auc_score, roc_curve

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / "src" / "engines"))
from _deepptmpred_runner import _extract_esm_features, _load_predict_module  # noqa: E402

ALPHAFOLD_API_URL = "https://www.alphafold.ebi.ac.uk/api/prediction/{acc}"
PDB_ENTRY_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]+$")

# Identico a DeepPTMPred/pred/train_PTM/data_loader.py::PTMDataLoader
# (ptm_aa_map/ptm_patterns) -- se duplica aqui a proposito (mismo criterio ya
# usado en _deepptmpred_runner.py para Settings.DEEPPTMPRED_PTM_TYPES): este
# script corre en el venv dedicado, nunca importa `src`.
PTM_AA_MAP = {
    "phosphorylation": ["S", "T"],
    "acetylation": ["K"],
    "ubiquitination": ["K"],
    "hydroxylation": ["P"],
    "gamma_carboxyglutamic_acid": ["E"],
    "lys_methylation": ["K"],
    "malonylation": ["K"],
    "arg_methylation": ["R"],
    "crotonylation": ["K"],
    "succinylation": ["K"],
    "glutathionylation": ["C"],
    "sumoylation": ["K"],
    "s_nitrosylation": ["C"],
    "glutarylation": ["K"],
    "citrullination": ["R"],
    "o_linked_glycosylation": ["S", "T"],
    "n_linked_glycosylation": ["N"],
}

PTM_PATTERNS = {
    "phosphorylation": r"phospho|Phospho",
    "acetylation": r"acety|Acety",
    "ubiquitination": r"ubiquit|Ubiquit",
    "hydroxylation": r"hydroxy|Hydroxy",
    "gamma_carboxyglutamic_acid": r"gamma.*carboxy|Gamma.*Carboxy",
    "lys_methylation": r"lys.*methyl|Lys.*Methyl",
    "malonylation": r"malonyl|Malonyl",
    "arg_methylation": r"arg.*methyl|Arg.*Methyl",
    "crotonylation": r"crotonyl|Crotonyl",
    "succinylation": r"succinyl|Succinyl",
    "glutathionylation": r"glutathion|Glutathion",
    "sumoylation": r"sumoy|Sumoy",
    "s_nitrosylation": r"nitrosyl|Nitrosyl",
    "glutarylation": r"glutaryl|Glutaryl",
    "citrullination": r"citrullin|Citrullin",
    "o_linked_glycosylation": r"O-.*glycos|o-.*glycos",
    "n_linked_glycosylation": r"N-.*glycos|n-.*glycos",
}


def _looks_like_pdb_entry(entry: str) -> bool:
    return bool(PDB_ENTRY_RE.match(entry))


def _filter_by_type(df: pd.DataFrame, ptm_type: str) -> pd.DataFrame:
    pattern = PTM_PATTERNS[ptm_type]
    target_aa = PTM_AA_MAP[ptm_type]
    filtered = df[df["ptm"].str.contains(pattern, case=False, regex=True)]
    filtered = filtered[filtered["mod_aa"].isin(target_aa)]
    return filtered


def _sample_sites(df: pd.DataFrame, ptm_type: str, n_per_class: int, seed: int) -> pd.DataFrame:
    filtered = _filter_by_type(df, ptm_type)
    filtered = filtered[~filtered["entry"].apply(_looks_like_pdb_entry)]
    pos = filtered[filtered["label"] == 1]
    neg = filtered[filtered["label"] == 0]
    pos_sample = pos.sample(n=min(n_per_class, len(pos)), random_state=seed) if len(pos) else pos
    neg_sample = neg.sample(n=min(n_per_class, len(neg)), random_state=seed) if len(neg) else neg
    return pd.concat([pos_sample, neg_sample], ignore_index=True)


def _fetch_protein(accession: str, structures_dir: Path) -> "tuple[str, Path] | None":
    """Descarga (con cache en disco) secuencia real + estructura real de AlphaFold DB.

    Devuelve (sequence, pdb_path) o None si la proteina no tiene modelo en
    AlphaFold DB (accessions no humanas/no revisadas, etc. -- se descarta ese
    sitio de la muestra, no se aborta el script completo).
    """
    seq_cache = structures_dir / f"{accession}.seq"
    pdb_cache = structures_dir / f"{accession}.pdb"

    if seq_cache.is_file() and pdb_cache.is_file():
        return seq_cache.read_text().strip(), pdb_cache

    # Fallos de red (timeouts, blips transitorios de conectividad) degradan a
    # None (se salta esa proteina) en vez de abortar la corrida entera --
    # confirmado real: un solo ConnectionError transitorio tumbo el script
    # completo antes de este fix. 3 intentos con backoff corto porque el
    # fallo observado fue transitorio (la misma URL funciona en el reintento).
    for attempt in range(3):
        try:
            resp = requests.get(ALPHAFOLD_API_URL.format(acc=accession), timeout=30)
            if resp.status_code != 200:
                return None
            entries = resp.json()
            if not entries:
                return None
            entry = next((e for e in entries if e.get("uniprotAccession") == accession), entries[0])
            sequence = entry.get("sequence")
            pdb_url = entry.get("pdbUrl")
            if not sequence or not pdb_url:
                return None

            pdb_resp = requests.get(pdb_url, timeout=60)
            if pdb_resp.status_code != 200:
                return None

            structures_dir.mkdir(parents=True, exist_ok=True)
            seq_cache.write_text(sequence)
            pdb_cache.write_bytes(pdb_resp.content)
            return sequence, pdb_cache
        except requests.exceptions.RequestException as exc:
            print(f"  [{accession}] fallo de red (intento {attempt + 1}/3): {exc}")
            if attempt < 2:
                time.sleep(5)
    return None


def _ensure_esm_cache(protein_id: str, sequence: str, esm_checkpoint: Path, custom_esm_dir: Path) -> None:
    """Pre-popula el cache que ``predict.py``'s ``load_custom_esm`` espera.

    A diferencia de ``_deepptmpred_runner.py`` (cache con hash de secuencia,
    ver docstring de ``_esm_cache_path``), ``predict.py`` (llamado
    DIRECTAMENTE aqui, no via el runner) espera el nombre plano
    ``{protein_id}_full_esm.npz`` sin hash -- se replica ese nombre a
    proposito, en un directorio de cache PROPIO de este script (nunca el
    ``DEEPPTMPRED_CUSTOM_ESM_DIR`` de produccion), para no interferir con el
    pipeline real.
    """
    custom_esm_dir.mkdir(parents=True, exist_ok=True)
    esm_path = custom_esm_dir / f"{protein_id}_full_esm.npz"
    if esm_path.is_file():
        return
    features = _extract_esm_features(sequence, esm_checkpoint)
    np.savez_compressed(esm_path, features=features, protein_id=protein_id,
                        sequence=sequence, length=len(sequence))


def _best_threshold(labels: np.ndarray, probs: np.ndarray) -> float:
    """Umbral que maximiza el estadistico J de Youden (tpr - fpr) sobre el ROC real."""
    fpr, tpr, thresholds = roc_curve(labels, probs)
    j = tpr - fpr
    return float(thresholds[np.argmax(j)])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calibra el umbral de probabilidad de DeepPTMPred con datos reales."
    )
    parser.add_argument("--dataset", default=str(PROJECT_ROOT / "DeepPTMPred" / "data" / "ptm_data.csv.gz"))
    parser.add_argument("--train-ptm-dir", default=str(PROJECT_ROOT / "DeepPTMPred" / "pred" / "train_PTM"))
    parser.add_argument("--esm-checkpoint",
                        default=str(PROJECT_ROOT / "DeepPTMPred" / "esm" / "checkpoints" / "esm2_t33_650M_UR50D.pt"))
    parser.add_argument("--structures-cache", default=str(PROJECT_ROOT / "DeepPTMPred" / "data" / "calibration_structures"))
    parser.add_argument("--esm-cache", default=str(PROJECT_ROOT / "DeepPTMPred" / "pred" / "calibration_esm"))
    parser.add_argument("--out-dir", default=str(PROJECT_ROOT / "DeepPTMPred" / "data" / "calibration"))
    parser.add_argument("--n-per-class", type=int, default=15,
                        help="Sitios positivos y negativos a muestrear por tipo de PTM (por defecto 15+15=30).")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ptm-types", nargs="*", default=list(PTM_AA_MAP.keys()),
                        help="Subconjunto de tipos a calibrar (por defecto los 17 soportados).")
    args = parser.parse_args()

    train_ptm_dir = Path(args.train_ptm_dir)
    # _load_predict_module ya inserta train_ptm_dir en sys.path e importa
    # 'predict', ademas de aplicar el parche real de 'K' (custom_objects)
    # que predict.py::PTMPredictor.__init__ no incluye por si solo -- ver
    # docstring de _deepptmpred_runner.py::_load_predict_module.
    predict = _load_predict_module(train_ptm_dir)

    structures_dir = Path(args.structures_cache)
    esm_cache_dir = Path(args.esm_cache)
    esm_checkpoint = Path(args.esm_checkpoint)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.dataset, compression="gzip" if args.dataset.endswith(".gz") else None)
    print(f"Dataset cargado: {len(df)} filas, {df['entry'].nunique()} entries unicos.")

    summary_rows = []

    for ptm_type in args.ptm_types:
        print(f"\n=== {ptm_type} ===")
        sample = _sample_sites(df, ptm_type, args.n_per_class, args.seed)
        if sample.empty:
            print("  sin datos disponibles, se omite.")
            continue
        print(f"  muestra: {len(sample)} sitios ({(sample['label'] == 1).sum()} positivos, "
              f"{(sample['label'] == 0).sum()} negativos), {sample['entry'].nunique()} proteinas unicas.")

        config = predict.PredictConfig(ptm_type=ptm_type, project_root=str(train_ptm_dir.parent.parent))
        config.custom_esm_dir = str(esm_cache_dir)
        try:
            predictor = predict.PTMPredictor(config)
        except Exception as exc:
            print(f"  no se pudo cargar el modelo real de {ptm_type}: {exc}")
            continue

        rows = []
        for accession, group in sample.groupby("entry"):
            fetched = _fetch_protein(accession, structures_dir)
            if fetched is None:
                print(f"  [{accession}] sin modelo en AlphaFold DB, se omite ({len(group)} sitio(s)).")
                continue
            sequence, pdb_path = fetched

            positions = group["pos"].tolist()
            mismatches = [
                (pos, expected, sequence[pos - 1] if 0 < pos <= len(sequence) else None)
                for pos, expected in zip(positions, group["mod_aa"])
                if not (0 < pos <= len(sequence)) or sequence[pos - 1] != expected
            ]
            if mismatches:
                print(f"  [{accession}] {len(mismatches)} posicion(es) no coinciden con la secuencia real de "
                      f"AlphaFold DB (desfase de numeracion/isoforma), se omiten: {mismatches[:3]}")
                bad_positions = {m[0] for m in mismatches}
                group = group[~group["pos"].isin(bad_positions)]
                positions = group["pos"].tolist()
            if not positions:
                continue

            try:
                _ensure_esm_cache(accession, sequence, esm_checkpoint, esm_cache_dir)
            except Exception as exc:
                print(f"  [{accession}] fallo calculando features ESM reales: {exc}")
                continue
            # Ver docstring del modulo: forzar una PyRosettaCalculator nueva
            # por proteina, el propio predict.py no lo hace por si solo.
            predictor.data_loader.rosetta_calc = None
            try:
                pred_df = predictor.predict_ptm_sites(accession, sequence, positions, pdb_path=str(pdb_path))
            except Exception as exc:
                print(f"  [{accession}] fallo la prediccion real: {exc}")
                continue

            merged = pred_df.merge(
                group[["pos", "label"]].rename(columns={"pos": "position"}),
                on="position", how="inner",
            )
            rows.append(merged)

        if not rows:
            print(f"  ninguna proteina de la muestra pudo procesarse para {ptm_type}.")
            continue

        result_df = pd.concat(rows, ignore_index=True)
        out_path = out_dir / f"{ptm_type}_calibration.tsv"
        result_df.to_csv(out_path, sep="\t", index=False)

        n_pos = int((result_df["label"] == 1).sum())
        n_neg = int((result_df["label"] == 0).sum())
        if n_pos == 0 or n_neg == 0:
            print(f"  {len(result_df)} sitios procesados pero solo una clase presente "
                  f"({n_pos} pos / {n_neg} neg) -- no se puede medir AUROC real, se omite del resumen.")
            continue

        auroc = roc_auc_score(result_df["label"], result_df["probability"])
        threshold = _best_threshold(result_df["label"].to_numpy(), result_df["probability"].to_numpy())
        print(f"  procesados {len(result_df)}/{ len(sample)} sitios muestreados -- "
              f"AUROC real={auroc:.3f}, umbral sugerido (Youden's J)={threshold:.3f}")
        summary_rows.append({
            "ptm_type": ptm_type, "n_sites": len(result_df), "n_pos": n_pos, "n_neg": n_neg,
            "n_proteins": result_df["protein_id"].nunique(), "auroc": auroc, "suggested_threshold": threshold,
        })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = out_dir / "summary.tsv"
        summary_df.to_csv(summary_path, sep="\t", index=False)
        print(f"\nResumen escrito en {summary_path}:")
        print(summary_df.to_string(index=False))
    else:
        print("\nNingun tipo de PTM produjo un resumen con AUROC real.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
