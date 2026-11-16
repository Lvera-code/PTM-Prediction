"""Centralized pipeline configuration: paths, thresholds and external credentials.

Every adjustable parameter is resolved from environment variables with
conservative defaults, so as to allow reconfiguration without touching
source code and without exposing credentials when publishing the repository.
"""

import os
import sys
from pathlib import Path


def _env_str(key: str, default: str) -> str:
    return os.environ.get(key, default)


def _env_float(key: str, default: float) -> float:
    raw = os.environ.get(key)
    if raw is None:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def _env_int(key: str, default: int) -> int:
    raw = os.environ.get(key)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _env_bool(key: str, default: bool) -> bool:
    raw = os.environ.get(key)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


class Settings:
    """The single source of truth for all pipeline configuration."""

    # --- Phase 1 / Orchestrator: the pipeline's input and output folders ---
    FASTA_INPUT_DIR: Path = Path(_env_str("FASTA_INPUT_DIR", "inputs"))
    FASTA_OUTPUT_DIR: Path = Path(_env_str("FASTA_OUTPUT_DIR", "outputs"))

    # --- Biological validation panel: local cache derived from dbPTM ---
    # See the decision of 2026-08-09 (vault). A regenerable cache via
    # scripts/import_dbptm_panel.py, never a source, so it is not versioned (.gitignore).
    DBPTM_DATA_DIR: Path = Path(_env_str("DBPTM_DATA_DIR", "data/dbptm"))

    # --- Phase 1.5: structure extraction (PDB/mmCIF via gemmi, LOCAL) ---
    # Chain selection strategy when the input file has more
    # than one protein chain (see `src/utils/structure_parser.py`). Never
    # implicit: the chosen chain is always logged together with the reason.
    #   'longest'  -> the chain with the most residues in its polymer is chosen.
    #   'explicit' -> PDB_EXPLICIT_CHAIN_ID is used (mandatory in that case).
    PDB_CHAIN_SELECTION_STRATEGY: str = _env_str("PDB_CHAIN_SELECTION_STRATEGY", "longest")
    PDB_EXPLICIT_CHAIN_ID: str = _env_str("PDB_EXPLICIT_CHAIN_ID", "")

    # --- Phase 2: DeepMVP (the only FASTA path engine, engine 1/2 of the PDB path) ---
    # Verified by reading github.com/bzhanglab/DeepMVP directly (README.md,
    # DeepMVP.py, lib/PTModels.py, lib/Metrics.py) on 2026-07-27. A real,
    # packaged repo (requirements.txt + environment.yml), GPL-3.0.
    DEEPMVP_HOME: Path = Path(_env_str("DEEPMVP_HOME", "DeepMVP"))
    DEEPMVP_SCRIPT_NAME: str = _env_str("DEEPMVP_SCRIPT_NAME", "DeepMVP.py")
    # Python 3.7.10 + TensorFlow==2.4.2 (see the repo's environment.yml): an
    # old stack, for which a dedicated venv is recommended (the same pattern as BepiPred-3.0 in
    # BCell-Epitope-Prediction). By default it uses the same interpreter that
    # runs pipeline.py; point this at a dedicated venv in production.
    DEEPMVP_PYTHON_BIN: str = _env_str("DEEPMVP_PYTHON_BIN", sys.executable)
    # Pre-trained weights NOT included in the repo (confirmed in the README):
    # a manual download from DEEPMVP_DOWNLOAD_URL (one .tar.gz per version),
    # to be unpacked inside this folder. It contains one subdirectory per PTM
    # type (acetylation_k/, glycosylation_n/, methylation_k/, methylation_r/,
    # phosphorylation_st/, phosphorylation_y/, sumoylation_k/, ubiquitination_k/),
    # each with its own model.json + an ensemble of .h5 models, confirmed by
    # reading lib/PTModels.py::ptm_prediction_for_multiple_ptms.
    DEEPMVP_MODEL_DIR: Path = Path(_env_str("DEEPMVP_MODEL_DIR", "DeepMVP/models"))
    DEEPMVP_DOWNLOAD_URL: str = "https://deepmvp.ptmax.org/"
    DEEPMVP_TIMEOUT_SECONDS: int = _env_int("DEEPMVP_TIMEOUT_SECONDS", 1800)
    # Task 2 = "PTM site prediction" in DeepMVP's CLI (task 1 is
    # mutation impact prediction, not used here). Confirmed by reading
    # DeepMVP.py directly, and not configurable via Settings: it is a fixed value
    # of the contract with the real CLI, not an adjustable pipeline parameter.
    DEEPMVP_TASK_PTM_SITE_PREDICTION: str = "2"
    # Confidence filter (Phase 3, the core, decision 2026-07-27): DeepMVP does
    # not publish a fixed probability cutoff per PTM type. Instead,
    # each weights folder ships its own validation 'site_prediction.tsv',
    # and the 'fpr' column it returns per row is the real FPR that
    # that specific model would have if that probability were used as the cutoff
    # (confirmed by reading lib/Metrics.py::add_confidence_metrics/calc_FPR).
    # It is a BETTER calibrated threshold than a fixed probability score: the
    # Phase 3 filter uses 'fpr <= DEEPMVP_MAX_FPR', not a threshold on
    # 'y_pred'. The raw score (y_pred) is always preserved unchanged.
    DEEPMVP_MAX_FPR: float = _env_float("DEEPMVP_MAX_FPR", 0.05)

    # --- Phase 2: DeepPTMPred (engine 2/2 of the consensus, PDB path only) ---
    # Verified by reading github.com/kuikui-wang/DeepPTMPred directly
    # (README.md, pred/train_PTM/predict.py, pred/train_PTM/e2_single_data.py,
    # pred/train_PTM/environment.yml) on 2026-07-27. A real repo, with weights
    # included (.h5 per type, ~19MB each). Licence RESOLVED
    # 2026-07-29: the repo has no declared LICENSE, but the paper's
    # corresponding author (CC BY-NC 4.0) confirmed
    # directly by email that the repo's code follows those same
    # terms, so non-commercial use (research/TFG, Scipion/CNB) is covered.
    #
    # AN IMPORTANT FINDING: unlike DeepMVP, neither predict.py nor
    # e2_single_data.py has a real CLI; ptm_type/pdb_path/protein_id/the ESM
    # checkpoint path are hardcoded inside their
    # `if __name__ == "__main__":` block. That is why this project does NOT invoke those
    # scripts directly by subprocess (unlike DeepMVP): it uses its own
    # runner (`src/engines/_deepptmpred_runner.py`) which imports their
    # parameterizable classes (PredictConfig, PTMPredictor) and REIMPLEMENTS the
    # ESM-2 feature extraction, because e2_single_data.py::extract_full_sequence_esm
    # has a confirmed real bug: it redefines 'custom_checkpoint_path' as a
    # LOCAL variable with a hardcoded absolute AutoDL path
    # (/root/autodl-tmp/...), ignoring any value passed or set at module level.
    DEEPPTMPRED_HOME: Path = Path(_env_str("DEEPPTMPRED_HOME", "DeepPTMPred"))
    DEEPPTMPRED_TRAIN_PTM_DIR: Path = Path(
        _env_str("DEEPPTMPRED_TRAIN_PTM_DIR", "DeepPTMPred/pred/train_PTM")
    )
    DEEPPTMPRED_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "DEEPPTMPRED_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_deepptmpred_runner.py"),
        )
    )
    # Python 3.10 + TensorFlow==2.15 + PyTorch 2.0 + PyRosetta + fair-esm (see
    # the repo's environment.yml): a heavy stack, DIFFERENT from DeepMVP's
    # (Python 3.7 + TF 2.4), which requires its own dedicated venv/conda, never
    # shared with DEEPMVP_PYTHON_BIN. tensorflow-addons (used by
    # predict.py for its loss function) has been archived/deprecated upstream
    # since 2024. RISK DISCARDED 2026-07-27: genuinely installed and imported
    # in an isolated venv against TF 2.15, without error (only an
    # expected deprecation UserWarning).
    DEEPPTMPRED_PYTHON_BIN: str = _env_str("DEEPPTMPRED_PYTHON_BIN", sys.executable)
    # ESM-2 checkpoint (fair-esm, ~2.5GB), a separate download (not included in
    # the repo, unlike the PTM .h5 weights, which are included).
    # IMPORTANT (verified 2026-07-27 by reading esm/pretrained.py directly):
    # fair-esm requires a COMPANION file next to the main checkpoint,
    # '<checkpoint>-contact-regression.pt' (same directory, same base
    # name). Its internal heuristic (_has_regression_weights) does not exclude
    # the esm2_* models, so it ALWAYS looks for it except for the esm1v/
    # esm_if/270K/500K variants. If it is missing, torch.load() fails with a
    # LOCAL FileNotFoundError (it is not a network attempt; see the docstring of
    # _deepptmpred_runner.py::_extract_esm_features for the complete
    # detail of why this IS 100% local despite this extra file).
    # When downloading the checkpoint, verify that this companion file
    # is also available (the same origin as the main checkpoint).
    DEEPPTMPRED_ESM_CHECKPOINT: Path = Path(
        _env_str("DEEPPTMPRED_ESM_CHECKPOINT", "DeepPTMPred/esm/checkpoints/esm2_t33_650M_UR50D.pt")
    )
    # Per-accession cache of ESM features (.npz), reused between runs
    # and between PTM types of the same accession (it avoids recomputing the
    # complete embedding 17 times, once per PTM type).
    DEEPPTMPRED_CUSTOM_ESM_DIR: Path = Path(
        _env_str("DEEPPTMPRED_CUSTOM_ESM_DIR", "DeepPTMPred/pred/custom_esm")
    )
    DEEPPTMPRED_TIMEOUT_SECONDS: int = _env_int("DEEPPTMPRED_TIMEOUT_SECONDS", 3600)
    # The 17 supported PTM types, confirmed by reading
    # predict.py::PredictConfig.ptm_aa_map (the repo's exact names, not
    # translated): unlike DeepMVP (task=2 predicts all the types in
    # a single invocation), DeepPTMPred predicts ONE type per invocation, so the
    # engine must invoke the runner once PER TYPE, per accession.
    # Confidence filter (Phase 3, the core): unlike DeepMVP,
    # DeepPTMPred exposes no 'fpr'-style calibration mechanism (neither in
    # predict.py nor in its README), only a raw 0-1 probability. Its
    # own script uses a hardcoded cutoff of 0.5 without documenting where it
    # comes from (not calibrated against any published validation set, unlike
    # DeepMVP's 'fpr').
    #
    # CALIBRATED 2026-08-01 (replacing the 2026-07-31 calibration, invalidated
    # by the plDDT bug, see below): real per-type thresholds, computed
    # by `scripts/generate_deepptmpred_calibration.py` +
    # `_rebuild_calibration_summary.py` over
    # `github.com/meilerlab/PTMPrediction/data/ptm_data.csv.gz`
    # (n=75 proteins/class), WITH two train/inference fixes already applied
    # in `_deepptmpred_runner.py::_load_predict_module`: (1) phi/psi=0.0
    # forced at inference (2026-07-31, matching the training
    # distribution) and (2) local_plDDT recomputed from the real CA B-factor
    # of the PyRosetta pose instead of `predict.py`'s `plDDT_values` array
    # (badly indexed per atom rather than per residue, 2026-08-01). The real
    # rule (bug-fixed with respect to the paper, not literal):
    # FPR<=0.18 TPR-max for phosphorylation, FPR<=0.20 TPR-max
    # for the rest; falling back to 0.5 if no operating point meets the
    # maximum FPR in the sample.
    # Exact source: `DeepPTMPred/data/calibration/summary.tsv` (gitignored,
    # regenerable with the two scripts above). With the plDDT fix, 16/17
    # types land in the solid-to-excellent range (AUROC 0.77-0.99), including
    # the 4 types that were previously "mediocre" (crotonylation 0.61->0.82,
    # citrullination 0.66->0.78, glutarylation 0.67->0.77, s_nitrosylation
    # 0.68->0.77), which confirms that the plDDT bug (not a model ceiling)
    # was the real cause. `n_linked_glycosylation` remains at AUROC 0.51 (worse than
    # chance, a real limit of the published weights), which is why
    # it is EXCLUDED from the consensus in `ptm_annotation.py` despite having a
    # calibrated threshold here (the threshold is not the problem, the model's
    # discriminative power is).
    DEEPPTMPRED_CALIBRATED_THRESHOLDS: dict = {
        "acetylation": 0.6350621,
        "arg_methylation": 0.34068727,
        "citrullination": 0.36854228,
        "crotonylation": 0.86312497,
        "gamma_carboxyglutamic_acid": 0.24807824,
        "glutarylation": 0.4747097,
        "glutathionylation": 0.466462,
        "hydroxylation": 0.35899624,
        "lys_methylation": 0.43064716,
        "malonylation": 0.41699925,
        "n_linked_glycosylation": 0.99802846,
        "o_linked_glycosylation": 0.2619363,
        "phosphorylation": 0.24020174,
        "s_nitrosylation": 0.5140331,
        "succinylation": 0.50403893,
        "sumoylation": 0.37326753,
        "ubiquitination": 0.5321394,
    }
    # Fallback if DEEPPTMPRED_CALIBRATED_THRESHOLDS does not have the type (which
    # should not happen with the 17 types of DEEPPTMPRED_PTM_TYPES, but it
    # keeps the pipeline functional if a new type is added without
    # recalibrating). It is also a manual override via env var for occasional
    # experimentation without editing code.
    DEEPPTMPRED_MIN_PROBABILITY: float = _env_float("DEEPPTMPRED_MIN_PROBABILITY", 0.5)

    @classmethod
    def deepptmpred_threshold_for(cls, ptm_type: str) -> float:
        """The calibrated threshold for ``ptm_type``, or the generic fallback if absent."""
        return cls.DEEPPTMPRED_CALIBRATED_THRESHOLDS.get(ptm_type, cls.DEEPPTMPRED_MIN_PROBABILITY)

    DEEPPTMPRED_PTM_TYPES: tuple = (
        "phosphorylation", "acetylation", "ubiquitination", "hydroxylation",
        "gamma_carboxyglutamic_acid", "lys_methylation", "malonylation",
        "arg_methylation", "crotonylation", "succinylation", "glutathionylation",
        "sumoylation", "s_nitrosylation", "glutarylation", "citrullination",
        "o_linked_glycosylation", "n_linked_glycosylation",
    )

    # --- Optional corroboration: MeToken (PDB path only, NEVER a consensus
    # engine; Decision 2 remains paused) ---
    # Verified by reading github.com/A4Bio/MeToken (ICLR 2025) on 2026-08-01: the
    # published checkpoint is a TYPE CLASSIFIER on ALREADY KNOWN sites
    # (confirmed in model_interface.py:40, where the "no-PTM" class is masked
    # when with_null_ptm=0, as the published checkpoint comes, and
    # verified empirically: it predicts types with high confidence also at
    # positions without a real PTM). That is why it NEVER decides pasa_umbral/consenso,
    # the same purely informative pattern as the secretory pathway
    # (src/structural/uniprot_localization_client.py): it only corroborates the TYPE
    # at sites that the consensus has ALREADY accepted (pasa_umbral=true); see
    # src/engines/ptm_annotation.py::annotate_pdb_path and
    # src/engines/metoken_engine.py.
    METOKEN_ENABLED: bool = _env_bool("METOKEN_ENABLED", True)
    METOKEN_HOME: Path = Path(_env_str("METOKEN_HOME", "MeToken"))
    METOKEN_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "METOKEN_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_metoken_runner.py"),
        )
    )
    # Python 3.10 + torch (CPU) + torch_scatter (compiled from source, as there is no
    # prebuilt wheel for this torch/CPU combination) + transformers +
    # biopython>=1.80 (see bug 1 patched in _metoken_runner.py): a dedicated
    # venv, never shared with DEEPMVP_PYTHON_BIN/DEEPPTMPRED_PYTHON_BIN.
    METOKEN_PYTHON_BIN: str = _env_str("METOKEN_PYTHON_BIN", sys.executable)
    # Real weights (release 1.0, pretrained_model.zip, ~88MB, verified as
    # downloadable on 2026-08-01): checkpoint.ckpt is a flat state_dict
    # (confirmed with torch.load + model.load_state_dict -> "<All keys
    # matched successfully>"), not a complete pytorch-lightning checkpoint
    # (that is lightning_checkpoint.ckpt, not used here).
    METOKEN_CHECKPOINT: Path = Path(
        _env_str("METOKEN_CHECKPOINT", "MeToken/pretrained_model/checkpoint.ckpt")
    )
    METOKEN_TIMEOUT_SECONDS: int = _env_int("METOKEN_TIMEOUT_SECONDS", 600)

    # --- Kinase Library (informative corroboration of kinase specificity,
    # both paths; biological coherence analysis 2026-08-07, point 5) ---
    # Johnson et al. 2023 Nature (303 Ser/Thr kinases) + Yaron-Barir et al. 2024
    # Nature (the complete Tyr kinome), packaged in 'kinase-library' (PyPI,
    # TheKinaseLibrary/kinase-library, licence CC-BY-NC-SA-3.0, the same
    # non-commercial family already accepted for EMNGly). It only needs the complete
    # sequence plus the 1-based position of the phosphoacceptor, so it applies to both paths.
    # It NEVER decides pasa_umbral/consenso, the same non-deciding pattern as
    # MeToken.
    KINASE_LIBRARY_ENABLED: bool = _env_bool("KINASE_LIBRARY_ENABLED", True)
    # A DEDICATED conda environment (numpy~=1.26.4/pandas~=2.2.3 pinned by the
    # package itself, incompatible with the versions pinned by THIS venv, see
    # requirements.txt), never the main venv. Verified 2026-08-07 with a
    # real run (kl.Substrate('PSVEPPLsQETFSDL').predict() over p53 S33,
    # top hit ATM, matching the real literature on the DNA damage response).
    KINASE_LIBRARY_PYTHON_BIN: str = _env_str(
        "KINASE_LIBRARY_PYTHON_BIN",
        "/home/enzo/miniconda3/envs/kinase_library/bin/python",
    )
    KINASE_LIBRARY_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "KINASE_LIBRARY_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_kinase_library_runner.py"),
        )
    )
    # CORRECTED 2026-08-07 (post-implementation robustness audit): the
    # original assumption of "much faster than MeToken, without a GPU" was
    # CORRECT per site but IRRELEVANT, because the real cost is
    # per BATCH, not per individual site, and the batch is not small. Measured
    # 2026-08-07 against Tau's complete sequence (758 aa): ~5s/site, so a
    # batch of 30 sites already took ~150s, MORE than the original timeout of 120s.
    # Real proteins of the validation panel itself accept far more
    # phosphorylation sites than that: the complete Tau has 115 sites
    # accepted by the consensus (see outputs/AF-P10636-F1-model_v4_ptm_sites.csv),
    # which at that rate would take ~10 min. Confirmed live the same day:
    # a real run of the validation panel (p53, 33 sites) took more than
    # 120s with the original default, so the old timeout guaranteed silent
    # failure (degradation to None) on any protein with more than ~25
    # accepted sites, defeating the engine's purpose. Generous by
    # default, with 900s of margin beyond the peak already seen in production,
    # because it is the number of sites per protein (not per run) that
    # scales here.
    KINASE_LIBRARY_TIMEOUT_SECONDS: int = _env_int("KINASE_LIBRARY_TIMEOUT_SECONDS", 1800)

    # --- EMNGly (the real consensus engine for 'n_linked_glycosylation', PDB path
    # only; it replaces CoNglyPred, decision 2026-08-06) ---
    # CoNglyPred (github.com/whm242446/CoNglyPred, the original Decision 2
    # candidate) was confirmed
    # DEFINITIVELY unavailable 2026-08-06 (re-verified via the GitHub API: 0
    # releases/tags/issues, last push 2024-08-15, ZERO .pth/.pt/.pkl files
    # in the whole tree despite its own README instructing one to load
    # 'best.pth'), and the email of 2026-08-01 to the corresponding author
    # never produced a reply with the weights.
    #
    # EMNGly (github.com/StellaHxy/EMNgly, Hou et al., Bioinformatics
    # 39(11):btad650, 2023) is the replacement: verified exhaustively
    # 2026-08-06 (Opus subagent + our own re-verification by cloning the real
    # repo) that it DOES have real, downloadable weights, unlike
    # CoNglyPred. Architecture: ESM-1b (site_emb 1280 + local_emb 1280,
    # ``model/get_esm_embedding.py``) + MIF (Microsoft's Masked Inverse
    # Folding, ``model/MIF/``, a real STRUCTURAL embedding (256-dim) over
    # the PDB's N/CA/C backbone coordinates), concatenated (2816-dim) and
    # classified with an SVM (``sklearn.svm.SVC``, rbf kernel,
    # probability=True). It preserves the design property already decided when
    # discarding MTPrompt-PTM as a replacement for DeepPTMPred: the second engine
    # of this type must consume real 3D structure, not sequence alone.
    #
    # Weights verified at BYTE level, not from the README (2026-08-06):
    # - MIF (``mif.pt``, 13.8MB): bundled directly in the repo, with no
    #   separate download (official alternative source: Microsoft's Zenodo
    #   10.5281/zenodo.6573779).
    # - SVM (``N-GlyDE.pickle``, 36MB): a public Google Drive folder
    #   (id ``1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw``), downloaded with a
    #   real HTTP Range request without authentication; the first bytes
    #   confirm pickle protocol 4, first object
    #   ``sklearn.svm._classes.SVC``, ``n_features_in_=2816`` (matching
    #   EXACTLY 1280+1280+256), ``_sklearn_version=1.1.1``.
    #   ``N-GlyDE.pickle`` was chosen (not ``N-GlyAltas_classifier.pkl``, also
    #   available in the same folder) because it is trained over the
    #   benchmark with NEGATIVES RESTRICTED TO THE N-X-[S/T] SEQUON, the
    #   correct benchmark regime for this role: a model trained
    #   against unrestricted negatives (dbPTM/MusiteDeep) can be right
    #   merely by learning the N-X-[S/T] motif, which this pipeline ALREADY applies
    #   via ``_nglyco_window``, so it would not contribute a real second opinion.
    #
    # It has NO declared LICENSE in any repo (the paper is CC-BY but that
    # does NOT cover the code). Unlike CoNglyPred, this is NOT
    # blocking: the weights are already downloadable, and the email pending to
    # the corresponding authors (ICT-CAS) is only to leave the licence
    # record clean (the same pattern that worked with
    # DeepPTMPred), not a critical dependency.
    EMNGLY_ENABLED: bool = _env_bool("EMNGLY_ENABLED", True)
    EMNGLY_HOME: Path = Path(_env_str("EMNGLY_HOME", "EMNgly"))
    EMNGLY_RUNNER_SCRIPT: Path = Path(
        _env_str(
            "EMNGLY_RUNNER_SCRIPT",
            str(Path(__file__).resolve().parent.parent / "engines" / "_emngly_runner.py"),
        )
    )
    # Python 3.x + fair-esm + torch + scikit-learn (a dedicated venv is
    # recommended, 'emngly', never shared with other engines). A REAL
    # FINDING 2026-08-06: the SVM pickles were generated with
    # scikit-learn==1.1.1, but the repo's 'environment.yml' pins 1.5.1,
    # and loading with a version different from the training one triggers
    # 'InconsistentVersionWarning' (it probably works, as SVC is stable,
    # but it is not guaranteed), so pin scikit-learn==1.1.1 in the venv,
    # not what the repo's environment.yml asks for.
    EMNGLY_PYTHON_BIN: str = _env_str("EMNGLY_PYTHON_BIN", sys.executable)
    # ESM-1b (esm1b_t33_650M_UR50S, ~7.4GB plus the contact regression
    # companion, the same pattern as Settings.DEEPPTMPRED_ESM_CHECKPOINT
    # above): a separate download, NOT included in the repo. Unlike
    # the original script (``model/get_esm_embedding.py``, which uses
    # ``torch.hub.load`` against ``dl.fbaipublicfiles.com`` on every run,
    # violating this project's local-only policy), this project's
    # runner ALWAYS requires a local ``.pt`` path, the same mechanism already
    # verified in ``_deepptmpred_runner.py`` (fair-esm enters the
    # local branch of ``pretrained.load_model_and_alphabet`` as soon as the name
    # ends in '.pt', and never calls the network branch).
    EMNGLY_ESM_CHECKPOINT: Path = Path(
        _env_str("EMNGLY_ESM_CHECKPOINT", "EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S.pt")
    )
    # MIF weights, bundled in the repo (see above), with no separate download.
    EMNGLY_MIF_WEIGHTS: Path = Path(_env_str("EMNGLY_MIF_WEIGHTS", "EMNgly/model/MIF/weights/mif.pt"))
    # An SVM trained on N-GlyDE (negatives restricted to the sequon, see
    # above): a manual download from Google Drive, see the README.md
    # Installation section.
    EMNGLY_SVM_CHECKPOINT: Path = Path(
        _env_str("EMNGLY_SVM_CHECKPOINT", "EMNgly/checkpoints/N-GlyDE.pickle")
    )
    # Per-accession cache of ESM-1b/MIF embeddings, the same pattern as
    # DEEPPTMPRED_CUSTOM_ESM_DIR (it avoids recomputing the complete
    # sequence/structure embedding on every run).
    EMNGLY_CACHE_DIR: Path = Path(_env_str("EMNGLY_CACHE_DIR", "EMNgly/cache"))
    EMNGLY_TIMEOUT_SECONDS: int = _env_int("EMNGLY_TIMEOUT_SECONDS", 1800)
    # Threshold 0.5 (decision 2026-08-06): the EMNGly paper does not publish a
    # fixed probability cutoff, only AUC/MCC at whole-dataset level,
    # and 0.5 is the native decision boundary of an SVC(probability=True) without
    # calibration against a specific operating threshold. It is NO longer provisional:
    # the 2 go/no-go checks PASSED
    # 2026-08-07, with the structure_emb alignment verified against real
    # GlyGen sites in a PDB with gaps, and MCC=0.8197 reproduced (and
    # exceeded) against the published 0.736, at this very threshold.
    EMNGLY_MIN_PROBABILITY: float = _env_float("EMNGLY_MIN_PROBABILITY", 0.5)

    # --- N-glycosylation consensus (PDB path): the 2-engine fusion rule ---
    # Decision 2026-08-06: it promotes EMNGly (together with StackGlyEmbed at the
    # time, removed from the project 2026-08-10 by project decision, as it
    # depended on a sibling project's venv, real friction for the
    # Scipion integration) to a real consensus engine for
    # 'n_linked_glycosylation'/'glycosylation_n' on the PDB path (the only place
    # where EMNGly can run, as it demands a real PDB). The reason: DeepPTMPred
    # is confirmed dead for this exact type (AUROC~=0.51,
    # CONSENSUS_EXCLUDED_TYPES in ptm_annotation.py), so without a second real
    # engine this type was left with DeepMVP alone, with no
    # consensus possible.
    #
    # The rule (provisional, a project decision if it is adjusted after real use):
    # of the engines that DID manage to evaluate the position (DeepMVP always;
    # EMNGly degrades to absent without raising if it is not installed or the
    # subprocess fails, see emngly_engine.py), 'pasa_umbral' = at least 1
    # engine passes its own threshold (generalizing the 2-engine OR rule already
    # used for the other types); 'consenso' = at least
    # NGLYCO_CONSENSUS_MIN_ENGINES (default 2) pass, which in practice
    # requires BOTH engines (DeepMVP and EMNGly) to pass. If only DeepMVP
    # managed to evaluate the position (EMNGly degraded), 'consenso' is ALWAYS
    # False, as there is no way to reach 2 with a single available engine,
    # the same behaviour (single engine, no consensus) that existed before
    # this improvement.
    NGLYCO_CONSENSUS_MIN_ENGINES: int = _env_int("NGLYCO_CONSENSUS_MIN_ENGINES", 2)

    # --- Informative biological coherence warnings (analysis 2026-08-07,
    # see README.md "Scope and interpretation"): both NEVER decide
    # pasa_umbral/consenso, the same pattern as MeToken.
    # Unlike it, they do not depend on pdb_path:
    # they run on both paths whenever they are enabled (they are
    # lightweight checks, one HTTP query to UniProt per accession, or pure
    # pandas logic, with no external installation to degrade). Both
    # default to True, and are toggleable in case a
    # deterministic/offline run is needed (tests, network unavailable).
    SECRETORY_PATHWAY_CHECK_ENABLED: bool = _env_bool("SECRETORY_PATHWAY_CHECK_ENABLED", True)
    PTM_CROSSTALK_CHECK_ENABLED: bool = _env_bool("PTM_CROSSTALK_CHECK_ENABLED", True)

    @classmethod
    def ensure_dirs(cls) -> None:
        """Creates all the configured input/output folders if they do not exist."""
        cls.FASTA_INPUT_DIR.mkdir(parents=True, exist_ok=True)
        cls.FASTA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        cls.DBPTM_DATA_DIR.mkdir(parents=True, exist_ok=True)
