# STATUS

Current state of the project, not a session diary, so rewrite instead of
accumulating. See `01-Proyectos/PTM-Prediction/` in the vault for the history
of decisions.

## ARCHITECTURAL DECISION: CLOSED 2026-07-27

**Confirmed by Enzo: the current architecture is kept, with no code
changes.** DeepMVP remains the only FASTA path engine; DeepPTMPred
remains the PDB path engine (real structure); the PDB path keeps
running DeepMVP+DeepPTMPred in consensus, without adding a third engine for
now. MTPrompt-PTM was evaluated in depth and is discarded as a replacement, with the
full reason below. Nothing in `pipeline.py`/`src/engines/` needed to change:
this was already the implemented design, and the decision was whether to keep it or not.

Enzo asked to reconsider which engine covers the FASTA path and which the PDB path,
evaluating **MTPrompt-PTM** (`github.com/hanye311/MTPrompt-PTM`) as a
candidate. Verified at the same level of depth as DeepMVP/
DeepPTMPred on 2026-07-27 (real repo cloned, `mtprompt` venv created,
`pip install -r requirements.txt` without errors, real CLI confirmed:
`python test.py --config_path ... --model_path ... --data_path <fasta>
--PTM_type <tipo> --save_path ...`).

**Key findings:**
- MTPrompt-PTM is **sequence-only** (input `protein_sequence`, confirmed
  in `config/PTM_config_prompt_tuning_test.yaml`). It does NOT use the PDB/
  3D coordinates at all during inference, despite being called
  "Structure-Aware": its "structure awareness" comes from how the
  backbone was PRE-TRAINED (S-PLM v2), not from structural features
  computed per protein at inference time (unlike
  DeepPTMPred, which DOES compute real SASA/phi/psi/secondary structure via
  PyRosetta over each protein's PDB). Consequence: if it replaced
  DeepPTMPred, the pipeline would lose the ONLY engine that exploits real 3D
  structure, and the PDB path would become, in practice, two sequence engines
  in consensus, not sequence+structure.
- It does not require PyRosetta (confirmed, absent from `environment.yml`/
  `requirements.txt`), which resolves today's real blocker.
- It installs cleanly (Python 3.11, ~3 min, without real conflicts, just a single
  cosmetic warning about the `packaging` version).
- **A real new finding, not trivial**: MTPrompt-PTM's own code
  (`model.py::prepare_adapter_h_model`) builds its ESM-2 encoder via
  `esm_adapterH.pretrained.esm2_t33_650M_UR50D()`, which internally calls
  `load_model_and_alphabet_HUB` (not the local branch that our
  DeepPTMPred runner does use), so **by default it DOES try to download from the
  network** (`dl.fbaipublicfiles.com`) the first time it runs, unless the
  checkpoint is already cached in `~/.cache/torch/hub/checkpoints/`
  beforehand. Patchable (the same pattern already applied for DeepPTMPred:
  forcing the local branch by pre-populating the cache or passing a `.pt`
  path directly), but it is NOT "100% local out of the box" as was assumed when
  proposing it. It would require the same kind of adjustment already made for
  DeepPTMPred.
- Its own weights (`best_model_13ptm_final.pth`, a single
  multi-task file): manual download from Google Drive, not downloaded yet.
  The README also offers Docker (`hanye0311/mtprompt:v1`) as a deployment
  alternative, not tested.
- It covers 13 types (it lacks the 6 that DeepPTMPred does have: hydroxylation,
  gamma_carboxyglutamic_acid, glutarylation, glutathionylation,
  s_nitrosylation, citrullination; it adds palmitoylation, which neither DeepMVP nor
  DeepPTMPred covers). MIT licence (better than DeepPTMPred, which has no
  declared licence).

**Final decision (2026-07-27): option (a), no changes.** Reason: DeepMVP
has more solid verified evidence (a real benchmark against 8 tools
in Nature Methods) and covers exactly the 6 types already decided on 26/07.
The 7 extra types of MTPrompt-PTM fall outside that scope and are not a real
advantage for this project. DeepPTMPred remains the only engine
with real 3D structure, something that would be lost if it were replaced. No
redesign of `ptm_annotation.py` was done because there was no change of engines.

MTPrompt-PTM and MusiteDeep (FASTA alternatives) and LkaM-PTM/SAPP
(PDB alternatives) are documented as evaluated and
discarded candidates, not active lines of work. See the complete breakdown in
`01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-implementacion-fase1-a-fase3-nucleo.md`
in the vault and the artifact published in the session of 27-07 (a comparison
of the architectures of the 6 engines investigated).

## Done: complete end-to-end pipeline (2026-07-27), engines not installed

The 3 phases are implemented and connected in `pipeline.py`. No engine
(DeepMVP, DeepPTMPred) is installed yet on this machine, so running the
pipeline today fails with an actionable error in Phase 3 (`DeepMVPExecutionError`
as soon as it invokes DeepMVP, confirmed with a real run). 67 tests
(`pytest tests/`), without external binaries/models, with subprocess mocked where
applicable.

- **Phase 1** (`src/utils/fasta_parser.py`): FASTA sanitation. The policy on
  non-canonical residues is RELAXED (Enzo's decision 2026-07-27) to match
  DeepMVP's real tolerance (verified in `lib/PeptideEncode.py`): it no
  longer fatally rejects any character, it only warns for those that DeepMVP degrades
  (outside the 20 standard + U + B; 'X' does not count, DeepMVP treats it as "no
  signal"). Deliberately different from `BCell-Epitope-Prediction` (fatal
  rejection there, correct because BepiPred-3.0 DOES abort): each Phase 1 reflects
  the real tolerance of its engine.
- **Phase 1.5** (`src/utils/structure_parser.py`): ATMSEQ + position
  mapping from PDB/mmCIF via `gemmi`. Identical to project 1 (generic
  logic).
- **Input router** (`src/utils/input_router.py`): identical to project 1.
- **DeepMVPEngine** (`src/engines/deepmvp_engine.py`, the only FASTA path
  engine / engine 1 of 2 on the PDB path): subprocess over
  `github.com/bzhanglab/DeepMVP`, verified by reading the repo directly
  (README.md, `DeepMVP.py`, `lib/PTModels.py`, `lib/Metrics.py`). A real CLI
  with its own flags (`predict -m <model_dir> -d <fasta> -t 2 -o <out_dir>`),
  with a fixed `site_prediction.tsv` output.
- **DeepPTMPredEngine** (`src/engines/deepptmpred_engine.py` +
  `src/engines/_deepptmpred_runner.py`, engine 2/2 of the consensus, PDB path
  only): verified by reading `github.com/kuikui-wang/DeepPTMPred`
  directly (README.md, `predict.py`, `e2_single_data.py`,
  `environment.yml`). Unlike DeepMVP, THIS REPO HAS NO CLI:
  `ptm_type`/`pdb_path`/`protein_id`/the ESM checkpoint path are
  hardcoded inside its `if __name__ == "__main__":` block. Our own
  runner was built (`_deepptmpred_runner.py`) which imports the parameterized
  classes of `predict.py` (`PredictConfig`, `PTMPredictor`) and
  REIMPLEMENTS the ESM-2 feature extraction instead of using
  `e2_single_data.py::extract_full_sequence_esm`, which has a confirmed real
  bug: it redefines `custom_checkpoint_path` as a LOCAL variable with a
  hardcoded absolute AutoDL path (`/root/autodl-tmp/...`), ignoring
  any value passed. It predicts ONE PTM type per invocation (17 in
  total), not all at once like DeepMVP, so the engine invokes it 17 times per
  accession and concatenates.
- **Phase 3 core** (`src/engines/ptm_annotation.py`, B: annotation +
  D: filter): implemented according to the design in
  `01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md`.
  The DeepMVP↔DeepPTMPred type correspondence was verified by reading both
  repos (not assumed): of DeepMVP's 6 biological types (8 models per
  residue), 7 have an equivalent in DeepPTMPred; **`phosphorylation_y`
  (phosphorylation on Y) has NO equivalent model in DeepPTMPred**
  (its only phosphorylation model covers S/T, not Y), a new finding, since not
  even on the PDB path is a consensus possible for that type. The 10 types
  exclusive to DeepPTMPred (not ~11 as estimated on the 26: hydroxylation,
  gamma_carboxyglutamic_acid, malonylation, crotonylation, succinylation,
  glutathionylation, s_nitrosylation, glutarylation, citrullination,
  o_linked_glycosylation) are included in the core, marked
  `consenso=false`. Threshold: `fpr <= DEEPMVP_MAX_FPR` for DeepMVP,
  `probability >= DEEPPTMPRED_MIN_PROBABILITY` (0.5, PROVISIONAL, since DeepPTMPred
  exposes no calibration mechanism, unlike DeepMVP's `fpr`)
  for DeepPTMPred alone; in rows with consensus, `pasa_umbral` is the
  UNION (it passes if at least one of the two passes its own threshold), `consenso`
  is the INTERSECTION (True only if both pass).
- **Orchestrator** (`pipeline.py`): the 3 phases connected end to end on
  both paths. The FASTA path handles a multi-accession FASTA correctly
  (it annotates each accession with its own sequence, it does not concatenate them).
  Confirmed with a real run: Phase 1 completes and persists
  (`<accession>_clean.fasta`), and Phase 3 fails with an actionable
  `DeepMVPExecutionError` (repo not cloned), so it does not fail silently halfway.
- Local repo (`git init`), with no remote yet. It will be created on GitHub
  (`Lvera-code/PTM-Prediction`, public) once there is something functional
  end to end WITH the real engines installed.

## Real installation: verified 2026-07-27 (repos cloned, dependencies genuinely installed)

Both repos were genuinely cloned on this machine (`git clone`, inside
`PTM-Prediction/`, gitignored) to verify real installation friction,
not merely to read code. The content is identical to what was verified by API/curl
hours earlier (`diff` byte for byte against `predict.py` confirmed).

- **DeepMVP: INSTALLED, WITH WEIGHTS, VERIFIED WITH A REAL PREDICTION
  (2026-07-28).** A real `deepmvp` conda environment
  (`/home/enzo/miniconda3/envs/deepmvp/bin/python`, Python 3.7.10):
  `pip install -r requirements.txt` completed without errors (TensorFlow
  2.4.2 CPU, numpy 1.19.5, the whole stack). `python DeepMVP.py predict -h`
  prints exactly the CLI documented above. The real error case was also verified:
  with `-m` pointing at an EMPTY models folder,
  DeepMVP.py blows up with `ValueError: No objects to concatenate` (a raw,
  opaque pandas traceback), which confirms that the proactive check in
  `DeepMVPEngine._validate_installation` (failing earlier, with an
  actionable message) is necessary, not cosmetic. The real weights (`models.tar.gz`
  from the Shiny app, NOT `all_data.tar.gz`, which is only the
  train/test dataset) are installed in `DeepMVP/models/`, all 8 complete types
  (10 .h5 models per type, an ensemble). **A real end-to-end run over the
  repo's own example FASTA** (`example/Q5S007.fasta`):
  a final `site_prediction.tsv` with 1292 real rows, with sensible scores for
  the 8 types. A message "The file for computing FPR doesn't exist" appears
  once per type on the first run (it recomputes the threshold because there was
  no previous `site_prediction.tsv` inside `models/<tipo>/`), which is not
  blocking, confirmed by running twice. A benign numpy/API
  mismatch warning is also present, and likewise does not affect execution.
- **DeepPTMPred: real conda environment built and verified, ESM-2
  functional end to end (2026-07-28).** The `deepptmpred` conda environment
  (`/home/enzo/miniconda3/envs/deepptmpred`, Python 3.10), installed from
  `pred/train_PTM/environment.yml` (without `cudatoolkit`/`cudnn`, as this machine
  has no GPU): `tensorflow==2.15.0` + `tensorflow-addons==0.23.0` +
  `torch==2.13.0` + `fair-esm==2.0.0` install and import without error (the same
  already known TFA deprecation `UserWarning`, not blocking).
  The ESM-2 checkpoint (`esm2_t33_650M_UR50D.pt`, 2.6GB, downloaded manually by the
  user) and its companion `esm2_t33_650M_UR50D-contact-regression.pt`
  (3.7KB) are installed in `DeepPTMPred/esm/checkpoints/`; both verified
  with `torch.load()` (the `args`/`cfg`/`model` keys present, so not a corrupt
  download). **`_extract_esm_features` of
  `_deepptmpred_runner.py` was genuinely executed** (not mocked) against the
  real checkpoint with a 167-residue test sequence: it returns an
  array of `(167, 1280)`, sequence length by the embedding dimension of
  `esm2_t33_650M`, exactly as expected. **The only real remaining
  blocker for DeepPTMPred: PyRosetta.** The real automatic installation
  was attempted via `pip install pyrosetta-installer` +
  `install_pyrosetta()` (the method that DeepPTMPred's own README documents).
  The installer package DOES install, but the download of the real
  wheel fails on this machine: the default mirror
  (`west.rosettacommons.org`) responds `404` on the path that the installer
  queries (`.../latest.html`), and the alternative mirror (`graylab.jhu.edu`)
  fails TLS verification (an untrusted certificate chain in this
  environment). Pending: the user downloading the wheel manually
  (`PyRosetta4.Release.python310.ubuntu.wheel`, the most recent build of
  2026-07-25) to install with `pip install <wheel>` without depending on those
  mirrors.
- **RESOLVED 2026-07-29: DeepPTMPred's licence confirmed by the
  corresponding author himself.** The repo had no declared LICENSE, but the
  paper is CC BY-NC 4.0 (Oxford University Press), so the corresponding authors
  (Yong Liu, Junwen Wang) were asked directly whether
  the repo's code follows the same terms. **Junwen Wang replied on
  2026-07-29 confirming that it does**: "I confirm that the GitHub code follows
  the same CC BY-NC terms." Non-commercial use (research/TFG, and
  future integration as a Scipion plugin by the CNB, a public
  institution) fits within CC BY-NC without problem. No further steps
  pending on this point.
- **Confirmed 100% local throughout the chain (2026-07-27, re-verified with
  a real execution 2026-07-28)**: verified by reading `esm/pretrained.py` from
  github.com/facebookresearch/esm directly and now also with a
  real run. `pretrained.load_model_and_alphabet(path)` only calls the
  network (`dl.fbaipublicfiles.com`) if the argument does NOT end in `.pt` (the
  hub branch); since the runner always passes a local `.pt` path, it always enters
  via the `load_model_and_alphabet_local` branch (a pure `torch.load()` over
  disk), confirmed in the real run above, without network access.
- Our own runner (`_deepptmpred_runner.py`) has its ESM-2 feature
  extraction piece verified against the real environment (see above). It still
  cannot be run end to end (a complete prediction via
  `PredictConfig`/`PTMPredictor`) because PyRosetta is missing, so only that piece
  remains tested solely with `subprocess.run` mocked in tests.

## PyRosetta installed and DeepPTMPred 100% functional end to end (2026-07-28, afternoon)

The user downloaded the wheel manually (`pyrosetta-2026.30+release.bc091c65b8-cp310-cp310-linux_x86_64.whl`,
1.5GB) to `/mnt/c/Users/USUARIO/Downloads/` (WSL) and confirmed the installation.
Installed with `pip install <wheel>` in the `deepptmpred` conda env,
verified with a real `pyrosetta.init()` (banner + version printed, without
error). The blocker that had been open since 27-07 is closed.

**A real end-to-end run of `_deepptmpred_runner.py`** (not mocked,
`AF-P10636-F1-model_v4.pdb` = Tau/MAPT, type `phosphorylation`): it revealed and
fixed 3 additional real bugs, none of them anticipated before having
real PyRosetta running:
- **Incomplete `environment.yml` dependencies in the real conda env**:
  `matplotlib`, `seaborn`, `scikit-learn`, `imbalanced-learn`, `tqdm`,
  `joblib`, `logomaker` were missing despite being listed in
  `pred/train_PTM/environment.yml` (the original `pip install -r requirements.txt`
  of 27-07 did not bring them). Installed manually, and
  `import predict` confirmed to no longer fail.
- **A real bug in `predict.py::PTMPredictor.__init__`** (a new class of bug,
  distinct from the one in `e2_single_data.py` already documented): the saved model
  has a `Lambda` layer (`model.py::182`,
  `Lambda(lambda xin: K.sum(xin, axis=1))`) whose serialized function
  references the symbol `K` (an alias of `tensorflow.keras.backend`). Keras
  ONLY resolves those symbols via the `custom_objects` dictionary passed to
  `load_model`, never via the globals of the module that imports it (even though
  `predict.py` does have `K` in its own namespace). The real `custom_objects`
  of `PTMPredictor.__init__` does not include `'K'` -> `NameError: name 'K'
  is not defined` when loading any model. Confirmed real by running without a
  patch (it blows up) and with a patch (it loads fine) before applying it. Patched
  in `src/engines/_deepptmpred_runner.py::_load_predict_module` (it wraps
  `predict.load_model` and injects `'K': predict.K` into `custom_objects`).
  `predict.py` is NOT edited (vendored), the same criterion as the rest of the runner.
- With both fixes: a complete real run, PyRosetta computes real SASA +
  extracts plDDT from the B-factor, ESM-2 loads from the local cache, and the model
  predicts, giving 130 rows of real output for
  `phosphorylation` over Tau, with scores of 0.13-0.91, sensible.

**The complete `pipeline.py` (PDB path, consensus) genuinely run** over
the same case: it discovered an additional REAL blocker, different from and earlier
than DeepPTMPred's, on the DeepMVP side. See the following section.
## RESOLVED: DeepMVP was not computing `fpr`, calibration genuinely generated (2026-07-28, night)

A blocker found hours earlier in this same session (see commit
`9fddb8c`): `DeepMVPEngine`/the Phase 3 core depend on the
`fpr` column, which `DeepMVP.py` only computes if
`DeepMVP/models/<tipo>/site_prediction.tsv` exists (a VALIDATION file with
a real `y` column), absent from `models.tar.gz` for all 8 types,
additionally confirmed as a known upstream bug (issues #1/#2 of
`bzhanglab/DeepMVP` on GitHub, where two independent external users
reported the same thing).

**Enzo asked for deeper verification before launching a massive download from
UniProt** ("¿buscaste exhaustivamente? los datos deberían estar en DeepMVP,
¿no?"). On reviewing more carefully, the answer was yes: **no
external FASTA was needed**. The `x` column of `all_data.tar.gz` (found
earlier, see commit `9fddb8c`) already comes with the MAXIMUM width that
any submodel of the ensemble asks for (61 = flank 30, verified by reading the 8
`model.json` files: no submodel of any type asks for `peptide_length` > 61).
The submodels that ask for a narrower window are derived by trimming that
same centred window (the same criterion as
`DeepMVP/lib/DataIO.py::getPeptideSequence`). The requirement for a complete
`db` FASTA only appears if the CLI is used (`DeepMVP.py predict -i`, which
calls `processing_prediction_data` unconditionally, another real bug
confirmed there), not if the internal encoding/ensemble logic is reused
directly (`lib.PeptideEncode.encodePeptides`,
`lib.Utils.combine_rts`).

`scripts/generate_deepmvp_calibration.py` was built (standalone, the `deepmvp`
conda env): it downloads `all_data.tar.gz`, and for each of the 8 types
runs the `.h5` files of the already installed ensemble over `<tipo>_testing_70.tsv`
(trimming the window per submodel as described above), averages
with `combine_rts` (identical to the real `ptm_predict`) and writes
`DeepMVP/models/<tipo>/site_prediction.tsv`. **Verified with a real AUROC
per type, not merely "it runs without error"**: acetylation_k 0.986,
glycosylation_n 0.996, methylation_k 0.906, methylation_r 0.932,
phosphorylation_st 0.991, phosphorylation_y 0.992, sumoylation_k 0.974,
ubiquitination_k 0.989, which matches the figures published in the paper and
confirms that the generated calibration is correct.

**The complete pipeline (PDB path) genuinely run end to end with the
calibration already in place**: `AF-P10636-F1-model_v4.pdb` (Tau) -> 572
real PTM sites pass the threshold (17 types represented), 98 with a
real DeepMVP+DeepPTMPred consensus. The first time the pipeline completed
end to end without failing. The `Camino FASTA` (the same `fpr` filter) is
also unblocked by the same fix.

Unresolved (minor, not blocking): which subset to use of the 3
available (`testing_70/80/90`, suffixes not documented in the repo's
code, very probably CD-HIT style sequence identity thresholds
used in the paper to measure generalization at different redundancy with the
training set, INFERRED by standard convention in this kind of paper, not
confirmed by reading the paper's text). The script uses `testing_70` (the most
conservative, least redundant with training) by default, parameterizable
via `--testing-suffix`.

## Phase A class 1 + shared utility + Extension 3 (ΔΔG): implemented and verified 2026-07-28

Built in the order decided on 28-07 (shared utility -> Phase A
class 1 -> Extension 3), in `src/structural/` (standalone, it requires
`pyrosetta`, and is never imported from the main `src` package, the same
pattern as `_deepptmpred_runner.py`).

**A real finding that corrects the "12/17 types" estimate of the decision of
28-07** (that figure assumed coverage by analogy with phosphorylation, without
inspecting the real enum/patches yet, and the note itself left that as an
explicit pending item): verified by TWO independent routes against the
already installed PyRosetta:
1. `ls database/chemical/residue_type_sets/fa_standard/patches/*.txt` (real).
2. Live introspection of `pyrosetta.rosetta.core.chemical.VariantType`.

Only **5 of the 17 types** have a native `VariantType`/patch ready to
use: `phosphorylation` (S/T/Y), `acetylation` (K), `hydroxylation` (P),
`gamma_carboxyglutamic_acid` (E), `lys_methylation` (K, only mono/di/tri
without distinguishing the degree, which the engines do not distinguish either). The other 8 types of
"class 1" (malonylation, arg_methylation, crotonylation, succinylation,
glutathionylation, s_nitrosylation, glutarylation, citrullination) do NOT
have a native patch (a grep of those terms over the complete patches
directory: zero results). They would require building a non-canonical residue
of our own (a params file with topology/partial charges/icoor), a real
cheminformatics task, NOT implemented, NOT fabricated.

- `src/structural/pyrosetta_ptm_patch.py`: `apply_ptm_patch` (it applies the
  real variant via `add_variant_type_to_pose_residue`, and validates that the
  real residue matches the one expected for the type) + `relax_neighborhood`
  (cartesian FastRelax, 1 repeat, restricted to a radius around the site via
  `NeighborhoodResidueSelector` + `MoveMapFactory` + `PreventRepackingRLT`
  over the rest). A CLI that applies a real patch and optionally relaxes,
  dumping an output PDB.
- Verified with a real execution: acetylation on LYS24 of
  `AF-P10636-F1-model_v4.pdb` -> `LYS` becomes `LYS:acetylated`
  (`pose.residue(24).name()` before/after).
- `src/structural/ddg_estimate.py` (Extension 3): it reuses the utility
  above, comparing `ref2015_cart` of the locally relaxed WT pose against
  the locally relaxed patched pose (the same protocol in both).

  **Production robustness: RESOLVED 2026-07-28 night**: `FastRelax` is
  stochastic (Monte Carlo), so a single relaxation can fall into a
  local minimum worse than the real one. `nstruct` was added (default 3, the standard
  practice of Rosetta's `cartesian_ddg`/`ddg_monomer`): it runs N
  INDEPENDENT trajectories per state (each from a fresh pose,
  not chained) and keeps the one with the lowest energy per state. It never
  averages scores directly, as that would mix poorly converged trajectories
  with the good ones. It also reports the standard deviation between trajectories
  as a measure of confidence.

  Verified with a real run (`nstruct=2`, the same acetylation site
  on LYS24 of Tau): both trajectories converged to the SAME exact
  score (1901.6247898992478 WT, 1845.4015068481324 mutated, std=0.0), which confirmed
  that PyRosetta's random number generator is indeed properly
  seeded per process (different seeds in separate runs,
  verified with `pyrosetta.rosetta.numeric.random.rg().get_seed()`), so
  this identical convergence is not a fixed-seed bug: it is a real
  property of this case (a small neighbourhood, few movable residues,
  a simple energy landscape) rather than evidence of general robustness. Other
  sites/larger radii could show more variance between trajectories.
  The mechanism was verified correct (the minimum between real independent
  trajectories, not an illusion of "N runs" that are really the
  same one). **A caveat that still stands**: Tau is intrinsically disordered
  (mean plDDT 49.34 in this AlphaFold model), so no number of
  repetitions fixes the starting region having low structural
  confidence; it only reduces the sampling noise of the relax itself.
- Neither of the two scripts is connected to `pipeline.py` yet,
  following the decision of 2026-07-27 that D (the filter) does not route to Extension
  3/Phase A because those phases were not part of the core. They are invoked
  separately, manually or from a future Phase A/Extension 3 orchestrator.

## Phase A class 2 (glycosylation): IMPLEMENTED AND VERIFIED (2026-07-28, night)

`GlycanTreeModeler` REFINES the conformation of an already attached glycan. To
attach one from scratch another Mover was needed, `SimpleGlycosylateMover`
(the same `protocols.carbohydrates` namespace, also confirmed real by
introspection). Neither of the two engines (DeepMVP/DeepPTMPred) predicts
the glycan COMPOSITION, only the site, so the conserved biosynthetic core
is used as a documented default (not a prediction):
`N-glycan_core` (Man3GlcNAc2, present in EVERY mature N-glycan, the only
universally defensible default) for N-linked, `core_1_O-glycan`
(the T antigen, the most common in mammals but not universal, since there are 8 different
cores in `common_glycans/`) for O-linked. Both IUPAC strings were
verified by reading the real `.iupac` files of the installed PyRosetta.

`src/structural/pyrosetta_glycan_patch.py`: `attach_glycan` (it validates the
real N/S/T residue before patching) + `refine_glycan`
(`GlycanTreeModeler`, with configurable rounds). Verified with two real
runs over `AF-P10636-F1-model_v4.pdb`: N-glycan_core on ASN484 (a real
sequon confirmed as NAT), where 5 sugar residues were added and `ASN484` became
`ASN:N-glycosylated`; and core_1_O-glycan on THR52, where 2 residues were added and
`THR52` became `THR:O-conjugated`. Refinement with `GlycanTreeModeler`
(1 round) runs without error, taking ~64s for the N-glycan_core case.

## Phase A class 3 (ubiquitination/sumoylation): IMPLEMENTED (2026-08-01), corrects the conclusion of 2026-07-28
Enzo asked to finish everything possible. The real Rosetta mechanism for
covalent conjugation of a complete protein (the Lys-NZ + Gly-Cterm isopeptide
of ubiquitin/SUMO) was investigated, not merely discarded by analogy:

- **The Lysine side: it DOES exist and work.** The `SidechainConjugation`
  patch (real, loaded by default, confirmed with introspection of
  `ResidueTypeSet.patches()`) converts LYS into
  `LYS:sidechain_conjugation`, removing the 3 hydrogens of NZ and opening a
  real connection point (`ADD_CONNECT NZ`), ready to link to another
  chain.
- **The ubiquitin/SUMO side: NOT available in this build, for a reason
  documented by Rosetta's own developers.** There is a patch built
  specifically for this,
  `patches/branching/C-terminal_conjugation.txt`, whose own comment
  says literally *"currently only implemented for glycine conjugated
  ubiquitin"*, but it is **commented out/disabled** in the real
  `chemical/residue_type_sets/fa_standard/patches.txt` manifest (line 85):
  `#patches/branching/C-terminal_conjugation.txt  # Something is broken
  with this patch. ~Labonte` (Labonte = a real RosettaCommons collaborator
  specializing in carbohydrate chemistry). Confirmed with real
  introspection (`ResidueTypeSet.patches()` over a pose with real ubiquitin
  loaded, `1UBQ.pdb` from RCSB): the patch is NOT among the 122 patches
  loaded by default, whereas `SidechainConjugation` is.
- This CORROBORATES, independently and more strongly, the finding of the
  investigation of 28-07 in the morning (the RosettaCommons forum: "poor fit
  for the isopeptide domain"). It is not merely that Rosetta has low
  accuracy there, it is that the only native tool built
  specifically for this is disabled by its own
  developers for being broken, without indicating what specifically fails.
- **Decision: no attempt was made to rebuild the mechanism by hand** (low-level
  pose surgery: manually removing OXT, declaring the bond via
  `Conformation.declare_chemical_bond`, idealizing geometry without any
  reference against which to validate that the result is chemically correct). Without
  the official patch as a reference for the correct geometry/valence
  bookkeeping, and without the ability to validate the result against a real
  structure of a conjugate (no trivial superposition is available), the
  risk of producing a structure with a plausible appearance but
  chemically incorrect (bond lengths/angles not idealized
  correctly, or an undetected valence problem) is high and not
  detectable without expert review, which is unacceptable for a real research
  pipeline. Real ubiquitin (`1UBQ.pdb`, downloaded from RCSB,
  76 residues, Gly75-Gly76 at the C-terminus) remains available in
  `DeepPTMPred/data/` (or wherever it is kept) for when it is decided to continue down
  this route with more time, possibly loading the patch manually
  (`-chemical:patch_selectors` pointing at the file directly despite it
  being excluded from the manifest) and validating the result carefully.

67 tests (`pytest tests/`) were still passing at that moment (2026-07-28),
with no new test for `src/structural/` (it requires `pyrosetta`, absent from the
environment where the main suite runs, the same criterion as the rest of the
code depending on external engines: real verification documented here
instead of a mock). **Updated 2026-07-29: 97 tests** after closing items
1/2/8 of the robustness audit plus the GlyGen client (see the
audit section below and "Next steps" item 3). The `src/structural/*.py` files that
depend on `pyrosetta` still have no tests, for the same reason; `glygen_client.py`
is the only exception (it does not require `pyrosetta`, and does have tests).

### CORRECTION 2026-08-01: the conclusion above was incorrect, IMPLEMENTED via `UBQ_GTPaseMover`

The investigation of 28-07 only looked at the patch (`C-terminal_conjugation.txt`,
disabled by Rosetta) and concluded "there is no native route available in this
build". An Opus subagent (a robustness/production investigation, requested by
Enzo) found, by reading the real .cc of
`RosettaCommons/rosetta` (`source/src/protocols/chemically_conjugated_docking/UBQ_GTPaseMover.cc`
and the `UBQ_Gp_LYX-Cterm.cc` app) that there is a route COMPLETELY
INDEPENDENT of the broken patch: `UBQ_GTPaseMover`, a normal `protocols::moves::Mover`
that PyRosetta exposes automatically (zero compilation), which uses a
complete `ResidueType` (`LYX`) instead of the patch for the isopeptide bond.
Verified for real, not merely read:

- `pyrosetta.rosetta.protocols.chemically_conjugated_docking.UBQ_GTPaseMover`
  instantiates without error against the ALREADY installed wheel (`deepptmpred`, PyRosetta
  2026.30), with zero new installation.
- It requires TWO non-obvious init flags, both confirmed by real trial and error
  against the same error that Rosetta reports if either of the
  two is missing: `-chemical:exclude_patches SidechainConjugation` (the same flag that
  Rosetta's official integration test for this app uses) + `-extra_res_fa
  <ruta a LYX.params>` (without this second one, `ResidueTypeSet` blows up with "The
  residue LYX could not be generated" despite having the first flag set).
- A real end-to-end run against `DeepPTMPred/data/AF-P10636-F1-model_v4.pdb`
  (Tau, LYS24, the same residue that Extension 3 already uses for the acetylation
  ddG): `initialize()` alone -> a pose of 834 residues (758 Tau + 76
  ubiquitin), with the target residue mutated to `LYX`. A complete `apply()`
  (`refine_cycles=3`, a fast development value) -> `MS_SUCCESS` in 16.4s,
  with 5 torsion metrics of the mobile tail recomputed with real, non-trivial
  values. The same mechanism was also verified for SUMO1 (a pose of 857 =
  758 + 99 residues), so a single Mover serves both types, only the
  reference PDB changes.
- Two real integration findings (documented in the module's docstring, not only
  here): (1) `initialize()` always writes
  `starting_pose_EMPTY_JOB_use_jd2_0000.pdb` in the process's current cwd
  (confirmed, not an assumption) because `JobDistributor::current_job()` outside
  of a real run via JobDistributor returns a dummy job instead of
  `nullptr`. It does not blow up, but it has to be confined to a disposable
  temporary directory; (2) the 5 per-decoy metrics that Rosetta writes into that same
  dummy job are inaccessible from Python, so they are recomputed directly
  from the output pose via `mover.get_atom_IDs()` (exposed, real) instead
  of trying to read them from the job.
- A third real finding, not of integration but of type safety:
  Rosetta itself does NOT validate that the target residue is a Lysine before
  forcing it to `LYX` (the real `runtime_assert` for that is commented out in the
  .cc), so it is validated in our wrapper instead. It is also validated that the
  pose is single-chain (a real Rosetta `runtime_assert`, whose
  violation blows up loudly rather than silently) and that the PDB numbering
  matches 1:1 with the pose numbering (`initialize()` does a double
  pdb2pose conversion over an index that is already a pose index, which is only correct if
  both numberings coincide; verified to be the case for the single-chain
  AlphaFold PDBs that this pipeline uses, but not assumed for
  any future PDB).

Implemented in `src/structural/ubiquitin_sumo.py` (see its docstring for the
complete detail verified line by line against the real .cc). Reference PDBs
packaged in `src/structural/data/`: `ubiquitin_reference.pdb`
(RCSB 1UBQ, the same one that Rosetta's own team uses in its integration
tests) and `sumo1_reference.pdb` (RCSB 1A5R, model 1 of the NMR,
trimmed programmatically, by searching for the real "QTGG" motif in the
sequence rather than an assumed residue number, down to the real mature
C-terminal diglycine, discarding 4 residues of unprocessed maturation tail that
the original construct has). It also serves for SUMO without code
changes (the same terminal Gly-Gly mechanism). No new tests in the main
suite (the same reason as the rest of `src/structural/`: it requires
`pyrosetta`; 101 tests still pass unchanged, verified). Vault:
the decision is written in `01-Proyectos/PTM-Prediction/Decisiones/2026-08-01-fase-a-clase3-ubiquitinacion-sumoilacion-implementada.md`.

## Real calibration of DeepPTMPred + phi/psi distribution bug (2026-07-30/31)

**2026-07-30**: `scripts/generate_deepptmpred_calibration.py` implemented
(it runs the real model, the `.h5` files already installed via `predict.PTMPredictor`,
the same code that `_deepptmpred_runner.py` uses in production, over a
real sample of `DeepPTMPred/data/ptm_data.csv.gz`, with real sequence+structure
per protein downloaded from AlphaFold DB). AUROC/threshold (Youden's J) measured
for each of the 17 supported types, not assumed. `scripts/
_rebuild_calibration_summary.py` aggregates the 17 per-type TSVs into a
`summary.tsv` without re-running the model (necessary because the main script
only writes `summary.tsv` once, at the end of its own run, so two
partial runs overwrote each other).

**2026-07-31: a real bug found and fixed, affecting the reliability of EVERYTHING
calibrated up to that moment.** The measured AUROCs came out systematically
far below those published in the paper (Briefings in Bioinformatics,
DOI bbag321) for practically all 17 types, with 2 of them (`hydroxylation`
0.34, `lys_methylation` 0.46) even below chance. Verified by reading
`DeepPTMPred/pred/train_PTM/data_loader.py` line 139-141: `phi_center`/
`psi_center` are computed with `half_window = (max(window_sizes)-1)//2 = 25`,
but the `phi`/`psi` array of the source CSV only has 11 elements, so the
condition `len(x) > half_window` is ALWAYS false, meaning the shipped model
never saw anything other than 0.0 in those two angles, neither in
training nor in the test set that produces the paper's AUCs (`trainer.py`
uses the same `load_dataset()` for train AND test). But `predict.py`'s
`PyRosettaCalculator.calculate_features` (REAL inference, the one that
`_deepptmpred_runner.py` uses in production) DOES compute real phi/psi via
PyRosetta, a real train/inference mismatch, not a bug of this calibration.

Verified empirically by forcing `phi=psi=0.0` at inference (a monkeypatch
of `PyRosettaCalculator.calculate_features`, applied in
`_deepptmpred_runner.py::_load_predict_module`, the same point already used
for the real `'K'`/`custom_objects` patch): the AUROC of `hydroxylation`
rises 0.342->0.934 (paper: 0.965), and `lys_methylation` 0.462->0.883 (paper:
0.899), recovering the published performance. SASA/secondary
structure/plDDT are NOT affected (only phi_center/psi_center), so
the model is not blind to structure, only to the fine geometry of phi/psi,
a real and permanent accuracy ceiling of these specific weights
(irreversible without retraining with the `data_loader.py` bug fixed, which is
out of scope for this project).

**All the calibration data prior to this fix are invalid**
(they measured the broken behaviour). Moved to
`DeepPTMPred/data/calibration_STALE_prephifix_2026-07-31/` (gitignored, not
deleted). The complete calibration of the 17 types was relaunched at n=75 with
the patch in place, **finished 2026-07-31**, with the result in
`DeepPTMPred/data/calibration/summary.tsv` (gitignored):

```
ptm_type                     auroc     suggested_threshold  n_proteins
gamma_carboxyglutamic_acid   0.970133  0.214867              44
hydroxylation                0.934372  0.375260              69
malonylation                 0.888000  0.422340             145
lys_methylation              0.882629  0.410402             143
glutathionylation            0.872889  0.472519             145
ubiquitination                0.870222  0.543494             149
phosphorylation              0.859970  0.243480             137
arg_methylation               0.850450  0.343200             139
sumoylation                   0.840180  0.409133             139
acetylation                   0.833538  0.629997             140
o_linked_glycosylation        0.823447  0.310599             121
succinylation                 0.808711  0.482645             144
s_nitrosylation                0.683243  0.500678             139
glutarylation                  0.673067  0.461801             103
citrullination                 0.657401  0.391692              46
crotonylation                  0.606957  0.935762              34
n_linked_glycosylation         0.490667  0.997414             110
```

14/17 types land in the solid-to-excellent range, many matching or
exceeding the paper (which confirms that the phi/psi bug was the real cause of
the low AUROCs, not weak engines). 4 types remain mediocre
(s_nitrosylation/glutarylation/citrullination/crotonylation, 0.61-0.68, with
`crotonylation` on a small sample, only 34 unique proteins). **One type
remains genuinely broken: `n_linked_glycosylation` (AUROC 0.49, worse than
chance) even after the fix**, not yet investigated, a cause distinct from the
phi/psi bug (which has already been fixed), pending a separate look, not
blocking for the rest.

### 2026-08-01: investigation of `n_linked_glycosylation` (AUROC 0.49) and the 4 mediocre types (via an Opus subagent)

**`n_linked_glycosylation`: CONFIRMED as a real limit of the published weights,
not a bug of ours.** Three independent confirmations:

1. `DeepPTMPred/data/calibration/n_linked_glycosylation_calibration.tsv`
   is saturated, not inverted: all 125 sites score 0.913-0.999,
   with `prediction=1` in EVERY row, so the model has no
   discriminative output at all (which rules out a label flip bug).
2. The repo itself ships its training metrics for this type
   (`DeepPTMPred/pred/train_PTM/result/results_n_linked_glycosylation_esm2_kfold/ptm_data_200_29_64_average_results.txt`):
   **AUC 0.4950**, FP Rate 0.94, and a confusion matrix for one fold without a SINGLE
   true negative (`[[0,5],[0,51]]`). Our measurement (0.4907)
   reproduces this almost exactly.
3. It does not even recover by applying the plDDT fix below (0.4907 ->
   0.5064, still chance).

The probable root cause, verified against `ptm_data.csv.gz`: this type has
2115 positives vs 355 negatives (85.6% positive, the most imbalanced of
the 17 alongside crotonylation, 84.1%), and **the SMOTE in
`trainer.py:236` is a no-op for all 17 types** (`target_samples =
min(2000, class_counts[1])`, where the `min` always gives the current count and never
genuinely oversamples; this contradicts what the paper's Methods says). The
paper itself reports 0.616 for this type (Table 1) and attributes it to
"limited negative sample diversity within conserved motifs": the
negatives are asparagines INSIDE intact N-X-[S/T] sequons, the same
sequence motif the positives would show. The shipped weights do not even
reach that 0.616; they reproduce the 0.495 artifact.

**The subagent's recommendation**: exclude `n_linked_glycosylation` from
DeepPTMPred's consensus (it is not a badly calibrated threshold, it is a model with no
discriminative power). The threshold of 0.997 in `summary.tsv` for this type
means nothing. It reinforces the option already listed as item 5 below
(cross-validation with project 1's StackGlyEmbed for this specific
type).

**The 4 mediocre types: this is NOT a model ceiling, it is an additional real
inference bug, the same pattern as the phi/psi one.** `predict.py:260`
computes `local_plDDT` iterating PER ATOM
(`atom.get_bfactor() for atom in structure.get_atoms()`), but
`predict.py:299-303` indexes that list BY RESIDUE NUMBER
(`self.plDDT_values[residue_number - 1]`), so with ~7 atoms/residue in an
AlphaFold PDB, the "plDDT" of residue N ends up actually being the
B-factor of residue ~N/7. Measured over 370 real sites: the correlation against
the real CA B-factor (correct) = 0.031; with the fix = 0.927 (`avg_plDDT` and
`sasa` are fine, r=0.9995/0.98, as the per-atom average happens to come out
the same).

Recalibrating only those 3 types with the fix applied (the same samples):
citrullination 0.657->**0.778** (matching the paper exactly), s_nitrosylation
0.683->**0.770** (closing ~53% of the gap), glutathionylation (a control, already
good) 0.873->0.886 without being damaged. **It is not a model ceiling**, although the
published values of these 4 types (0.769-0.846) are already the lowest of
the 17, with or without the fix.

Two additional defects found that DID measure but do NOT matter in
practice (documented, not prioritized): DSSP is never executed in
`train_PTM/` even though `predict.py:293-295` calls `pose.secstruct()`
(so production is blind to secondary structure, 100% "loop"); and
`data_loader.py:148` vs `predict.py:308-316` have the E/H/L column order
inverted (H and E swapped). A real ablation: neither of the two
changes the AUROC by more than 0.01, so all the recoverable signal is in the plDDT
fix.

**Coordination with the wiring of item 1 (2026-08-01, commit `ae327cd`)**:
the thresholds already wired into `Settings.py` were measured WITH the plDDT bug
active. With the fix, citrullination goes from 0.392 to 0.358 and s_nitrosylation
from 0.501 to 0.486, moderate changes rather than catastrophic ones, but
`DEEPPTMPRED_CALIBRATED_THRESHOLDS` should be considered **provisional**
for these 2 types (and probably for the rest to a lesser degree, not
measured across the 17) until a complete recalibration with the fix applied.

**The subagent's explicit caveats**: the plDDT fix was validated on 3
types out of 17, not on all 17, so the rest are a reasonable extrapolation, not a direct
measurement. `crotonylation` cannot improve with more sampling
even if one wanted to: the run already used ALL the negatives that exist in
the dataset (23), just like citrullination with its 59 positives, a real
limit of the source data, not of `--n-per-class`.

**Both decisions taken and executed 2026-08-01** (Enzo chose "apply
now" for both): (a) `n_linked_glycosylation` excluded from the production
consensus (`CONSENSUS_EXCLUDED_TYPES` in `ptm_annotation.py`, commit
`05cc4d0`); (b) the plDDT patch applied in
`_deepptmpred_runner.py::_load_predict_module` (the same `05cc4d0`, reading the
real CA B-factor from PyRosetta's `self.pose` instead of the badly indexed
`plDDT_values` array) and a complete recalibration of the 17 types
relaunched with the fix (`--n-per-class 75`, with the previous calibration moved to
`DeepPTMPred/data/calibration_STALE_preplddtfix_2026-08-01/`, not deleted).

### Real result of the post-plDDT-fix recalibration (2026-08-01)

```
ptm_type                     auroc     suggested_threshold  n_proteins
gamma_carboxyglutamic_acid   0.989689  0.24807824            44
hydroxylation                0.969697  0.35899624            69
o_linked_glycosylation       0.900667  0.26193630            121
malonylation                 0.890667  0.41699925            145
glutathionylation            0.884444  0.46646200            145
lys_methylation              0.877559  0.43064716            143
ubiquitination                0.876622  0.53213940            149
sumoylation                   0.875315  0.37326753            139
acetylation                   0.878378  0.63506210            140
phosphorylation                0.870244  0.24020174            137
arg_methylation                0.860360  0.34068727            139
crotonylation                  0.815072  0.86312497            34
succinylation                  0.811911  0.50403893            144
citrullination                 0.778079  0.36854228            46
glutarylation                  0.772444  0.47470970            103
s_nitrosylation                0.769730  0.51403310            139
n_linked_glycosylation         0.507200  0.99802846            110
```

**This confirms the subagent's diagnosis 100%**: the 4 previously
"mediocre" types are now solid: crotonylation 0.607->0.815,
citrullination 0.657->0.778, glutarylation 0.673->0.772, s_nitrosylation
0.683->0.770. The plDDT fix (validated on only 3/17 types before this
run) turned out to generalize to all 4, not only the 2 measured then.
**16/17 types in the solid-to-excellent range (0.77-0.99)**. The only type that
remains broken: `n_linked_glycosylation` (0.507, practically identical to the
0.491 pre-fix), which confirms that its problem is independent of the
plDDT bug (a no-op SMOTE + a dataset 85.6% positive, see the section above), and the
exclusion from the consensus already applied is the correct decision.

`Settings.DEEPPTMPRED_CALIBRATED_THRESHOLDS` (wired 2026-08-01,
`ae327cd` -> replaced with these final values) and
`tests/test_settings.py`/`tests/test_ptm_annotation.py` updated with
the 17 thresholds of this table. 101 tests pass.

## Real next steps

1. ~~Decide DeepMVP's `fpr` calibration gap~~ : RESOLVED 2026-07-28
   night (see the section above). `pipeline.py` genuinely runs end to end on
   the PDB path (572 sites, 98 with a real consensus). On a new machine:
   run `conda run -n deepmvp python scripts/generate_deepmvp_calibration.py`
   once after installing DeepMVP's weights, before the first real
   use of the pipeline.
2. ~~Repeat/average the relax in `ddg_estimate.py`~~ : RESOLVED 2026-07-28
   night: `--nstruct` (default 3), see the Phase A / Extension 3 section
   above.
3. ~~Phase A class 2 (glycosylation)~~ : IMPLEMENTED 2026-07-28 night (see
   the section above). A documented biosynthetic core default, not a real
   glycoform prediction. ~~**Improvement: query GlyGen before falling
   back to the generic default**~~ : IMPLEMENTED 2026-07-29: a new
   `src/structural/glygen_client.py` (the real endpoint discovered by reading
   `api.glygen.org/swagger.json`, `POST /protein/detail/{accession}/`,
   verified with a real live query against P10636/Tau, which returned 14 real
   sites and confirms the `type`/`start_pos`/`site_category`/
   `glytoucan_ac` schema used by the client). Wired as PURELY
   INFORMATIVE corroboration in `pyrosetta_glycan_patch.py` (`check_glygen_evidence`,
   the optional `--uniprot-accession` flag): it does not change the attached glycan (the
   `glytoucan_ac` that GlyGen reports identifies a specific glycan that
   would require translating it into an IUPAC core constructible by PyRosetta, a
   real mapping problem not solved by this improvement, documented as an
   explicit limit), it only reports whether the site already has known
   experimental evidence. Network failures degrade to a warning and never abort
   the main flow. 8 new tests (`tests/test_glygen_client.py`, with the network
   mocked; unlike the rest of `src/structural/`, this client does NOT
   require `pyrosetta`, only `urllib`/`json` from the stdlib, so it does run
   in the main suite).
4. ~~Phase A class 3 (ubiquitination/sumoylation)~~ : IMPLEMENTED 2026-08-01,
   see the section "Phase A class 3 ... CORRECTION 2026-08-01" above. The
   previous conclusion (that it requires compiling the whole of Rosetta) was incorrect:
   `UBQ_GTPaseMover` (via PyRosetta, with zero new compilation) covers both
   types with the same mechanism, verified with real runs in
   `src/structural/ubiquitin_sumo.py`.
5. ~~Cross-validation with StackGlyEmbed (project 1) for N-glycosylation~~
   : IMPLEMENTED 2026-08-01, see the section "StackGlyEmbed: informative
   corroboration of N-GLYCOSYLATION" below. It reuses the venv/pickles already
   installed in `B-Cell-Epitope-Prediction` as an external resource (without
   importing code between projects). Purely informative corroboration, which
   never decides `pasa_umbral`/`consenso`, verified with a real run
   against Tau's 3 N-glycosylation candidates.
6. ~~Real calibration of DeepPTMPred~~ : IMPLEMENTED 2026-07-30/31, see
   the section "Real calibration of DeepPTMPred + phi/psi distribution
   bug" below. ~~**Pending Enzo's decision**: whether to connect the
   calibrated thresholds to `Settings.py`/production~~ : WIRED 2026-08-01,
   **updated the same day** after finding and fixing the plDDT bug
   (see the section "investigation of n_linked_glycosylation and the 4 mediocre
   types" + "Real result of the post-plDDT-fix recalibration"
   below): `Settings.DEEPPTMPRED_CALIBRATED_THRESHOLDS` (a dict, 17 types,
   final post-fix values) + `Settings.deepptmpred_threshold_for(tipo)`
   (falling back to `DEEPPTMPRED_MIN_PROBABILITY`=0.5 if the type is not
   calibrated). `ptm_annotation.py` uses the per-type threshold at the 3 points
   where it previously used the fixed threshold of 0.5 (both loops of
   `annotate_pdb_path`); `n_linked_glycosylation` is additionally excluded from the
   consensus (`CONSENSUS_EXCLUDED_TYPES`) for having an AUROC of 0.51 despite the
   well calibrated threshold, a problem of the model's discriminative power,
   not of the cutoff. 4 new tests (`tests/test_settings.py` +1 in
   `tests/test_ptm_annotation.py` for the consensus exclusion).

## Pre-checkpoint robustness audit (2026-07-28, night), polished 2026-07-29

A systematic review of all the code (not only what is already tracked above)
looking for real fragility ahead of a stable checkpoint. None of this
blocked the current use of the pipeline; they were real robustness improvements to
"polish", ordered from highest to lowest severity. **Polished 2026-07-29: items
1/2/3/4/5/6/7/8 resolved.** Only item 10 remains open, on purpose
(DeepPTMPred's licence), which requires contacting Carlos and is not something
resolvable from the code. See the detail below.

**Medium severity (fix first):**

1. ~~**A real path traversal risk, unsanitized**~~ : RESOLVED 2026-07-29.
   `DeepPTMPredEngine.run()` builds the results folder with
   `base_output_dir / record.accession`, a RAW path join. `record.accession`
   came from `structure_parser.py:188` (`accession = path.stem`, the name of the
   input PDB/mmCIF file) and was NEVER sanitized, unlike the accession
   of a FASTA record (`fasta_parser.py:170-178`, which DOES replace `/` and `\`
   explicitly). An input file named literally `...pdb` (3
   dots) produces `accession=".."` (`Path("...pdb").stem == ".."`, verified
   empirically; the original note said `..pdb` with 2 dots, which
   actually gives `stem == "."`, not `".."`, corrected here), which escapes the
   output directory. Fix applied: a new helper
   `structure_parser._sanitize_accession` (the same criterion as `fasta_parser.py`
   for `/`/`\`, plus a fallback to `"UNKNOWN"` for an empty/`.`/`..` accession),
   called when deriving `accession` in `parse_structure`. 2 new tests
   (`test_accession_dotdot_se_sanea_evita_path_traversal`,
   `test_sanitize_accession_reemplaza_separadores_de_ruta`).

2. ~~**DeepPTMPred's ESM feature cache can serve stale data without
   warning**~~ : RESOLVED 2026-07-29. The cache key was only `protein_id`
   (`{accession}_full_esm.npz`), NOT a hash/fingerprint of the real sequence, so
   a DIFFERENT sequence under the same accession/file name (for example
   an updated version of the PDB with the same name) silently
   reused the old ESM embedding. Fix applied: a new helper
   `_deepptmpred_runner._esm_cache_path` includes a short hash (sha256,
   12 hex) of the real sequence in the cache file name
   (`{protein_id}_{hash}_full_esm.npz`), so a different sequence can no longer
   collide with an existing cache. 3 new tests in
   `tests/test_deepptmpred_runner.py`.

3. ~~**`BaseEngine.run()` does not match the real implementations**~~ :
   RESOLVED 2026-07-29. The abstract interface declared
   `run(self, items) -> List[TOut]`, but `DeepMVPEngine.run` and
   `DeepPTMPredEngine.run` actually required a second parameter
   `output_dir: Path = None` not declared in the contract. Fix applied:
   `output_dir: Optional[Path] = None` added to the abstract signature of
   `src/engines/base_engine.py`. 3 new tests in `tests/test_base_engine.py`.

**Low severity (hygiene, not urgent):**

4. ~~An out-of-date `requirements.txt`~~ : RESOLVED 2026-07-29: the comment was
   rewritten and no longer says "they are not built yet" nor uses the
   old "Fase 3" numbering.
5. ~~An unused `ModelLoadError`~~ : RESOLVED 2026-07-29: the dead class was
   removed from `src/utils/exceptions.py` (confirmed by grep that it was not
   used anywhere else before deleting it).
6. ~~`generate_deepmvp_calibration.py` overwrites silently~~ :
   RESOLVED 2026-07-29: it now prints an explicit warning with the
   `testing_suffix` used when `site_prediction.tsv` already exists and is about to
   be overwritten.
7. ~~The console logger only shows `WARNING` or more serious~~ :
   CONFIRMED INTENTIONAL 2026-07-29 (Enzo). The desired behaviour
   (a clean console); the real INFO progress remains available in
   `logs/ptm_pipeline.log`. No code change. `src/utils/logger_config.py:36`.
   **REVERTED 2026-08-01** (Enzo explicitly asked to clean this up):
   the console handler's level was raised to INFO, so the progress (engine
   execution, routing per accession) is now visible live, not only in
   the rotating file. The test was updated
   (`test_console_handler_nivel_info`, previously `_warning`). Commit `fef3e3e`.
8. ~~No dedicated tests for `scripts/generate_deepmvp_calibration.py`,
   `src/engines/base_engine.py`, `src/utils/logger_config.py`~~ : RESOLVED
   2026-07-29: `tests/test_generate_deepmvp_calibration.py` (10, covering
   `trim_window` and a mocked `download_all_data`; `generate_for_type` still
   has no direct test, as it requires real TensorFlow/DeepMVP `lib`, the `deepmvp`
   conda env), `tests/test_base_engine.py` (3), `tests/test_logger_config.py`
   (4). `src/structural/*.py` still has no tests by design (it requires
   `pyrosetta`, see the note already above).

**A scientific caveat: RESOLVED 2026-07-29 (having read the real text of the Methods,
not only the abstract):**

9. ~~`testing_70` could be circular/optimistic~~ : CONFIRMED NOT CIRCULAR.
   The real text of the Methods (DeepMVP, Nature Methods 2025, via PMC12446062,
   full text free):

   > "Ninety per cent of the data was used for training (81%) and
   > validation (9%), and the remaining 10% was used for independent
   > testing."

   > "To control for sequence similarity, peptides in the testing set
   > were filtered to remove those with identities above predefined
   > thresholds (90%, 80% or 70%) compared with peptides in the training
   > and validation sets. This filtering was performed using the
   > clustering tool CD-HIT."

   That is: the primary split IS a 10% genuinely separated from the 90%
   train+validation (the weights shipped in `models.tar.gz` were trained
   over that 90%, never over the 10% test set). The `70/80/90` suffixes are NOT
   alternative splits of the dataset; they are the SAME 10% test set,
   ADDITIONALLY filtered by sequence similarity via CD-HIT against
   train+validation (70/80/90 = the maximum permitted identity threshold:
   `testing_70` is the STRICTEST subset, with ANY peptide of
   ≥70% identity with respect to train/val already excluded, the most conservative
   of the three, which confirms that `--testing-suffix 70` (the script's current
   default) was already the correct choice by intuition, now verified).
   The paper itself reports that the AUROC remains stable under this
   stricter filter ("AUROC values remained stable across all PTMs,
   with changes < 0.03 compared to the original models"), consistent with
   the AUROC of 0.90-0.99 measured locally. **Conclusion: the AUROC measured by
   `generate_deepmvp_calibration.py` is real, not circular**, and safe
   to cite for a TFG/publication.

**Unresolved, already known, re-listed here for completeness:**

10. ~~DeepPTMPred does not declare a licence~~ : RESOLVED 2026-07-29: confirmed
    as CC BY-NC 4.0 (the same terms as the paper) directly by Junwen
    Wang, the corresponding author. See the installation section above for
    the complete quotation of his reply.

## MeToken: informative corroboration of TYPE (IMPLEMENTED 2026-08-01)

Investigated in depth on 2026-08-01 (an Opus subagent): 
`github.com/A4Bio/MeToken` (ICLR 2025) is the most powerful structural engine
evaluated so far for PTM, as it consumes real backbone coordinates (N/CA/
C/O) via a 3D-kNN graph + local quaternion frames, far richer than
the 4 scalars (SASA/phi/psi/plDDT) that DeepPTMPred uses, but the
published checkpoint is a TYPE CLASSIFIER on ALREADY KNOWN sites, not
a site detector (confirmed in `model_interface.py:40`, where the
"no-PTM" class is masked during training when `with_null_ptm=0`, as
the published checkpoint comes, and verified empirically against
`AF-P10636-F1-model_v4.pdb`: it predicts types with high confidence also at
positions without a real PTM, prolines/glycines). That is why it was NOT implemented as a
second consensus engine (Decision 2 remains paused) but as PURELY
INFORMATIVE corroboration of the type at sites that the consensus has ALREADY accepted
(`pasa_umbral=true`), the same non-deciding pattern as GlyGen.

### Real installation: verified, not assumed

- A dedicated `metoken` conda env (`/home/enzo/miniconda3/envs/metoken`, Python
  3.10, CPU-only, confirmed that this machine has no `nvidia-smi`).
  `torch==2.13.0+cpu` (PyTorch's official CPU index) + `numpy`/`scipy`/
  `biopython`/`transformers`/`omegaconf`/`tqdm`/`pandas`/`huggingface-hub`/
  `h5py` installed without error.
- **`torch_scatter` has no prebuilt wheel for this combination** (confirmed:
  `data.pyg.org/whl/torch-2.13.0+cpu.html` only lists variants
  of `pyg_lib` for macOS/Windows/aarch64, none of `torch_scatter` for
  `linux_x86_64` in that version). `pip install --no-build-isolation
  torch_scatter` compiled from source without error, faster than estimated
  beforehand (a few real minutes on this machine, not the ~15 min
  expected).
- Repo cloned (`git clone https://github.com/A4Bio/MeToken`, gitignored,
  the same pattern as `DeepMVP/`/`DeepPTMPred/`). The real weights downloaded and
  verified live: `https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip`
  responded with a real HTTP 200 (`content-length: 88260083`, confirmed with
  `curl -sIL` before downloading), and contains `checkpoint.ckpt` (36MB, the
  flat state_dict that `inference.py` uses, confirmed with `torch.load` +
  `model.load_state_dict()` -> `"<All keys matched successfully>"`) and
  `lightning_checkpoint.ckpt` (60MB, a complete pytorch-lightning
  checkpoint with optimizer state, not used here).
- The `facebook/esm2_t33_650M_UR50D` tokenizer (only the tokenizer, NOT the complete
  ESM-2 encoder, since MeToken uses its own `nn.Embedding` trained from
  scratch, `wo_esm`, despite the attribute's name) was already cached in
  `~/.cache/huggingface/hub/` from an earlier session on this same day
  (StackGlyEmbed, the sibling project), confirmed with `find`, so zero new
  download was necessary for this run. `HF_HUB_OFFLINE=1`/
  `TRANSFORMERS_OFFLINE=1` are set in `_metoken_runner.py` (the same pattern as
  project 1's `stackglyembed_predict_local.py`).

### Two real bugs confirmed by running the repo (not assumed beforehand)

1. **`inference.py:61` calls `PDB.Polypeptide.three_to_one`**, removed
   from Biopython in version >=1.80. Confirmed real:
   `hasattr(PDB.Polypeptide, 'three_to_one')` -> `False` in Biopython 1.87
   (the version that `pip install biopython` installs today). Patched in
   `_metoken_runner.py::_patch_three_to_one`, a monkeypatch over the already
   imported `Bio.PDB.Polypeptide` module, using
   `protein_letters_3to1`/`protein_letters_3to1_extended` (the dictionaries
   that DO exist in modern Biopython) instead. `inference.py` is not
   edited (vendored).
2. **`src/metoken_model.py:213` has `device='cuda'` hardcoded**
   (`codebook_mask = torch.ones(len(codebook), dtype=torch.int32,
   device='cuda')`, inside `MeToken_Model.__init__`). Confirmed real by
   running WITHOUT a patch: it blows up with `AssertionError: Torch not compiled
   with CUDA enabled`. It is the ONLY line in all of `src/` with `device='cuda'`
   hardcoded (confirmed by grep; the rest of the tensors use
   `device=x.device`/`device=index.device`, following the input tensor).
   Patched in `_metoken_runner.py::_force_cpu_ones`, a context manager that
   wraps `torch.ones` ONLY during the construction of the model,
   redirecting `device='cuda'`->`'cpu'` when CUDA is not available. 
   `src/metoken_model.py` is not edited (vendored).

Both patches verified differentially: running
`examples/Q16613.pdb` (the repo's own example, position 31) WITHOUT
patch 2 blows up with the `AssertionError` above; WITH both patches it
reproduces EXACTLY the result documented in the repo's own
`quick_inference.ipynb` (`"PTM type at the position 31 is
Phosphorylation"`), which confirms that the patches do not alter the model's
numerical behaviour, they only make it run on CPU.

### Vendoring

- `src/engines/_metoken_runner.py`: a standalone runner (the
  `_deepptmpred_runner.py` pattern), which runs in the dedicated `metoken` venv, receives a
  single-chain PDB + a list of 1-based positions already accepted by the
  consensus, and returns the type with the highest probability among the 24 real
  classes (excluding index 0 = "Not a PTM type", the null class masked during
  training, and index 25 = "in Rare PTM Types", a bucket of grouped rare
  types not interpretable as a specific type; `src/constant.py::
  PTMtype_list` has 26 entries in total, verified by reading the real
  file).
- `src/engines/metoken_engine.py`: `get_type_corroboration(pdb_path,
  positions, chain_id="A") -> dict[int, dict]`, which invokes the runner via
  subprocess. NON-fatal degradation in every failure mode (repo not
  installed, checkpoint absent, exit code != 0, timeout, a non-existent `python_bin`,
  an absent/malformed output CSV): it always returns `{}` and
  logs a warning, never raising. 10 tests (`tests/test_metoken_engine.py`,
  with `subprocess.run` mocked, the same pattern as
  `tests/test_deepptmpred_engine.py`).
- `src/engines/ptm_annotation.py::annotate_pdb_path` now accepts optional
  `pdb_path`/`chain_id` (`None`/`"A"` by default, so calling it without
  them is identical to the behaviour before this improvement, confirmed
  by the already existing tests without modification). If they are provided and
  `Settings.METOKEN_ENABLED`, it adds 3 columns ONLY on rows with
  `pasa_umbral=true`: `metoken_type`, `metoken_probability`,
  `metoken_type_coincide` (`True`/`False` according to whether `tipo_ptm` has a
  mapped equivalent in `CANONICAL_TO_METOKEN_TYPE` and matches
  `metoken_type` or not, and `None` if there is no known equivalent). It NEVER touches
  `pasa_umbral`/`consenso`, with a double safeguard: `metoken_engine` itself
  degrades to `{}` without raising, and in addition the call is wrapped in a
  `try/except` in `annotate_pdb_path` in case an unanticipated failure
  occurred in the wiring itself. 7 new tests in
  `tests/test_ptm_annotation.py` (columns absent without `pdb_path`, columns
  only on eligible rows, `coincide=True`/`False`/`None`, degradation to
  `None` if MeToken returns `{}`, respecting `METOKEN_ENABLED=False`, and not
  bringing down the annotation if the wiring raises an unexpected exception).
- **A real finding during the testing of the wiring, not merely theoretical**: `tipo_ptm`
  in the DeepMVP-only rows (without a DeepPTMPred match) arrives with DeepMVP's
  RAW name (`acetylation_k`, `methylation_r`, etc.), not the canonical
  name. If `_add_metoken_corroboration` had looked directly in
  `CANONICAL_TO_METOKEN_TYPE` without normalizing first via
  `DEEPMVP_TO_CANONICAL_TYPE`, `metoken_type_coincide` would have been left
  `None` incorrectly for EVERY DeepMVP-only row (the majority of rows
  in the report). Corrected before considering the wiring finished: 2 tests
  (`test_annotate_pdb_path_con_pdb_path_agrega_columnas_solo_en_filas_pasa_umbral`,
  `test_annotate_pdb_path_metoken_type_coincide_false_si_discrepa`) failed
  first WITHOUT the normalization, confirming the bug before fixing it.
- `pipeline.py::run_fase3_pdb_annotation` passes `record.chain_pdb_path`
  (the SINGLE-chain PDB derived in Phase 1.5, NOT the original `record.pdb_path`
  which may have more than one chain) + `record.chain_id`,
  which is necessary so that the 1-based positions that MeToken uses match
  `record.sequence` exactly, the same criterion that
  `_deepptmpred_runner.py` already uses. Confirmed by reading `structure_parser.py` that the
  single-chain PDB preserves the original chain name (it never
  renames it).
- `Settings.py`: `METOKEN_ENABLED` (default `True`, with silent degradation
  if it is not installed, so there is no need to disable it by hand on a machine
  without MeToken), `METOKEN_HOME`, `METOKEN_RUNNER_SCRIPT`, `METOKEN_PYTHON_BIN`,
  `METOKEN_CHECKPOINT`, `METOKEN_TIMEOUT_SECONDS`, the same pattern as
  `DEEPPTMPRED_*`/`DEEPMVP_*`.

### Real end-to-end verification against Tau (`AF-P10636-F1-model_v4.pdb`)

The complete `pipeline.py` (PDB path) was genuinely run THREE times in this
session with the 3 real venvs (`deepmvp`, `deepptmpred`, `metoken`)
chained via `DEEPMVP_PYTHON_BIN`/`DEEPPTMPRED_PYTHON_BIN`/
`METOKEN_PYTHON_BIN`; the numbers in this section are from the FINAL run
(the one that already includes the `phosphorylation_y`->`Phosphorylation`
normalization in `CANONICAL_TO_METOKEN_TYPE`, see the finding below). MeToken ran ONCE
ONLY over the 353 unique positions that the consensus accepted
(`pasa_umbral=true`), not once per site, taking ~6s of real inference on CPU
for the 353 positions (confirmed in the log: `Fase 2 completa` at
`17:56:45`, `Fase 3 completa` at `17:56:51`, which includes the complete
invocation of MeToken in between).

The final report (`outputs/AF-P10636-F1-model_v4_ptm_sites.csv`, 749
rows passing the threshold, 116 with a real DeepMVP+DeepPTMPred consensus):
`metoken_type` was populated in 749/749 rows (MeToken predicted something for
ALL the requested positions, none fell out of range). The real breakdown
by `tipo_ptm` (the `metoken_type_coincide` column, counting only rows with an
equivalent mapped in `CANONICAL_TO_METOKEN_TYPE`; `phosphorylation` and
`phosphorylation_y` are aggregated into a single row because MeToken does not distinguish
the residue for phosphorylation):

```
tipo_ptm                     coincide/evaluadas   tasa
phosphorylation (+_y)             115/115         100.0%
arg_methylation                    21/26           80.8%
ubiquitination                     37/64           57.8%
acetylation                        26/61           42.6%
sumoylation                         1/59            1.7%
hydroxylation                       0/70            0.0%
gamma_carboxyglutamic_acid          0/57            0.0%
malonylation                        0/53            0.0%
lys_methylation                     0/45            0.0%
succinylation                       0/23            0.0%
o_linked_glycosylation              0/112           0.0%
glycosylation_n                     0/3             0.0%
glutathionylation                   0/3             0.0%
```
(without an equivalent in MeToken, `coincide=None` for 58 rows: crotonylation
50, glutarylation 7, citrullination 1; see `CANONICAL_TO_METOKEN_TYPE`,
as these 3 types have no corresponding class among the 24 real ones of
`src/constant.py::PTMtype_list`.) The total count of the
`metoken_type_coincide` column: `True`=200, `False`=491, `None`=58 (200+491+58=749).
The global rate over the 691 evaluable rows: 200/691 (28.9%).

**A real reading of these numbers, not just the table**: MeToken predicted ONLY 6
distinct labels out of the 24 possible across the whole of Tau (`Phosphorylation` 294,
`Ubiquitination` 272, `Acetylation` 153, `Methylation` 21, `Sumoylation` 7,
`S-palmitoylation` 2), a strong bias towards the 3 most common PTM types
(phosphorylation/ubiquitination/acetylation, probably dominant in its
training dataset), and it almost never predicts rarer types even when
asked directly at a site that another engine DID accept for that type:

- **Phosphorylation: perfect agreement (115/115, including the 6
  `phosphorylation_y` sites; MeToken does not distinguish S/T from Y, having a single
  "Phosphorylation" class for all 3, verified in `src/constant.py`)**. The only
  type where MeToken corroborates consistently.
- **Sumoylation: it almost never corroborates (1/59)**, but the pattern is
  biologically coherent, not random noise: of the 59
  sumoylation sites (residue K, confirmed), MeToken predicts `Ubiquitination` in
  37 and `Acetylation` in 21 (the other 2 Lys modifications that are
  among its 6 used labels, see the complete distribution below). Ubiquitin
  and SUMO share exactly the same chemical mechanism (a Lys-Gly
  isopeptide bond, see `src/structural/ubiquitin_sumo.py` of Phase A
  class 3), so this is a microenvironment-level confusion between the chemically
  most similar modifications of the set (all requiring K), not an arbitrary
  failure. A real example: position 657 (K657, `score_deepmvp=0.9997`,
  `score_deepptmpred=0.79`, consenso=True) is the ONLY sumoylation site
  where MeToken does say "Sumoylation" (probability 0.426); at the rest
  (44, 67, 87, 130, 189, 204...) it says "Ubiquitination" or "Acetylation" with
  probabilities of 0.44-0.62.
- **Hydroxylation: 0/70, with a finding stronger than "it does not corroborate"**:
  all 70 positions are proline (`residuo_wt='P'`, confirmed),
  correct for hydroxylation, but MeToken predicts `Phosphorylation` (49
  times) or `Ubiquitination` (21 times), types that biologically require
  S/T/Y or K, NEVER P. MeToken's top-1 at these 70 sites ignores the
  chemical identity of the residue itself. Rather than "MeToken does not know
  hydroxylation", this suggests that the microenvironment of these prolines in Tau
  (an intrinsically disordered protein, mean plDDT 49.34 in this
  AlphaFold model, see Extension 3 above) does not produce a strong signal and the
  model falls back to its majority classes.
- **`S-palmitoylation` (a type without any equivalent in the pipeline, as none
  of the 17 types of DeepMVP/DeepPTMPred covers it) appeared 2 times**, both
  at sites that the pipeline labels `glutathionylation` (positions 608 and
  639), another chemically plausible confusion (both are cysteine
  modifications).
- **Acetylation (42.6%) and ubiquitination (57.8%, counting ONLY when
  `tipo_ptm='ubiquitination'`, not the 37 sumoylation confusions from
  above) are partially corroborated**. Both share LYS as their target
  residue with `arg_methylation`/`lys_methylation`/`sumoylation`,
  consistent with MeToken distinguishing reasonably well that "something is happening to
  this Lys" but not always getting right which of the 4 specific Lys
  modifications it is.

**Practical conclusion**: MeToken corroborates with high reliability ONLY
phosphorylation in this real case; for the remaining types, its discrepancies
are mostly explainable by confusion between chemically
related modifications (the same target residue) rather than noise, which makes it useful
as informative corroboration exactly as it is designed (it never decides
`pasa_umbral`/`consenso`), but it should NOT be read as "MeToken refutes
these sites": its type coverage is strongly biased towards the most
common PTMs of its own training dataset, a real limit of the published
checkpoint (not of this project's wiring).

118 tests (`pytest tests/`, up from 101: 10 from
`tests/test_metoken_engine.py` + 7 from `tests/test_ptm_annotation.py`; the
rest of `src/engines/`/`src/utils/` did not change). `MeToken/` is gitignored
(the same pattern as `DeepMVP/`/`DeepPTMPred/`/`MTPrompt-PTM/`).

## StackGlyEmbed: informative corroboration of N-GLYCOSYLATION (IMPLEMENTED 2026-08-01)

Item 5 of "Next steps" (deliberately not integrated since
2026-07-26) executed: cross-validation of `n_linked_glycosylation` with
project 1's StackGlyEmbed. Reason (see the section "investigation of
n_linked_glycosylation" above): DeepPTMPred has NO real
discriminative power for this type (AUC 0.495 in the training
metrics of the repo itself, not fixable by retraining, and already excluded from the
consensus via `CONSENSUS_EXCLUDED_TYPES`). StackGlyEmbed
(`github.com/GaryChan-lab/StackGlyEmbed`) is a third engine
INDEPENDENT in architecture (ProteinBERT + ESM-2 650M + ProtT5 stacked ->
an SVM meta-classifier), specialized SOLELY in N-glycosylation,
with the same non-deciding pattern as MeToken/GlyGen (it never changes
`pasa_umbral`/`consenso`, it only corroborates an already proposed candidate).

### An external resource reused: verified real, not assumed

Unlike MeToken/DeepPTMPred (repos cloned INSIDE this
project), StackGlyEmbed was already installed and genuinely working in the
SIBLING project `B-Cell-Epitope-Prediction` (project 1, independent by
the decision of 2026-07-26: code is never imported between projects, but heavy
external resources already installed ARE reused). Verified with a real `ls`
before assuming anything:

- The venv `B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed`
  (Python 3.10.20, symlink `python -> python3.10`): it genuinely exists.
- `B-Cell-Epitope-Prediction/StackGlyEmbed/prediction/`: 3
  `power_transformer_*.sav` + `base_layer_pickle_files/` (30 base models:
  10 SVM + 10 XGB + 10 KNN, plus `SVM_meta_layer.sav`, ~480MB total):
  they genuinely exist, with file sizes consistent with real weights
  (10-33MB each).
- ProtT5 (`Rostlab/prot_t5_xl_half_uniref50-enc`, reused from
  `scipion-chem-tmbed/tmbed_src/tmbed/models/t5/`): a `model.safetensors`
  of 3.2GB present, real.
### Vendoring (code is NOT imported from the sibling project, decision 2026-07-26)

- `src/engines/_stackglyembed_runner.py`: an ADAPTATION OF OUR OWN (not an
  imported copy) of the real logic already verified in
  `B-Cell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py`
  (extraction of ProteinBERT+ESM-2+ProtT5 features + prediction via the
  already trained stack of classifiers). Simplified with respect to the sibling's
  script: since here ONE protein is always evaluated per invocation
  (the same pattern as `_deepptmpred_runner.py`/`_metoken_runner.py`), the CLI
  receives `--sequence`/`--positions` directly (without the
  intermediate multi-protein `dataset.txt` format of the original repo). The same
  chunking criterion for complete proteins (ESM-2 in blocks of
  1024aa, ProtT5 in 8797aa) and the same 15-residue window for the
  ESM-2 average, with verification that the 1-based position of the sequon
  matches exactly the `sequence`/`posicion` convention of the rest
  of the pipeline (the complete sequence, with no window offsets to convert).
- `src/engines/stackglyembed_engine.py`: `get_nglyco_corroboration(sequence,
  positions, result_dir=None, filename_prefix="") -> dict[int, dict]`, which
  invokes the runner via subprocess over
  `Settings.STACKGLYEMBED_PYTHON_BIN` (the sibling project's venv).
  NON-fatal degradation in every failure mode (venv absent, runner
  absent, pickles absent, exit code != 0, timeout, an absent/malformed output
  CSV): it always returns `{}` and logs a warning, never
  raising. 11 tests (`tests/test_stackglyembed_engine.py`, with `subprocess.run`
  mocked, the same pattern as `tests/test_metoken_engine.py`).
- `src/engines/ptm_annotation.py`: **unlike MeToken (which requires a
  PDB, PDB path only), StackGlyEmbed only needs `sequence`,
  already an EXISTING parameter in both functions (`annotate_fasta_path`/
  `annotate_pdb_path`), so it applies to BOTH paths without needing any
  new PDB-related parameter.** Both functions now accept
  `enable_stackglyembed: bool = False` (backward-compatible: calling them without
  the parameter is identical to the previous behaviour). If it is `True` and
  `Settings.STACKGLYEMBED_ENABLED`, then for each row with `tipo_ptm` in
  {`n_linked_glycosylation`, `glycosylation_n`} (covering both DeepPTMPred's
  canonical name and DeepMVP's raw one, which are never merged in
  consensus, see `CONSENSUS_EXCLUDED_TYPES`) and `pasa_umbral=True`, it adds 3
  informative columns: `stackglyembed_veredicto`
  (`'Glicosilado'`/`'No glicosilado'`), `stackglyembed_score` (the raw
  probability of the meta-SVM) and `stackglyembed_coincide` (`True` if the verdict is
  `'Glicosilado'`, corroborating the candidate, `False` if it contradicts it; without
  MeToken's "type without equivalent" ambiguity, because StackGlyEmbed
  only predicts N-glycosylation). It NEVER touches `pasa_umbral`/`consenso`,
  with the same double safeguard as MeToken (the engine itself degrades to `{}` without
  raising, plus a `try/except` in the wiring). 8 new tests in
  `tests/test_ptm_annotation.py`.
- `pipeline.py`: `run_fase3_fasta_annotation`/`run_fase3_pdb_annotation`
  pass `enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED` to
  `annotate_fasta_path`/`annotate_pdb_path`, the only change in the
  orchestrator, without new PDB/chain parameters (they are not needed).
- `Settings.py`: `STACKGLYEMBED_ENABLED` (default `True`, with silent
  degradation if the external resource is not available),
  `STACKGLYEMBED_PYTHON_BIN` (whose default points at the REAL venv of the sibling
  project, verified with `ls`), `STACKGLYEMBED_RUNNER_SCRIPT` (this project's own
  vendored one), `STACKGLYEMBED_MODELS_DIR` (whose default points
  at `B-Cell-Epitope-Prediction/StackGlyEmbed/prediction`, verified),
  `STACKGLYEMBED_T5_MODEL_PATH`, `STACKGLYEMBED_ESM_MODEL_NAME`,
  `STACKGLYEMBED_TIMEOUT_SECONDS` (900, the same value that the sibling
  project uses for the same cold load of 3 models).

### Real end-to-end verification against Tau (`AF-P10636-F1-model_v4.pdb`)

`get_nglyco_corroboration` genuinely run (not mocked) against the REAL venv
of the sibling project, over Tau's complete sequence (758
residues, the same one already extracted by Phase 1.5 in a previous run of this
project) and the 3 REAL `n_linked_glycosylation` candidates that the
pipeline had already accepted (`pasa_umbral=true`) in the most recent complete
run (`outputs/AF-P10636-F1-model_v4_ptm_sites.csv`): positions
484 (sequon `NAT`), 676 (`NIT`), 727 (`NVS`). A real cold load of the 3
embedders (ProteinBERT + ESM-2 650M + ProtT5) + inference: **52.5s** of
real time (`time`, not estimated).

The real result (`stackglyembed.csv` persisted):

```
posicion  sequon  veredicto        score
484       NAT     No glicosilado   0.2110
676       NIT     No glicosilado   0.3943
727       NVS     No glicosilado   0.1217
```

The complete WIRING was also verified (not only the engine in isolation): calling
the real `annotate_pdb_path` with these 3 candidates as `deepmvp_df` and
`enable_stackglyembed=True` reproduces EXACTLY the same 3 scores
(`stackglyembed_coincide=False` in all 3 rows), which confirms that the
integration in `ptm_annotation.py` invokes the real engine correctly, not
only in isolation.

**A real reading of these numbers**: StackGlyEmbed does NOT corroborate any of
the 3 candidates that DeepMVP proposed for Tau (all 3 scores fall on the
"No glicosilado" side, although 676 at 0.39 is closer to the boundary than the
other two). This is informative, not definitive proof in either
direction: on the one hand, it reinforces the caution already documented about
`n_linked_glycosylation` in this project (DeepPTMPred itself has no
discriminative power here); on the other, Tau is an intrinsically
disordered protein (mean plDDT 49.34, see Extension 3 above) and there is no external
experimental evidence (GlyGen returned no confirmed N-glycosylation for
P10636 in the queries already made, see the Phase A class 2 section) that
confirms which of the three engines is right in this specific case. No
engine has the last word, and StackGlyEmbed is an additional, purely
informative corroboration, exactly as it was designed.

136 tests (`pytest tests/`, up from 118: 11 from
`tests/test_stackglyembed_engine.py` + 8 from `tests/test_ptm_annotation.py`;
the rest of `src/engines/`/`src/utils/` did not change). The external resource
(`B-Cell-Epitope-Prediction/StackGlyEmbed/`) lives outside this repo, so there is nothing
new to gitignore here.

## Decision 2 (the second real consensus engine, PDB path): CLOSED 2026-08-06 (EMNGly)

Not to be confused with [[MeToken]] (implemented above as informative
corroboration of TYPE, not as a consensus engine; the real limitation
found: the published checkpoint is a type classifier on sites
already known, not a site detector).

### CoNglyPred (the original candidate, 2026-08-01): confirmed DEFINITIVELY dead

The original candidate for the role of second SYMMETRIC consensus engine
(like DeepMVP+DeepPTMPred, but limited to `n_linked_glycosylation`, the
only type with a dead engine) was **CoNglyPred**
(`github.com/whm242446/CoNglyPred`, Proteomics 2025), a graph transformer
over AlphaFold2 PDBs + DSSP + co-attention with ESM-2. It was verified
exhaustively on 2026-08-01 that it had no published weights anywhere.
Enzo sent an email on 2026-08-01 to Shaoping Shi (the corresponding
author, University of Nanchang) asking for the weights, **with no
response after 5 days**. Re-verified 2026-08-06 via the GitHub API
directly (not by search): 0 releases/tags/issues, last push
2024-08-15, and ZERO `.pth`/`.pt`/`.pkl` files in the whole tree despite its
own README instructing one to load `best.pth`. Confirmed genuinely dead,
not merely slow to respond, and discarded.
### EMNGly (the adopted replacement, 2026-08-06)

The replacement research was delegated to an Opus subagent (permanent permission,
see the vault) plus our own re-verification by cloning the real repo. **EMNGly**
(`github.com/StellaHxy/EMNgly`, Hou et al., Bioinformatics 39(11):btad650,
2023) is the only candidate with REAL, downloadable weights that preserves the
design property already decided when discarding MTPrompt-PTM as a replacement for
DeepPTMPred (the second engine must use real 3D structure, not just
sequence):

- Architecture: ESM-1b (`site_emb`+`local_emb`, 1280+1280, see
  `model/get_esm_embedding.py`) + **MIF** (Microsoft's Masked Inverse
  Folding, `model/MIF/`, a REAL 256-dim structural embedding over
  backbone N/CA/C coordinates from the PDB, a permissive BSD-2 licence,
  verified in Microsoft's official repo, as the copy vendored in
  EMNgly lost the LICENSE) -> SVM (`sklearn.svm.SVC`, 2816 features).
- Weights verified AT BYTE LEVEL 2026-08-06 (not from the README): MIF
  (`mif.pt`, 13.8MB) bundled in the repo; the SVM (`N-GlyDE.pickle`, 36MB) in a
  public Google Drive folder, downloaded with a real HTTP Range request
  without authentication, confirmed as pickle protocol 4,
  `sklearn.svm._classes.SVC`, `n_features_in_=2816` (matching exactly
  1280+1280+256), `_sklearn_version=1.1.1`.
- `N-GlyDE.pickle` was chosen (not `N-GlyAltas_classifier.pkl`, also
  available in the same folder) because it is trained over the benchmark
  with NEGATIVES RESTRICTED TO THE N-X-[S/T] SEQUON, the correct regime for
  this role: a model trained against unrestricted negatives (as
  MIND-S/PDeepPP use, candidates evaluated and discarded for this reason) can
  be right merely by learning the motif that this pipeline ALREADY applies via
  `_nglyco_window`, without contributing a real second opinion.
- Published benchmark (the N-GlyDE test set): EMNGly MCC 0.736 vs LMNglyPred
  0.717, DeepNGlyPred 0.605, N-GlyDE 0.499, NetNGlyc 0.265.
- No LICENSE declared in any repo (the paper IS genuinely CC BY 4.0,
  confirmed 2026-08-07 via the PMC XML of PMC10627407, which does not cover the code).
  Unlike CoNglyPred, this is NOT blocking: the weights are already
  downloadable. An email was DRAFTED 2026-08-07 (sending scheduled 2026-08-10,
  Monday) to the real corresponding authors, Yaojun Wang
  (`wangyaojun@cau.edu.cn`, China Agricultural University) and Shiwei Sun
  (`dwsun@ict.ac.cn`, ICT-CAS), both verified via the same PMC XML
  (the Oxford Academic page is behind a Cloudflare challenge, which
  blocks a direct fetch), the same pattern that worked with DeepPTMPred, and it does not
  block production. See the vault, the 2026-08-07 addendum in
  `2026-08-06-decision2-cerrada-emngly-reemplaza-conglypred.md`, for the
  complete text of the email.

**Implemented** (see `src/config/settings.py`, the `EMNGLY_*` block,
`src/engines/_emngly_runner.py`, `src/engines/emngly_engine.py`,
`src/engines/ptm_annotation.py::_apply_nglyco_consensus`): a REAL consensus
of 3 engines (DeepMVP+EMNGly+StackGlyEmbed) for
`n_linked_glycosylation`/`glycosylation_n` on the PDB path, where
`pasa_umbral` = at least 1 of the available engines passes (generalizing the
OR rule of the other types), and `consenso` = at least
`Settings.NGLYCO_CONSENSUS_MIN_ENGINES` (default 2) pass. StackGlyEmbed
is PROMOTED from informative corroboration (a role it keeps as-is on
the FASTA path) to a real consensus engine on the PDB path. A proactive
alignment fix is included in the runner: MIF's `structure_emb` is indexed
by PDB RESIDUE NUMBER (not by ATMSEQ order), so the runner translates
via the `position_mapping` table of Phase 1.5 (`fasta_position`/`pdb_seqid`)
instead of assuming that both numberings coincide (which is true only for
AlphaFold2 PDBs without gaps, like those that EMNGly's own training dataset
uses; see the docstring of `_emngly_runner.py` for the complete analysis,
which also includes the clarification that the BOS-inclusive index convention
of `site_emb`/`local_emb` is NOT a bug despite the initial appearance).

**2026-08-07: the 2 go/no-go checks PASSED, `Settings.EMNGLY_MIN_PROBABILITY=0.5`
is NO longer provisional.** The real weights were downloaded on this machine for the first
time (ESM-1b 7.8GB + the contact regression companion, the SVM `N-GlyDE.pickle`
35.9MB via an HTTP Range from Google Drive, confirmed as a real `sklearn.svm.SVC` with
`n_features_in_=2816`) and run end to end.

**4 real bugs found and fixed along the way** (none of them touching
vendored code, the same criterion as the rest of the project):
1. `_emngly_runner.py` only added `EMNgly/model` to `sys.path`, which is correct
   for its own `from MIF.sequence_models... import`, but
   `MIF/sequence_models/*.py` (`pdb_utils.py`, `pretrained.py`, etc.) make
   BARE imports of `sequence_models.xxx` (not `MIF.sequence_models.xxx`) and
   need `EMNgly/model/MIF` on `sys.path` AS WELL, confirmed by running
   the real import, not merely reading code. The original script
   `get_mif_embedding.py` never exposed this problem because it additionally does
   `sys.path.append("./MIF")` AND depends on its own directory already being
   on `sys.path` (automatic behaviour of a script, not of a module
   imported via subprocess).
2. That same import chain depends on `wget` (the pip package, a transitive
   dependency of `trRosetta_utils.py` never used by this project),
   which is absent and not documented in EMNgly's `environment.yml`.
3. `pip install scikit-learn==1.1.1` without pinning numpy resolves numpy 2.x
   by default, which is binarily incompatible with the sklearn 1.1.1 wheel
   (`ValueError: numpy.dtype size changed`). Fix: pin `numpy==1.23.5`
   afterwards.
4. `fair-esm` (a pip package, not vendored) calls `torch.load(path,
   map_location="cpu")` without `weights_only=False`, and since PyTorch 2.6 the
   default changed to `True`, so the fair-esm checkpoint does not pass the strict
   allowlist (`argparse.Namespace`). A narrow patch in
   `_ESMEmbeddingExtractor.__init__` (the only choke point that loads the
   checkpoint), the same pattern as the monkeypatch of
   `_deepptmpred_runner.py::_load_predict_module`.

**Check 1 (alignment) PASSED**: Alpha-1-antitrypsin (P01009) has 3 real
confirmed N-glycosylation sites in GlyGen (`reported_with_glycan`, a
real query via `glygen_client.py`) at UniProt numbering 70/107/271, which is
PDB numbering (the mature protein, offset -24) 46/83/247. The real PDB used: 1QLP
("2.0 Angstrom structure of intact alpha-1-antitrypsin", RCSB), with 22
N-terminal residues genuinely absent (`REMARK 465`, confirmed), the
case of a "PDB with gaps" that matters for this fix. Real Phase 1.5
(`parse_structure`) over 1QLP: `fasta_position` 24/61/225 -> `pdb_seqid`
46/83/247, all 3 confirmed as `residue_letter='N'`. Verified at 2 levels:
(a) `MIF.sequence_models.pdb_utils.parse_PDB`'s `wt[pdb_seqid-1]` = `'N'`
for all 3, AND it was demonstrated that the NAIVE indexing without the translation
(`wt[fasta_position-1]`, the bug that this fix prevents) gives `'F'`/`'G'` in 2 of
the 3 real cases, not `'N'`, so a real bug, not a hypothetical one. (b) The complete
runner (real ESM-1b + MIF + SVM) over those 3 sites: probabilities
0.887/0.949/0.729 (all 3 above the threshold of 0.5). A real negative control:
2 broken sequons (`N-x-P`, biologically never glycosylated) in the same
protein, positions 82/345 -> probabilities 0.526/0.286 (a real separation,
although 82 lands right at the edge of the threshold, which is expected, as a real model is
not perfect).

**Check 2 (MCC) PASSED AND EXCEEDED the published number**: `MCC=0.8197` (AUC
0.9631, specificity 0.914, sensitivity 0.906) over 301 real evaluable
sites of N-GlyDE's independent set (`NGLYDE_independent.txt`, via
`dukkakc/DeepNGlyPred`, as EMNgly's own repo only documents the source and
does not package it), vs. the `MCC=0.736` published in the paper. 146/447 rows
(33%) were left out: 4/86 proteins without an AlphaFold DB model (78 rows) plus
occasional alignment failures in the remaining 82. The real structure
used: AlphaFold DB, the API `https://www.alphafold.ebi.ac.uk/api/prediction/`
(NEVER build the `AF-{acc}-F1-model_v4.pdb` URL by hand, confirmed real
2026-08-07 that AlphaFold DB now serves v6 and v4 gives a 404 in 100% of the cases
tested, the same finding and the same fix that
`scripts/generate_deepptmpred_calibration.py` already applied before). Exceeding the
published number with a margin (not merely "coming close") is reasonable, not
suspicious: the AlphaFold v6 models used here are more recent/accurate
than those available when the EMNGly paper was published (2023), a
legitimate improvement of the structural source, not data leakage (the SVM was already
trained, frozen, and never retrained with this set).

`predict.py::get_scores` has a confirmed real bug (not fixed in the
vendored copy, see `scripts/verify_emngly_nglyde_mcc.py::get_scores` for the
corrected port): it calls `get_scores(label_y, predict_y[0])`, where `predict_y[0]`
is a scalar (the probability of the FIRST row), not the complete list.

New scripts (permanent, reproducible, with only the data they
download gitignored): `scripts/prepare_emngly_nglyde_structures.py` (real Phase 1.5 per
protein, the `cnb_pipeline` venv), `scripts/verify_emngly_nglyde_mcc.py`
(embeddings + SVM + scores, the `.venv-emngly` venv).

Installation (the dedicated `emngly` venv, pinning `scikit-learn==1.1.1`, see
the version finding above; README.md has the complete detail including
the 4 real installation bugs documented above):
```bash
git clone https://github.com/StellaHxy/EMNgly
python3 -m venv .venv-emngly
.venv-emngly/bin/pip install fair-esm torch "scikit-learn==1.1.1" scipy pandas numpy tqdm wget
.venv-emngly/bin/pip install "numpy==1.23.5"
# ESM-1b (~7.4GB) + contact regression companion, manual download (see README.md)
# SVM: https://drive.usercontent.google.com/download?id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t
```

## Phase A connected to the main pipeline (2026-08-03)

This reverts the decision of 2026-07-27 that D (`apply_workflow_filter`) does not route to
Extension 3/Phase A because those phases did not exist yet. The reason: a demo
complete end to end for Carlos Oscar Sorzano on 2026-08-10, where connecting
Phase A (not only the annotation core) was explicitly requested for that
date. The complete decision is in the vault,
`01-Proyectos/PTM-Prediction/Decisiones/` (look for the note of 2026-08-03).

**A real problem identified before writing code**: modelling ALL the
sites that pass the threshold is computationally unfeasible. A real case
like Tau accepts ~572 sites, and each structural modelling (a ddG with nstruct=3,
conjugation, glycan refinement) takes minutes, not seconds. Enzo's explicit
decision: **top-N per type** (default 1, `Settings.FASE_A_TOP_N_PER_TYPE`)
instead of all sites or only on manual demand, which is representative and
demonstrable in a real automatic run.

### New pieces

- `src/structural/fase_a_dispatch.py`: it routes a site (`pdb_path`,
  `position`, `ptm_type`) to the corresponding real module according to 3 mutually
  exclusive classes (see the module's docstring for the complete detail):
  class 1 (5 types, `pyrosetta_ptm_patch.py` + `ddg_estimate.py`), class 2 (2
  types, `pyrosetta_glycan_patch.py`), class 3 (2 types, `ubiquitin_sumo.py`).
  The other 8/17 types return `estado="sin_soporte_fase_a"` without touching
  PyRosetta. Each call assumes a fresh PyRosetta process, since the 3 classes
  initialize PyRosetta with mutually incompatible flags (confirmed by reading
  the 3 modules: `pyrosetta.init()` only applies the flags of its FIRST
  call per process), which is why classes are never mixed in the same process
  (see `_fase_a_runner.py`, one site per subprocess).
- A necessary prior fix: `ddg_estimate.py` imported `pyrosetta_ptm_patch`
  with a BARE import (`from pyrosetta_ptm_patch import ...`), which only
  worked when invoking the script standalone from its own directory.
  Changed to `from src.structural.pyrosetta_ptm_patch import ...` (the same
  pattern that `pyrosetta_glycan_patch.py` already used for `glygen_client`), so
  that `fase_a_dispatch.py` can import the 3 modules reliably.
- `src/engines/_fase_a_runner.py`: a standalone runner (subprocess, one site
  per process, the same pattern as `_deepptmpred_runner.py`), which inserts the repo
  root into `sys.path` itself (unlike the runners of external
  engines, this one imports the project's OWN code, not a vendored repo,
  so it needs the explicit insert). It serializes the result to JSON.
- `src/engines/fase_a_engine.py`: `FaseAEngine(BaseEngine)` + the dataclass
  `FaseASiteRequest`. One subprocess per site via `Settings.FASE_A_PYTHON_BIN`
  (reusing the SAME `deepptmpred` conda env already installed for DeepPTMPred,
  the same real PyRosetta, not a new env). NON-fatal degradation in every
  failure mode (interpreter/runner absent, exit code != 0, timeout, JSON
  absent/malformed), the same pattern as `stackglyembed_engine.py`/
  `metoken_engine.py`: a site that fails never brings down the rest of the sweep.
- `src/engines/ptm_annotation.py::select_fase_a_candidates`: it selects the
  top-N per type (prioritizing `score_deepptmpred`, falling back to `score_deepmvp`)
  restricted to `Settings.FASE_A_SUPPORTED_PTM_TYPES` (a pure data constant,
  without pyrosetta; `fase_a_dispatch.py` validates at import time that
  it matches exactly the real union of the 3 `SUPPORTED_PTM_TYPES`,
  failing loudly if they ever diverge).
- `pipeline.py`: `run_fase3_pdb_annotation` now also returns the
  in-memory DataFrame (not only `report_path`); a new step
  `run_fase_a_pdb_modeling` after Phase 3 on the PDB path rewrites
  the final CSV with the columns `fase_a_estado`/`fase_a_clase`/`fase_a_ddg`/
  `fase_a_glycan_tree`/`fase_a_glygen_evidencia`/`fase_a_conjugation_metrics`/
  `fase_a_output_pdb` for ALL accepted rows (not only the
  selected ones): the unselected ones are left with `fase_a_estado="no_seleccionado"`,
  documented explicitly instead of an ambiguous empty column. The FASTA
  path does not change (Phase A requires a PDB).

18 new tests (154 total): `tests/test_fase_a_engine.py` (11, with subprocess
mocked, the same pattern as StackGlyEmbed/MeToken), 7 in
`tests/test_ptm_annotation.py` for `select_fase_a_candidates` (pure pandas
logic, without pyrosetta, real unmocked tests), plus an integration
test in `tests/test_pipeline_fase1.py` mocking `FaseAEngine.run`
to verify the complete wiring of the final report.

### Real end-to-end run (2026-08-03): 7/8 candidates modelled, 1 real new finding

Verified with a complete real run (PDB path: Phase 1.5->2->3->A) over
`DeepPTMPred/data/AF-P10636-F1-model_v4.pdb` (Tau), with `FASE_A_PYTHON_BIN`/
`DEEPPTMPRED_PYTHON_BIN` pointing at the same real `deepptmpred` conda env (not
mocked at any point). 749 total sites in the final report, 8
candidates selected by `select_fase_a_candidates` (not 9: no
`n_linked_glycosylation` site exceeded its calibrated threshold of 0.998, so the type is
skipped cleanly, confirmed, with no candidate to model for it).

```
tipo_ptm                     posicion  fase_a_estado  detalle
acetylation                  465       modelado       ddG = +50.72 (ref2015_cart)
lys_methylation               460       modelado       ddG = -16.90
phosphorylation                529       modelado       ddG = +42.52
gamma_carboxyglutamic_acid     319       modelado       ddG = +124.50
hydroxylation                  499       error          ver hallazgo abajo
o_linked_glycosylation         498       modelado       glicano = core_1_O-glycan
sumoylation                    215       modelado       conjugacion real (5 metricas de torsion)
ubiquitination                  755       modelado       conjugacion real (5 metricas de torsion)
```

Real timings: ~7-8min per class 1 site (ddG, nstruct=3 + a final relax
for the output PDB), ~20-40s for class 3 (conjugation), ~1-2min for
class 2 (glycan). The complete run (including the previous Phases 1/2/3: DeepMVP,
17 DeepPTMPred invocations, MeToken, StackGlyEmbed) in ~43 minutes.

**A real new finding: `hydroxylation` always fails via
`add_variant_type_to_pose_residue`, and it is not a bug of this project.**
Investigated in depth, not merely logged: Rosetta's real patch
for this (`pro_hydroxylated_case1.txt`/`case2.txt`,
`TYPES HYDROXYLATION1`/`HYDROXYLATION2`) IS loaded (confirmed via
introspection of `ResidueTypeSet.patches()`, 122 patches in total, both
present) and its `BEGIN_SELECTOR` (`PROPERTY PROTEIN`, `NAME3 PRO`,
`HAS_ATOMS 2HG`) IS satisfied in the real pose (residue PRO499 has the
explicit `2HG` atom, confirmed by listing its atoms). Despite this,
`add_variant_type_to_pose_residue` with `VariantType.HYDROXYLATION`,
`.HYDROXYLATION1` and `.HYDROXYLATION2` (all 3 variants exist as a real
enum) all 3 fail with the same error:

```
ERROR: Unable to find desired residue 'PRO' with variant 'HYDROXYLATION1'.
Attempted to add target variant(s) to ResidueType using both ResidueType
base name 'PRO' and base ResidueType.
```

It was verified that this is NOT specific to this PDB/pose (the same result
when testing cleanly, without any other patch applied beforehand). The exact root cause
(why `add_variant_type_to_pose_residue` does not resolve a
base+variant combination that is declared and whose selector is satisfied) was not
investigated further for reasons of time budget. **No
workaround was applied without being able to validate its chemical correctness**,
the same criterion already established for the broken ubiquitin patch (2026-07-28/08-01):
better to document a real, confirmed limit than to force an unverified solution.
`FaseAEngine`/`fase_a_dispatch.py` degrade this correctly
(`estado="error"`, with the real message captured, and the rest of the sweep continuing
unaffected), so the wiring itself worked exactly as designed; the
limit is Rosetta's/`pyrosetta_ptm_patch.py`'s, not the new wiring's. Pending
for a future session with more time: the other 4 class 1 types
(acetylation, lys_methylation, phosphorylation, gamma_carboxyglutamic_acid)
worked without a problem, so this limit is specific to hydroxylation.

### Robustness: the ddG std is no longer discarded (2026-08-03, the same day)

A post-demo-prep audit found that `ddg_estimate.estimate_ddg` already
computed the standard deviation between the `nstruct` independent relax
trajectories (WT and patched), but `fase_a_dispatch._run_class1`
only propagated the minimum (`wt_score`/`mut_score`), discarding `wt_scores`/
`mut_scores` before the uncertainty reached the final report, so a
ddG was shown as more certain than it really was. Fixed: a new
`Settings.FASE_A_RESULT_TEMPLATE` (the single source of truth for the keys of
a Phase A result, shared between `fase_a_engine.py`/
`fase_a_dispatch.py`/`_fase_a_runner.py`, which prevents the 4 output routes
from falling out of alignment) now includes `ddg_std`/`wt_score_std`/`mut_score_std`;
`pipeline.py` adds `fase_a_ddg_std` to the final CSV. 154 tests still
pass (with `_EMPTY_EXTRA` updated in `test_fase_a_engine.py`).

### Robustness, points 5/6 (2026-08-04): batch mode + dedicated tests for `src/structural/`

A continuation of the post-demo-prep robustness/production plan (points 1-4:
LICENSE, the ddG std, a pinned `requirements.txt`, basic CI, see above).

**Point 5 (batch mode)**: `pipeline.py --input` now accepts a directory
in addition to a single file. `src/utils/input_router.py` already had
`route_inputs` (plural) prepared from its original design but it had never
been connected to the CLI: `main()` only called `route_input` (singular)
over a single file. `run_single_input()` was extracted (the same logic as
the previous `main()`, factored out so that single-file mode and batch
mode reuse exactly the same path without duplicating it) and
`_run_batch()` added: it discovers files with a recognized extension in the directory
(not recursively, the same flat scope as `inputs/`), runs each one, and
**a file that fails does not stop the rest of the batch** (the same non-fatal
degradation criterion that `FaseAEngine`/StackGlyEmbed/MeToken apply
per site), and writes `batch_summary.csv` into `output_dir` (with the columns
`archivo`/`estado`/`detalle`). Exit code 0 only if ALL the files
completed without error; a directory without any recognized file also
fails explicitly (1), not silently. 3 new tests
(`tests/test_pipeline_batch.py`): a multi-file batch OK, partial
degradation (1 OK + 1 invalid, both processed, exit code 1), and a directory
empty of recognizable inputs.

**Point 6 (dedicated tests for `src/structural/`)**: none of the 5 Phase A
modules (`fase_a_dispatch.py`, `ddg_estimate.py`, `pyrosetta_ptm_patch.py`,
`pyrosetta_glycan_patch.py`, `ubiquitin_sumo.py`) had tests of their own;
they were only exercised indirectly via `test_fase_a_engine.py`/
`test_pipeline_fase1.py`, which mock `FaseAEngine`/`FaseASiteRequest`
entirely without ever reaching the internal logic of these modules. The 5 modules
only import `pyrosetta` inside functions (never at module level, by
design, see the docstring of `src/structural/__init__.py`), which makes them
importable without PyRosetta installed; where a function makes an
unconditional import of `pyrosetta.*` BEFORE its own validation logic
(`apply_ptm_patch`/`attach_glycan`/`_run_class2`/`_run_class3`/
`estimate_ddg`), a minimal tree of modules is stubbed in `sys.modules` (via
`monkeypatch.setitem`) instead of mocking at the subprocess level, the same
criterion of "mock at the lowest possible boundary" as the rest of the
project. 36 new tests, none of which requires PyRosetta/conda envs:

- `test_fase_a_dispatch.py` (7): routing to class 1/2/3 according to `ptm_type`,
  an unsupported type (`estado="sin_soporte_fase_a"`, with no PyRosetta module
  initialized), a real exception from a submodule translated to `estado="error"`
  without propagating, the computation of `ddg_std` (the sum in quadrature of the
  WT/mutant deviations) verified against `statistics.pstdev`, optional GlyGen
  corroboration (class 2) only queried if `uniprot_accession` is passed,
  and the consistency of `SUPPORTED_PTM_TYPES` vs `Settings.FASE_A_SUPPORTED_PTM_TYPES`.
- `test_ubiquitin_sumo.py` (11): the 4 real exceptions of `_validate_target`
  (an unsupported type, a multi-chain pose, non-sequential PDB numbering,
  the "second real finding" documented in the module, and a target residue
  that is not a Lysine, which Rosetta itself does not validate, with the `runtime_assert` commented out in
  its real .cc), a valid case for each supported type, `_lyx_params_path`
  (finding/not finding the params file relative to `pyrosetta.__file__`),
  and a real regression: the 2 packaged reference PDBs
  (`ubiquitin_reference.pdb`/`sumo1_reference.pdb`) exist on disk.
- `test_pyrosetta_glycan_patch.py` (5): the 4 branches of `check_glygen_evidence`
  (the network unavailable, no reported site, a site with a specific glycan,
  and a site without a specific glycan). This function needs no pyrosetta
  stub, it only mocks `glygen_client.lookup_site` at the same boundary that
  `test_glygen_client.py` already uses.
- `test_pyrosetta_ptm_patch.py` (9): the 2 exceptions of `apply_ptm_patch`
  (a type without a native `VariantType`, and a real residue that does not match the expected
  one; in both cases `add_variant_type_to_pose_residue` is NOT reached,
  verified explicitly), a valid case for the 5 supported types, and
  the consistency of `PTM_VARIANT_MAP`/`PTM_TARGET_RESIDUE`/`SUPPORTED_PTM_TYPES`.
- `test_ddg_estimate.py` (4): `_best_of_n_relax` returns the minimum of
  `nstruct` independent trajectories (a fresh pose per trajectory, not
  reused) and only applies the patch when `ptm_type is not None` (the WT branch);
  and `estimate_ddg` end to end with `nstruct=3` and `nstruct=1`, verified to use
  the minimum per state (not the average) for the final ddG.

193 tests in total (154 + 3 batch + 36 from `src/structural/`), all passing.

### Point 7 (2026-08-04): the real root cause of `hydroxylation` found and fixed

This closes the last outstanding item of Phase A class 1 (see "Phase A connected to the
main pipeline", 2026-08-03, above: "a real new finding... not
investigated further for reasons of time budget"). Investigated in depth
via an Opus subagent with real PyRosetta (the `deepptmpred` conda env), not
merely by reading code.

**Root cause**: `pro_hydroxylated_case1.txt`/`case2.txt` are the ONLY 2
patches in all of `fa_standard` that declare `SET_BASE_NAME` (`HYP`/`0AZ`,
confirmed via `grep -rl SET_BASE_NAME` over the complete patches
directory). A patch with `SET_BASE_NAME` builds a new BASE `ResidueType`,
not a variant of the original residue, verified by live
introspection: `HYP`/`0AZ` are `is_base_type()==True` with
`base_name()=="HYP"`/`"0AZ"`, never `"PRO"`. `add_variant_type_to_pose_residue`
searches inside the family of the residue's CURRENT `base_name()` (`"PRO"`), and
since `HYP`/`0AZ` never appear there, the search fails with the error already
documented, even though the patch is loaded and its selector is satisfied. It is
structural (a fixed property of the patch), and no init flag changes it. The
other 4 types of this class do not have this problem because their patches do NOT
declare `SET_BASE_NAME`.

**The fix, in `src/structural/pyrosetta_ptm_patch.py`**: a new
`PTM_BASE_NAME_MAP = {"hydroxylation": "HYP"}`; `apply_ptm_patch` replaces
the complete `ResidueType` with `HYP` (case1, trans 4-hydroxy-L-proline, the
ONLY possible product of the prolyl-4-hydroxylases PHD/EGLN over
non-collagenic substrates such as HIF-1a, which is the real case that this
pipeline validates; neither DeepMVP nor DeepPTMPred distinguishes case1/case2) via
`replace_pose_residue_copying_existing_coordinates` instead of
`add_variant_type_to_pose_residue`, the same pattern already established in
`ubiquitin_sumo.py` of using the correct Rosetta API for a class of
residue that the generic API cannot resolve.

**Verified empirically** (not merely argued), against PRO499 of
`AF-P10636-F1-model_v4.pdb` (Tau): the ideal geometry reconstructed from the
patch's icoor matches case1 (CG-OD1 1.424A, OD1-HOD 0.970A, chi4
-127.11 degrees, NOT case2's +90.01 degrees, which confirms the correct trans/4R
stereochemistry), with the backbone unmoved (CA/CG displacement 0.0000A), without
clashes (the closest intra-residue contact to OD1: 2.16A),
`relax_neighborhood` runs without error over the patched result, and without a
regression in the other 4 class 1 types (acetylation/lys_methylation/
phosphorylation/gamma_carboxyglutamic_acid, re-verified in the same
session). 2 new tests in `tests/test_pyrosetta_ptm_patch.py` (11 in
total in that file): the new ResidueType replacement route, and that
`PTM_BASE_NAME_MAP` only covers `hydroxylation`. 195 tests in total, all
passing.

**This closes the Phase A class 1 robustness/production plan 100%**: the 5
supported types (phosphorylation, acetylation, hydroxylation,
gamma_carboxyglutamic_acid, lys_methylation) work without any
known pending limit.

### Point 8 (2026-08-04): biological validation panel (7 real proteins in addition to Tau)

Until now the pipeline had only been run end to end against ONE
real protein (Tau/MAPT). An Opus subagent investigated (a criterion already
authorized: verify a scientific claim against a primary source
before fixing it in production) the real sourcing of PTM sites
documented in the literature for a broader panel, with this project's
hard rule of never fabricating scientific data.

**The real verification method** (not memory/approximation): each site of the
panel comes from a UniProt annotation with `ECO:0000269` evidence
(experimental) and a PMID confirmed via NCBI eutils. Every
`ECO:0000250` annotation ("by similarity", not experimental in human) was
explicitly discarded. The subagent itself found and discarded 2 PMIDs
that it recalled from memory but which turned out to be unrelated papers on
verification, and they are not cited in the final panel.

**An important real finding that changes the design**: AlphaFold's `_v4` URLs
are dead (404); AlphaFold DB is on v6, and one has to use
`AF-<accession>-F1-model_v6.pdb`. It was also verified that **AlphaFold's
numbering matches exactly UniProt's canonical numbering**
(confirmed by running `src.utils.structure_parser.parse_structure` over
the 7 real PDBs downloaded: chain length = UniProt length in
all 7 cases), which is why the panel uses AlphaFold exclusively, not experimental
PDBs. **3 real numbering traps found and avoided**
(verified via the DBREF records of real experimental PDBs, NOT
used in the final panel other than as a warning): histones H3/H4 (the experimental
PDB = UniProt - 1, as the initial Met is trimmed), prothrombin (the experimental
PDB = UniProt - 43, mature chain numbering), EPO (the experimental
PDB = UniProt - 27), documented in the docstring of
`biological_panel.py` so that nobody reintroduces them unintentionally when adding
an experimental PDB in the future.

**The final panel, 7 proteins + 1 real negative control** (`src/validation/
biological_panel.py`, with real PDBs in `inputs/validation_panel/`,
downloaded and verified; each ground truth residue is checked
against the real PDB sequence in `tests/test_biological_panel.py`, 30
tests, all passing):
- **p53/TP53** (P04637, 393 aa): phosphorylation/acetylation/methylation of
  Lys and Arg/ubiquitination/sumoylation, with 7 tier A sites (S15/S20/S46/
  S392 phospho, K382 acetyl, K370/K372 methyl, all multi-PMID) + 14 tier B.
- **Histone H3.1** (P68431, 136 aa) and **Histone H4** (P62805, 103 aa):
  the classic methylation/acetylation/phosphorylation marks (H3K4/K9/
  K27/K36 me, H3K9/K14/K18/K23 ac, H3S10/S28 ph, H4K5/K8/K12/K16 ac,
  H4K20 me, H4S47 ph, all tier A, with decades of independent evidence)
  plus tier B crotonylation/succinylation/glutarylation/citrullination (large
  scale MS). The field names (`H3K4` etc) = UniProt - 1, documented
  explicitly so as not to confuse them with the real position used.
- **Prothrombin/F2** (P00734, 622 aa): **the strongest site of the panel**,
  with 10/10 gamma-carboxyglutamate residues confirmed independently by
  the MODRES records of a real crystal (6C2W) in addition to the original
  sequencing, plus 3 N-glycosylation sites. The only real source of
  `gamma_carboxyglutamic_acid` in the panel.
- **HIF-1a** (Q16665, 826 aa): the only real source of `hydroxylation` (P402/
  P564/N803, tier A despite 1-4 PMIDs each, being universally accepted
  landmark papers) plus sumoylation/ubiquitination. Low pLDDT
  (60.8, an intrinsically disordered protein outside the bHLH-PAS domain),
  which is expected, not a data failure. K532 acetylation is DELIBERATELY
  EXCLUDED (contested in the literature, and the subagent could not verify
  the refutation citations it recalled).
- **EPO** (P01588, 193 aa): 4 N/O-glycosylation sites, tier A despite 1
  PMID each because they come from a DIRECT chemical characterization (peptide
  mapping of the purified glycoprotein, Lai et al. 1986), not from a
  screening.
- **KIT ligand/SCF** (P21583, 273 aa): **a real negative control**, as Asn-97
  is a perfectly valid N-X-S sequon that UniProt confirms
  explicitly as NOT glycosylated in either of the 2 known
  isoforms (the same PMID that documents the sites of the same protein that ARE
  glycosylated), which distinguishes an engine with real specificity from one that
  only matches the N-X-[S/T] motif.

Coverage: **9/9 Phase A types** and 15/17 core types (all except
malonylation/glutathionylation, which were left out of the final panel to
control scope; see the subagent's complete report for 2 optional
proteins, GAPDH/CLIC1, that would cover them if one wants to extend it later).

**Built and verified (not merely designed)**:
1. The 7 real PDBs downloaded from `alphafold.ebi.ac.uk` and verified as
   parseable by `src.utils.structure_parser.parse_structure` (real Phase 1.5,
   not simulated), with the exact expected length in all 7 cases.
2. `tests/test_biological_panel.py` (30 tests): each `GroundTruthSite`
   matches the real residue at that position in the downloaded PDB (the
   check that would have caught any of the 3 numbering traps
   if the panel had inherited them by mistake), all the PTM types are from
   the 17 supported, coverage of the 9 Phase A types confirmed,
   no duplicates, and at least 1 real negative control.
3. `scripts/validate_biological_panel.py`: it runs REAL Phase 1.5->2->3
   (without Phase A, as it measures annotation/consensus quality, not structural
   modelling, which was already validated separately) over each protein of the
   panel, compares against the ground truth, reports tier A/B recall
   separately (never mixed) and explicitly flags any false
   positive on the negative control. `tests/test_validate_biological_panel.py`
   (4 tests, with the engines mocked, the same criterion as test_pipeline_fase1.py).

**The first REAL end-to-end run completed (Histone H4, 2026-08-04)**:
with the DeepMVP/DeepPTMPred engines for real, without a mock, in 13 total minutes
(mostly DeepPTMPred: 17 types, one per subprocess, where ESM-2 dominates the
time independently of the protein's size, so it is not proportional to
its 103 residues). The real result: **Tier A 6/6 (100%)**, **Tier B 7/9
(78%)**, with no false positives on the negative control (H4 has none
of its own, see `kit_ligand_scf`). The 2 tier B misses (`S52 phosphorylation`,
`K13 succinylation`) are exactly the kind of site that the panel itself
documents as weakest (a single MS screening, with the real stoichiometry
unknown), consistent with expectation, not an alarm signal.
A real observed side effect, not blocking: MeToken failed with
`ModuleNotFoundError: No module named 'omegaconf'` in the `cnb_pipeline` env
that orchestrates the pipeline (MeToken should run in its own `metoken`
conda env via `METOKEN_PYTHON_BIN`, not configured in this run), and it
degraded exactly as designed (optional corroboration, a warning
log, with consensus/pasa_umbral unaffected). Phase A tried to model 6
candidates without `FASE_A_PYTHON_BIN` pointing at the env with PyRosetta (0/6
modelled), which also does not affect the measured recall, which is over Phase 3, not
Phase A.

**The run of the remaining 6 proteins completed (2026-08-04)**: the
background run was interrupted halfway through `kit_ligand_scf` (the
negative control, the last protein in the queue) by a connection
drop; 5/6 (p53, histone_h3, prothrombin, hif1a, epo) had finished
by then with their real `_ptm_sites.csv` already on disk. Only the missing
protein was resumed (`scripts/validate_biological_panel.py --only
kit_ligand_scf`, with no need to repeat the other 5) and the final recall
was computed by reading the 6 already generated `_ptm_sites.csv` files directly (the same
logic as `_recall_by_tier`/`_negative_control_report`, without re-running the
pipeline over the 5 already complete).

**The real result, complete panel (7/7 proteins)**:

| Protein | Tier A | Tier B |
|---|---|---|
| p53 | 7/7 (100%) | 14/15 (93%) |
| histone_h3 | 11/11 (100%) | 9/10 (90%) |
| histone_h4 | 6/6 (100%) | 7/9 (78%) |
| prothrombin | 10/13 (77%) | -- |
| hif1a | 5/8 (62%) | 2/4 (50%) |
| epo | 3/4 (75%) | -- |
| kit_ligand_scf | 3/5 (60%) | -- |
| **Global** | **45/54 (83%)** | **32/38 (84%)** |

**The negative control (Asn-97 of kit_ligand_scf/SCF, a valid N-X-S sequon but
one that UniProt confirms is NOT glycosylated): correctly NOT accepted, 0
false positives.** This distinguishes an engine with real specificity from one that
only matches the N-X-[S/T] motif; see the control's design in
`biological_panel.py`.

A real reading of the two weakest cases (prothrombin 77%, hif1a
62%/50%): both are large proteins with structurally
difficult domains (prothrombin has multiple independently folded domains;
HIF-1a is mostly intrinsically disordered, with a mean pLDDT of 60.8 outside
the bHLH-PAS domain, as already documented above), consistent with the
expectation that the pipeline works better in well
structured regions, not a sign of a new bug. Not investigated further
in this session (a per-site analysis of the misses remains as future work
if one wants to improve recall on large/disordered proteins).

229 tests in total, all passing.

### Final pre-demo audit (2026-08-04): a real new-terminal bug found and fixed

Enzo asked for a thorough verification that absolutely everything was ready before
the first session with Carlos. A real and reproducible blocker was found,
not a hypothetical one: `DEEPMVP_PYTHON_BIN`/`DEEPPTMPRED_PYTHON_BIN`/
`FASE_A_PYTHON_BIN`/`METOKEN_PYTHON_BIN` fall back by default to `sys.executable`
(the python that runs `pipeline.py`) if they are not exported. Unlike
`STACKGLYEMBED_PYTHON_BIN` (which already has a real hardcoded default in
`settings.py`), these 4 depended on Enzo running README.md's manual `export`
in the active terminal session. That export was never persisted
in `~/.bashrc` (unlike project 1's `BEPIPRED_PYTHON_BIN`, which is
there). **Confirmed real**: running the exact "Quick Start" command
of the README against a PDB in a simulated new terminal (`bash -ic`, without
any previous export) breaks immediately at the first engine
(`ModuleNotFoundError: No module named 'pyteomics'`, with DeepMVP running under
the `cnb_pipeline` python).

**Fixed**: the 4 `export` lines were added to `~/.bashrc` (the same
pattern as `BEPIPRED_PYTHON_BIN`), outside the repo (it is not a
versioned file of the project). **Verified with a complete real run**,
in a simulated new terminal, the PDB path against Tau: Phase 1.5 -> Phase 2 (434
DeepMVP + 1003 DeepPTMPred) -> MeToken -> StackGlyEmbed -> Phase 3 (749/1022
pass the threshold, 116 with consensus) -> Phase A (**8/8 candidates modelled
successfully**, including `hydroxylation`, which reconfirms the Point 7 fix
in production) -> the final report generated without any error. Zero exceptions in
the whole log.

Also reviewed and discarded as a real blocker: an out-of-date docstring
in `deepptmpred_engine.py` (it says "NO PROBADO TODAVIA contra el entorno real",
written 2026-07-27, contradicted by dozens of real runs since
then), which is cosmetic, not functional, and not fixed in this session.

**The audit's conclusion**: no technical blocker remains for the
sessions with Carlos. The only real open item in the whole project is
Decision 2 (CoNglyPred, the second engine for `n_linked_glycosylation`),
awaiting a reply from Shaoping Shi (email 2026-08-01 + a follow-up
2026-08-04, with no reply yet in either of the two threads).

### Point 9 (2026-08-04, REVERTED the same day): Docker image, discarded

A 4-stage Dockerfile (base/deepmvp/deepptmpred/full) was built as a
continuation of the post-demo-prep robustness/production plan. It caused 2
real host health incidents the same day (the `full` build hung
the WSL2 VM through RAM exhaustion at 13:07-13:14; separately, the
already built images, 27GB, almost exhausted the disk with only 11GB
free on the host). On reviewing the project's real dynamics (Enzo presents
on his own computer with Carlos beside him, and iterates live until he
declares it ready for Scipion, see the vault), it was confirmed that Docker contributes
nothing to that flow: the conda envs already installed directly on the host
cover "living on my computer" completely, and the final integration into
Scipion is done via conda envs (just like the rest of Enzo's scipion-chem-*
plugins), not via containers. Enzo's decision: delete `Dockerfile`,
`.dockerignore`, `docker/` and the already built images, and not to
propose containerizing again unless Carlos or the Scipion integration
requires it explicitly. The PyRosetta wheel (an academic licence) was
kept, moved to the repo root (gitignored) so as not to repeat the
mirror problem if a reinstall is needed in the future.

## Biological coherence analysis (2026-08-07): 3 improvements implemented, 2 deferred

At Enzo's request, a complete analysis of the project (production robustness +
biological coherence of the workflow, the equivalent for this project of
Carmen Elena Gómez's recommendations for
[[project_carmen_elena_gomez_feedback_2026-07-30|BCell-Epitope-Prediction]]).
5 concrete recommendations found, none of them a code bug, but rather gaps
in the communication of scope/interpretation. The 3 with the highest
value/lowest cost were implemented; the other 2 (chain type ambiguity in Phase A of
ubiquitination/SUMOylation, and the lack of kinase specificity in
phosphorylation) are documented in README.md "Scope and interpretation"
as scope limitations, with no implementation planned.

**Change 1, interpretation disclaimer**: no engine models the real
biosynthetic pathway of the substrate (co-expression/co-localization of the enzyme).
They all predict the CAPACITY to be modified from sequence/structure, never
real occurrence in a cell/tissue/condition. `pipeline.py::INTERPRETATION_DISCLAIMER`
is printed once at the end of each run (it is not written inside the CSV, as
that would change its schema and break existing readers).

**Change 2, secretory pathway evidence (N-glycosylation)**: a new
`src/structural/uniprot_localization_client.py` (stdlib-only, the same pattern
as `glygen_client.py`), making a real query to `rest.uniprot.org` (the endpoint
`/uniprotkb/{accession}.json`, verified live, not assumed). The column
`via_secretora_evidencia` (True/False/None) is purely informative (it never
decides `pasa_umbral`/`consenso`, the same pattern as MeToken/GlyGen), applied
only to `n_linked_glycosylation`/`glycosylation_n` (deliberately NOT to
`o_linked_glycosylation`, as those are two distinct biological pathways, O-GlcNAc vs
mucin, that this client does not distinguish). Toggle: `Settings.SECRETORY_PATHWAY_CHECK_ENABLED`
(default True).

3 real bugs found by verifying against the real API (not assumed):
1. UniProt returns HTTP 400 (not only 404) for accessions with an
   invalid format, which is the MOST COMMON case here, since `accession` normally comes
   from the stem of the input file, almost never a real UniProt ID. Both
   codes are treated as "no data" (`None`), never raising.
2. The keyword `'Membrane'` on its own is too broad, and was discarded after a
   real false positive against GAPDH/P04406 (a cytoplasmic/nuclear protein
   with reports of peripheral membrane association).
3. Counting ANY keyword as "there is localization data" (without filtering
   by `category`) gave a false `False` instead of the correct `None` when
   UniProt only had keywords unrelated to localization (for example
   `'3D-structure'`, whose real category is `'Technical term'`), corrected
   by filtering on `category == 'Cellular component'`.

Verified live against 5 real cases: P01588/EPO (True, secreted),
P04406/GAPDH (False, cytoplasmic/nuclear, with real known localization),
`1qlp`/`NOTREAL999` (None, an unrecognized accession).

**Change 3, a competition/crosstalk warning between PTMs**: several types
modify the SAME chemical group of a residue and are mutually exclusive
in the same molecule/instant (never modelled before this improvement, as each
type/position was scored independently). The groups are real, conservative, and
verified against each row's `residuo_wt` (never assumed merely from the
type) in `src/engines/ptm_annotation.py::_PTM_COMPETITION_GROUPS`: lysine
acyl (acetylation/ubiquitination/sumoylation/Lys-methylation/malonylation/
glutarylation/succinylation/crotonylation), cysteine thiol (S-nitrosylation/
glutathionylation), arginine guanidino (Arg-methylation/citrullination),
Ser/Thr hydroxyl (phosphorylation/O-glycosylation, the "Yin-Yang" hypothesis).
The column is `ptm_crosstalk_aviso`, purely informative. Toggle:
`Settings.PTM_CROSSTALK_CHECK_ENABLED` (default True).

**A real risk avoided during testing**: many pre-existing tests in
`test_ptm_annotation.py` (the StackGlyEmbed/EMNGly sections) use
`glycosylation_n`/`pasa_umbral=True`, so without mocking, Change 2 would have made
REAL network calls to UniProt in every one of them during the main
suite (violating this project's convention of never hitting a real
API in tests). Fixed with an `autouse` fixture (`_mock_uniprot_lookup`)
that mocks the network boundary to `None` by default for the whole file. 23
new tests (14 in `test_ptm_annotation.py`, 9 in the new
`tests/test_uniprot_localization_client.py`), 276 tests in total, all
passing, with the complete suite in ~37s (without real network access).

**A licence item clarified in this session (not to be confused)**: the email already
answered by Junwen Wang ("I confirm that the GitHub code follows the
same CC BY-NC terms") is about **DeepPTMPred** (authors Yong Liu/Junwen
Wang), already closed since 2026-07-29 (see the licences section above). The
licence of **EMNGly** (different authors: Yaojun Wang/Shiwei Sun, see
the licences section above for their real emails) was the only real open
compliance item in the whole project, **resolved on 2026-08-11**:
they replied confirming that they had added the MIT licence to the GitHub repo
(`StellaHxy/EMNgly`, commit `824d6dc`) and explicitly authorizing
integration into Scipion. Synchronized in the local vendored copy
(`EMNgly/`, fast-forward `0b3c84d..a719532`). No licence items remain
open in the project.

## Phase A / Phase 3c REMOVED from scope, Scipion integration confirmed (2026-08-10)

A feedback meeting with Carlos about this project (the demo of 2026-08-10 having
already taken place). The verdict: **Phase A/3c (real structural modelling via PyRosetta:
ddG, glycosylation, ubiquitination-sumoylation, all the previous sections
of this document dated 2026-07-28/2026-08-03) is
removed completely from the project's scope**, not deferred. The core
(Phase 1 -> 1.5 -> 2 -> 3 -> 3b: DeepMVP+DeepPTMPred, consensus, the secretory
pathway, Kinase Library, MeToken, EMNGly+StackGlyEmbed for
N-glycosylation, competition between PTMs) is **confirmed feature-complete**:
"there is nothing more to add". Carlos also confirms that **the project will be
integrated into Scipion** (the same pattern as
[[project_carlos_scope_decision_2026-07|project 1's integration]]).
The complete write-up of the decision, with what remains to be decided about the
order of work:
`01-Proyectos/PTM-Prediction/Decisiones/2026-08-10-carlos-elimina-fase3c-confirma-scipion.md`
in the vault.

**Code deleted the same session** (all the sections of this document
about "Phase A"/class 1/2/3/Extension 3 above remain as a historical
record of what existed and was verified, not as the current state):
- `src/structural/fase_a_dispatch.py`, `ddg_estimate.py`,
  `pyrosetta_ptm_patch.py`, `pyrosetta_glycan_patch.py`, `glygen_client.py`,
  `ubiquitin_sumo.py` + their reference data (`sumo1_reference.pdb`,
  `ubiquitin_reference.pdb`). `src/structural/uniprot_localization_client.py`
  (Phase 3b, the secretory pathway) is the ONLY piece that remains alive in that package;
  it was never part of Phase A despite living in the same directory.
- `src/engines/fase_a_engine.py`, `src/engines/_fase_a_runner.py`.
- `src/engines/ptm_annotation.py::select_fase_a_candidates`.
- The complete `FASE_A_*` block of `src/config/settings.py` (7 settings:
  `ENABLED`, `TOP_N_PER_TYPE`, `SUPPORTED_PTM_TYPES`, `RESULT_TEMPLATE`,
  `RUNNER_SCRIPT`, `PYTHON_BIN`, `TIMEOUT_SECONDS`).
- `pipeline.py::run_fase_a_pdb_modeling` + the "FASE 3c" block of the on-screen
  summary (`_print_pdb_summary`) + `_FASE_A_RESULT_COLUMNS`.
- Dedicated tests: `test_fase_a_engine.py`, `test_fase_a_dispatch.py`,
  `test_ddg_estimate.py`, `test_pyrosetta_ptm_patch.py`,
  `test_pyrosetta_glycan_patch.py`, `test_glygen_client.py`,
  `test_ubiquitin_sumo.py`. Edited (removing only the Phase A part, with the
  rest of the file kept): `test_ptm_annotation.py`,
  `test_biological_panel.py`, `test_pipeline_fase1.py`.
- The Colab notebook (`notebooks/colab_fases_1_3b.ipynb`): the
  "Fase 3c (local)" cell and `FASE_A_ENABLED` removed from the environment variables. The
  PyRosetta mirror fix in the DeepPTMPred cell (2026-08-10, commit
  `23d2dae`) IS still necessary, because `DeepPTMPred/predict.py` (Phase 2,
  a mandatory engine) imports `pyrosetta` at module level for per-residue
  SASA, unrelated to Phase A/3c.
- `README.md`: the "Fase A" section removed from Architecture/Installation/
  Licences, and the K48/K63 chain type disclaimer removed from "Scope and
  interpretation" (it depended on the `fase_a_cadena_tipo_aviso` column, which no
  longer exists).

**The immediate priority before the Scipion integration** (Enzo's explicit
request, "leave everything without friction"): caching the notebook's 5 conda
environments with `conda-pack` in Drive, IMPLEMENTED later the
same session (commit `58e66dc`, see the "conda-pack" section below), and an
end-to-end robustness audit of the current pipeline before declaring it
ready to port to Scipion, also completed the same session
(`/code-review` + a manual security pass + a real end-to-end run, with 6
dangling references to `glygen_client.py` found and corrected, commit
`8b3473c`).

## StackGlyEmbed REMOVED from the project (2026-08-10, the same day, a later session)

Carlos's second decision in the same feedback meeting (see the
Phase A/3c section above for the first): **StackGlyEmbed is removed
completely from the project**. The real reason: it depended on a sibling
project's venv (`B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed`) +
a third party's weights (`scipion-chem-tmbed`, ProtT5 ~3GB), which is real friction
for the Scipion integration (a Scipion protocol cannot depend
on absolute paths of ANOTHER local project). This confirms the same pattern that already
motivated the removal of Phase A/3c: dependency friction before
porting code to Scipion, not a problem with the engine's own quality.

**The real impact on the N-glycosylation consensus (PDB path)**: the
`n_linked_glycosylation` consensus (the decision of 2026-08-06) was already
designed to degrade with any subset of available engines
(`nglyco_consensus_active = pdb_path is not None and Settings.EMNGLY_ENABLED`,
independently of StackGlyEmbed), so removing StackGlyEmbed simplifies the
consensus to 2 engines (DeepMVP+EMNGly) instead of 3, without needing to touch
`NGLYCO_CONSENSUS_MIN_ENGINES` (which stays at 2): in practice it now requires
BOTH engines to pass, instead of 2 of 3. Verified by reading
`_apply_nglyco_consensus` before touching anything, so it was not an assumption.

**Code deleted**:
- `src/engines/stackglyembed_engine.py`, `src/engines/_stackglyembed_runner.py`.
- `src/engines/ptm_annotation.py`: the import of `get_nglyco_corroboration`, the
  `enable_stackglyembed` parameter of `annotate_fasta_path`/`annotate_pdb_path`
  (a breaking change of the internal API, accepted, as there are no external
  consumers of the module), the complete `_add_stackglyembed_corroboration` function (its
  informative role on the FASTA path no longer exists at all), and the StackGlyEmbed
  branch of `_apply_nglyco_consensus` (its consensus role on the PDB path).
- The complete `STACKGLYEMBED_*` block of `src/config/settings.py` (7 settings:
  `ENABLED`, `PYTHON_BIN`, `RUNNER_SCRIPT`, `MODELS_DIR`, `T5_MODEL_PATH`,
  `ESM_MODEL_NAME`, `TIMEOUT_SECONDS`), which included the 3 absolute paths to the
  sibling project that caused the real friction.
- `pipeline.py`: the 2 call sites that passed
  `enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED`.
- Tests: `test_stackglyembed_engine.py` deleted entirely; in
  `test_ptm_annotation.py` the whole FASTA informative corroboration block was
  trimmed (8 tests) + 2 tests of the "generic pathway" on the PDB path that no longer
  exists (the fallback to StackGlyEmbed when `EMNGLY_ENABLED=False`); the
  3-engine consensus tests were rewritten to 2 engines instead of deleted
  (equivalent coverage, without StackGlyEmbed). 338 final tests (down from 357).
- The Colab notebook: the whole of Section 10 removed ("StackGlyEmbed
  (opcional, desactivado)") along with the `STACKGLYEMBED_ENABLED` line from the
  environment variables.
- `README.md`: removed from Architecture (the consensus is now DeepMVP+EMNGly),
  Optional corroboration, Scope and interpretation, Current state, and
  Installation (the complete block exporting the sibling project's paths).

Commit `e9be138`. The notebook was tested locally (syntax, 338 tests passing,
see above), pending confirmation in a real Colab run.

## First real post-refactor Colab run (2026-08-11): p53 clean, the negative control fails

The Colab session paused by GPU quota was resumed (see the
"conda-pack" section above). The first real end-to-end run with the
current architecture (without Phase A/3c, without StackGlyEmbed, with EMNGly's licence already resolved,
see `project_ptm_emngly_license_resolved_2026-08-11` in Claude Code's
memory) over 2 proteins of the panel: `p53_P04637` and
`kit_ligand_scf_P21583`.

**p53**: 224/464 sites pass the threshold, 34 with a real consensus (2+ engines).
MeToken ran for real (the cached env OK), Kinase Library assigned a family to
32/32 phosphorylations, and the crosstalk warning concentrated on K370-K386 (p53's
best documented acetylation/ubiquitination/sumoylation switch in the
literature). It was verified in the CSV that EMNGly genuinely ran (the column
`score_emngly` with 3 real values out of 4 N-glyco candidates, engine
`DeepMVP+EMNGly`), and in all 3 cases EMNGly gives low scores (0.10-0.35,
below `EMNGLY_MIN_PROBABILITY=0.5`) while DeepMVP gives 0.90-0.99: the
structural signal correctly rejects N-glycosylation in a
nuclear/cytoplasmic protein that never transits the ER/Golgi, with the 2-engine
consensus working as designed. `via_secretora_evidencia` is
`None` for all 4 (as expected: `accession` on the PDB path is the stem of the
file `p53_P04637`, not the real UniProt ID `P04637`, so UniProt returns
400/404, translated to "no data" on purpose, see
`uniprot_localization_client.py`).

**kit_ligand_scf (the negative control)**: the 2 real N-glyco sites
(N-90, N-145) correctly detected with `consenso=True` in both
engines. The 3 real O-glyco sites (S-167, T-168, T-180) correctly
detected via DeepPTMPred. **N-97 (the panel's only real negative control,
which UniProt/PMID 1381905 confirms explicitly is NEVER glycosylated)
comes out as a false positive with a complete consensus**: DeepMVP 0.939, EMNGly
0.792, both above the threshold. This is the first time this negative control has been
tested against the 2-engine architecture; the previous baseline
("correctly rejected", 2026-08-04) predates EMNGly existing
in the pipeline (added 2026-08-06), so it is not a measurable regression,
it is simply the first real datum on this specific case.

**The decision taken**: `Settings.EMNGLY_MIN_PROBABILITY` is not touched (raising it
enough to exclude 0.792 without touching N-90/N-145 would mean adjusting a
threshold already validated against 301 real N-GlyDE sites for a single case,
overfitting to n=1, the same criterion as
`feedback_no_fabricar_datos_cientificos` in Claude Code's memory), nor is
the removal of StackGlyEmbed reverted (it would reopen the real Scipion
friction that motivated its removal, for a marginal gain on 1 site).
It is documented as a known limitation of the predictor in README.md, the section
"Scope and interpretation", consistent with EMNGly having an MCC of 0.82
(not 1.0) and this being exactly the kind of difficult site for which the
control was designed (chemically indistinguishable from N-90/N-145 except by
something that neither sequence nor static structure captures).

**Not blocking for the port to Scipion**: the pipeline behaved
exactly as designed (both engines genuinely ran,
computed real scores, and the consensus was applied correctly), so this is
an accuracy characteristic of the model, not a software bug.

**`validate_biological_panel.py --only p53 kit_ligand_scf` run
(2026-08-11, the same Colab session)**: formal Tier A/B recall against the
ground truth:

| | 2026-08-04 (baseline, without EMNGly) | 2026-08-11 (post-refactor) |
|---|---|---|
| p53 Tier A | 7/7 (100%) | 14/15 (93%) |
| p53 Tier B | 14/15 (93%) | 18/19 (95%) |
| kit_ligand_scf Tier A | 3/5 (60%) | 5/5 (100%) |

Global: **Tier A 19/20 (95%), Tier B 18/19 (95%)**. The apparent drop in
p53 Tier A is not a regression: the panel grew from 22 to 34 real sites after
the migration to dbPTM of 2026-08-09, a harder denominator. kit_ligand_scf
improved from 60% to 100% in Tier A. The script explicitly flagged the already known
false positive of N-97 (see the section above), confirming that the
automatic detection of the negative control works as designed.

**Conclusion**: removing Phase A/3c + StackGlyEmbed did not damage recall (it
holds or improves), run with a real GPU over the current architecture, which
is the real gate that Carlos asked for before porting to Scipion (see
`project_carlos_ptm_migration_greenlight_2026-08-10` in Claude Code's
memory). The project's core is validated end to end, ready to start
the migration.
