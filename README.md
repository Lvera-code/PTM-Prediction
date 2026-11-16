# PTM-Prediction

[![tests](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml/badge.svg)](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml)

Pipeline for predicting post-translational modification (PTM) zones from
FASTA or PDB/mmCIF. Second CNB project (commissioned by Carlos
Oscar Sorzano), independent of
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction).

## Architecture

Phase numbering aligned with project 1
([BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
Phase 1 -> Phase 1.5 (PDB only) -> Phase 2 (engines) -> Phase 3 (core) -> Phase 3b
(informative cross-references, PDB only).

- **FASTA path**: Phase 1 (sanitation) -> Phase 2: **DeepMVP** (single engine,
  6 PTM types: phosphorylation, acetylation, methylation, sumoylation,
  ubiquitination, N-glycosylation) -> Phase 3.
- **PDB path**: Phase 1.5 (ATMSEQ sequence extraction + position
  mapping via `gemmi`) -> Phase 2: **DeepMVP + DeepPTMPred** consensus (17
  PTM types; DeepPTMPred requires a mandatory `pdb_path`, with no
  sequence-only mode, hence the asymmetry between paths) -> Phase 3 -> **Phase 3b**.
- **Phase 3 (core)**: annotation/filtering + flow decision logic over
  the raw Phase 2 predictions (`src/engines/ptm_annotation.py`).
  It merges consensus where DeepMVP and DeepPTMPred agree on type+position.
  `n_linked_glycosylation` is deliberately EXCLUDED from THAT specific
  fusion (a DeepPTMPred model without real discriminative power for this
  type, AUROC 0.51, see STATUS.md), although both engines still
  report it separately. Instead, the PDB path has a different REAL consensus
  for this specific type: **DeepMVP + EMNGly** (decision
  2026-08-06, replacing the original CoNglyPred, which has no published weights, see
  STATUS.md "Decision 2"), `pasa_umbral` = at least 1 of the available
  engines passes, `consenso` = both pass.
- **Optional corroboration** (these never decide `pasa_umbral`/`consenso`, they only
  inform about sites that the consensus has ALREADY accepted): **MeToken** (type, PDB
  path), **secretory pathway evidence via UniProt**
  (the `via_secretora_evidencia` column, N-glycosylation, both paths) and
  **a competition/crosstalk warning between PTMs** (the
  `ptm_crosstalk_aviso` column, when 2+ types competing for the same
  chemical group of a residue pass at the same position), see "Scope and
  interpretation" below.

## Scope and interpretation

Biological coherence analysis 2026-08-07 (the equivalent, for this
project, of Carmen Elena Gómez's recommendations for
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
scope communication gaps found, none of them a code bug.

- **Predicted capability, not observed occurrence.** No engine in this
  pipeline (DeepMVP/DeepPTMPred/EMNGly) models the real biosynthetic
  pathway of the substrate. They all predict from sequence/
  structure whether a site IS modifiable in principle, never whether that PTM
  actually occurs in a specific cell/tissue/condition (that depends
  on the real enzyme being co-expressed and co-localized with the substrate,
  something no sequence predictor can capture). A result from
  this pipeline must be read as "potentially modifiable site", not
  as "this protein IS modified". Printed as a warning at the end of each
  run (`pipeline.py::INTERPRETATION_DISCLAIMER`).
- **No secretory pathway filter for N-glycosylation.** N-glycosylation
  occurs chemically in the ER/Golgi lumen, so a purely
  cytoplasmic/nuclear protein is almost never glycosylated at a sequon,
  regardless of the consensus score. The `via_secretora_evidencia` column
  (via `src/structural/uniprot_localization_client.py`, REAL evidence
  from UniProt, never invented) warns when the consensus has accepted an
  N-glyco site without known evidence of secretory localization. It NEVER
  rejects the candidate (the same criterion as Carmen
  Elena's recommendation not to discard N-glyco candidates in project 1). Its scope is
  limited to `n_linked_glycosylation`; `o_linked_glycosylation` has two
  distinct biological pathways (cytoplasmic/nuclear O-GlcNAc vs mucin in
  the secretory pathway) that this client does not distinguish.
- **No modelling of competition between PTMs of the same residue.** Several
  types modify the same chemical group of a residue and are mutually
  exclusive in the same molecule/instant (lysine acyl: acetylation/
  ubiquitination/sumoylation/methylation/malonylation/glutarylation/
  succinylation/crotonylation; cysteine thiol: S-nitrosylation/
  glutathionylation; arginine guanidino: methylation/citrullination;
  Ser/Thr hydroxyl: phosphorylation/O-glycosylation, the "Yin-Yang"
  hypothesis). The `ptm_crosstalk_aviso` column warns when 2+ types in
  real competition pass at the same position, see
  `src/engines/ptm_annotation.py::_PTM_COMPETITION_GROUPS`.
- **Kinase specificity for phosphorylation.** Neither DeepMVP nor DeepPTMPred
  distinguishes WHICH kinase family phosphorylates a site; both predict
  "phosphorylatable in general". Unlike the type of
  polyubiquitin chain (K48 = proteasomal degradation; K63 = signalling/
  repair, NOT degradative, which depends on the real E2/E3 conjugating additional
  ubiquitins inside the cell, a later cellular event not observable
  in sequence/structure), kinase specificity IS a local
  sequence property, and a real published source exists:
  [Kinase Library](https://kinase-library.phosphosite.org) (Johnson et al.
  2023 *Nature*, 303 Ser/Thr kinases + Yaron-Barir et al. 2024 *Nature*,
  the complete Tyr kinome). For each phosphorylation site that the consensus
  accepts, the columns `kinase_library_top_kinase`/
  `kinase_library_top_family`/`kinase_library_percentile`/
  `kinase_library_top3_kinases` (purely informative, they never decide
  `pasa_umbral`/`consenso`) report the most probable kinase/family according to
  the published specificity matrices. Verified 2026-08-07 against a
  real site (p53 S33): top hit ATM, which matches the real literature on
  the DNA damage response.
- **Known false positive in the panel's negative control (N-97 of
  `kit_ligand_scf`/SCF, 2026-08-11).** UniProt (PMID 1381905) confirms that
  this valid N-Y-S sequon is NEVER glycosylated in either of the 2 isoforms
  of SCF. It is the only real negative control of the validation panel. With
  the 2-engine consensus (DeepMVP+EMNGly, after removing StackGlyEmbed on
  2026-08-10) both engines accept it with `consenso=True`
  (DeepMVP 0.939, EMNGly 0.792), the first time this negative control has been
  tested against the current architecture (EMNGly did not exist in the pipeline on
  2026-08-04, when this same site was correctly rejected with an earlier
  version without EMNGly). A conscious decision: the
  `Settings.EMNGLY_MIN_PROBABILITY=0.5` threshold is NOT adjusted for this single case
  (it would move a threshold already validated against the 301 real sites of
  N-GlyDE, MCC 0.82, overfitting to n=1), nor is the removal of StackGlyEmbed
  reverted (real Scipion friction already resolved by that decision). It is documented
  as a known limitation of the predictor (MCC 0.82 implies that SOME false
  positive is to be expected, and this is precisely the type of site designed
  to be difficult: chemically indistinguishable from N-90/N-145 except by
  something that neither the sequence nor the static structure captures),
  and additionally covered by the general disclaimer of "predicted
  capability, not observed occurrence" above. Detail in STATUS.md.

## Current state (2026-08-10)

The complete end-to-end pipeline (Phase 1/1.5/2/3/3b) is implemented, installed and
verified with REAL runs (not merely planned) on this machine: DeepMVP
(weights + real calibration), DeepPTMPred (PyRosetta + ESM-2 + real calibration
of the 17 types, 2 train/inference bugs found and fixed,
phi/psi and plDDT, see STATUS.md) and MeToken. Confirmed feature-complete and without
Phase A/3c (real structural modelling via PyRosetta, removed from scope) or
StackGlyEmbed (it depended on a sibling project's venv, real friction
for the Scipion integration, also removed from scope) following
feedback from Carlos, see STATUS.md.

**Decision 2 (second consensus engine for N-glycosylation) CLOSED**:
CoNglyPred (the original candidate) confirmed as having no published weights anywhere,
replaced by **EMNGly** (real weights, verified at byte level),
implemented and wired into the `n_linked_glycosylation` consensus together
with DeepMVP. See STATUS.md "Decision 2" for the complete
detail, including the 2 go/no-go checks (both PASSED, 2026-08-07): the
alignment of `structure_emb` was verified against real GlyGen sites
in a PDB with gaps (Alpha-1-antitrypsin, 1QLP), and the published MCC
(0.736) was reproduced, and in fact exceeded, at 0.8197, by running the
real pipeline over the independent N-GlyDE set. The
`EMNGLY_MIN_PROBABILITY=0.5` threshold is no longer provisional.

## Usage

```bash
python pipeline.py --input inputs/mi_proteina.fasta
python pipeline.py --input inputs/mi_estructura.pdb   # PDB path: engines + consensus + informative cross-references
```

## Installation

```bash
pip install -r requirements.txt
pytest tests/
```

DeepMVP and DeepPTMPred require their own dedicated venvs/conda envs
(mutually incompatible stacks), never the pipeline's main venv:

```bash
git clone https://github.com/bzhanglab/DeepMVP
conda create -n deepmvp python=3.7.10 -y
conda run -n deepmvp pip install -r DeepMVP/requirements.txt
export DEEPMVP_PYTHON_BIN=$(conda run -n deepmvp which python)
# Weights: manual download from https://deepmvp.ptmax.org/ -> DeepMVP/models/

git clone https://github.com/kuikui-wang/DeepPTMPred
conda create -n deepptmpred python=3.10 -y
# Install per DeepPTMPred/pred/train_PTM/environment.yml (TensorFlow 2.15 +
# tensorflow-addons + torch + fair-esm) + pyrosetta-installer (see the repo's
# README) + ESM-2 checkpoint (~2.5GB, plus its *-contact-regression.pt companion)
export DEEPPTMPRED_PYTHON_BIN=$(conda run -n deepptmpred which python)
```

MeToken (OPTIONAL and informative corroboration of the type, PDB path, see
STATUS.md 2026-08-01, NEVER a consensus engine, it can be omitted without
losing functionality of the main pipeline; `Settings.METOKEN_ENABLED`
degrades with only a warning if it is not installed):

```bash
git clone https://github.com/A4Bio/MeToken
conda create -n metoken python=3.10 -y
conda run -n metoken pip install --index-url https://download.pytorch.org/whl/cpu torch
conda run -n metoken pip install numpy scipy biopython transformers omegaconf tqdm pandas huggingface-hub h5py
conda run -n metoken pip install --no-build-isolation torch_scatter  # no prebuilt wheel, compiles from source
export METOKEN_PYTHON_BIN=$(conda run -n metoken which python)
# Weights: download release 1.0 (~88MB) and unpack into MeToken/pretrained_model/
curl -L -o /tmp/pretrained_model.zip https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip
python -c "import zipfile; zipfile.ZipFile('/tmp/pretrained_model.zip').extractall('MeToken')"
```

EMNGly (the real consensus engine for `n_linked_glycosylation`, PDB path
only, replacing CoNglyPred, decision 2026-08-06, see STATUS.md
"Decision 2": CoNglyPred confirmed as having no published weights anywhere).
Optional (`Settings.EMNGLY_ENABLED` degrades with only a warning; the
N-glycosylation consensus falls back to DeepMVP alone, with no consensus possible, if it
is not available):

```bash
git clone https://github.com/StellaHxy/EMNgly
python3 -m venv .venv-emngly
.venv-emngly/bin/pip install fair-esm torch "scikit-learn==1.1.1" scipy pandas numpy tqdm wget
# 'wget' (the pip package, not the CLI binary): a real transitive dependency of
# MIF/sequence_models/trRosetta_utils.py (imported in a chain by
# pretrained.py even though this project never uses that code path). Without
# it the MIF import fails with ModuleNotFoundError (confirmed real
# 2026-08-07, not documented in EMNgly's environment.yml).
# pip resolves numpy to the latest 2.x by default, incompatible at runtime with
# the compiled wheel of scikit-learn==1.1.1 (ValueError: numpy.dtype size
# changed, confirmed real 2026-08-07). Reinstall afterwards, pinned:
.venv-emngly/bin/pip install "numpy==1.23.5"
export EMNGLY_PYTHON_BIN=$(pwd)/.venv-emngly/bin/python

# MIF (structural embedding, mif.pt) already comes bundled in the clone, with no
# separate download.

# ESM-1b (~7.4GB) + contact regression companion (fair-esm ALWAYS looks
# for it, the same mechanism as the DeepPTMPred ESM-2 checkpoint above):
mkdir -p EMNgly/esm/checkpoints
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S.pt \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt \
  https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt

# SVM (verified 2026-08-06 with a real HTTP Range, see STATUS.md): chose
# N-GlyDE.pickle (negatives restricted to the sequon, the correct regime),
# not N-GlyAltas_classifier.pkl (same folder, trained on the other
# benchmark).
mkdir -p EMNgly/checkpoints
curl -L -o EMNgly/checkpoints/N-GlyDE.pickle \
  "https://drive.usercontent.google.com/download?id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t"
```

Kinase Library (OPTIONAL and informative corroboration of kinase
specificity for phosphorylation, BOTH paths, FASTA and PDB, biological
coherence analysis 2026-08-07 point 5, see the "Scope and
interpretation" section above, NEVER a consensus engine;
`Settings.KINASE_LIBRARY_ENABLED` degrades with only a warning if it is not
installed). A DEDICATED environment: the `numpy~=1.26.4`/`pandas~=2.2.3` that
the package itself pins are incompatible with the versions pinned in
this pipeline's requirements.txt:

```bash
conda create -n kinase_library python=3.10 -y
conda run -n kinase_library pip install kinase-library
export KINASE_LIBRARY_PYTHON_BIN=$(conda run -n kinase_library which python)
```

> The SVM weights were generated with `scikit-learn==1.1.1`. The repo's
> `environment.yml` pins 1.5.1, but loading the pickles with a
> version different from the training one triggers `InconsistentVersionWarning`
> (not guaranteed). Pin `scikit-learn==1.1.1` in the venv, not what
> the repo's `environment.yml` asks for.

See `STATUS.md` for the complete detail of what has already been verified (all the
repos/envs cloned and tested on this machine, plus EMNGly's 2 go/no-go
checks, both PASSED 2026-08-07) and what is still missing.

## Licences

This repository (own code, see `LICENSE`): CC BY-NC 4.0, non-commercial
use, consistent with the most restrictive real dependency (see below).

- **DeepMVP**: GPL-3.0 (declared in the original repo's `DeepMVP/LICENSE`, corrected
  2026-08-03; an earlier version of this README said MIT by mistake, never
  verified against the real file until now).
- **DeepPTMPred**: the repo (`github.com/kuikui-wang/DeepPTMPred`, vendored
  in `DeepPTMPred/`) does not declare its own licence, but the associated paper is
  CC BY-NC 4.0 (Oxford University Press). The corresponding authors
  (Yong Liu, Junwen Wang) were contacted directly to confirm that the
  code follows the same terms. **Junwen Wang confirmed by email on
  2026-07-29**: *"I confirm that the GitHub code follows the same CC BY-NC
  terms."* Non-commercial use (research/TFG + future integration as a
  Scipion plugin by the CNB, a public institution) fits without problem
  within CC BY-NC. Complete detail in `STATUS.md`.
- **MeToken**: MIT (declared in the original repo's `MeToken/LICENSE`).
- **EMNGly**: **MIT** (resolved 2026-08-11). The `StellaHxy/EMNgly` repo did not
  declare a LICENSE; the corresponding authors (Yaojun
  Wang, `wangyaojun@cau.edu.cn`, China Agricultural University; Shiwei Sun,
  `dwsun@ict.ac.cn`, ICT-CAS) were contacted on 2026-08-10. They replied confirming that
  they had added the MIT licence to the GitHub repo and explicitly authorizing
  integration into Scipion. The MIT LICENSE was verified in commit `824d6dc`
  (copyright Xiaoyang Hou, Yu Wang) and synchronized in the local vendored copy
  (`EMNgly/`, fast-forward `0b3c84d..a719532`). The paper remains CC BY
  4.0 (PMC10627407). MIF (vendored inside EMNgly, `model/MIF/`) is from
  Microsoft (`microsoft/protein-sequence-models`), a permissive BSD-2 licence
  verified in the official repo. EMNgly's copy lost its LICENSE when
  vendoring it, but it is covered all the same by the MIT licence just added to the
  repo that contains it.
- **PyRosetta** (the per-residue SASA feature inside DeepPTMPred): an
  academic/non-commercial RosettaCommons licence, already covered by this
  project's research/TFG use.

## Architecture decisions

Documented in the vault, not in this repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
