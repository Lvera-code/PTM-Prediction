"""External data clients used by Phase 3b (informative cross-references).

Phase 3c (real structural modelling via PyRosetta) was removed from scope
2026-08-10 (project supervision feedback), so this package contains only
``uniprot_localization_client.py``, queried from
``src/engines/ptm_annotation.py`` (the main process, with no dependency on
pyrosetta) for the secretory pathway corroboration in N-glycosylation.
"""
