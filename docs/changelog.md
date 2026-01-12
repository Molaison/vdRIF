# Changelog

## 2026-01-13

- Built **full-scale** vdXform libraries for MTX-required CGs under `processed/03_vdxform_full/`:
  - `coo` (348,142), `conh2` (303,932), `ccn` (145,248), `ph` (417,170), `bb_cnh` (314,744) entries.
- Refactored candidate generation to be tractable on 1-core:
  - Prunes by **geometric satisfaction + sidechain-facing** before running expensive full-atom ligand clash checks.
  - Uses a deterministic `heapq` top-K per site (stable tie-breaking by `vdm_id`).
- Verified end-to-end determinism on full libs: `processed/99_harness/mtx_det_full_vdxform_det_top200/report.json`.
- Verified a real MTX motif output (8 residues) with full coverage and zero severe clashes:
  - `processed/99_harness/mtx_full_vdxform_single/MTX_motif_top200.pdb`
  - `processed/99_harness/mtx_full_vdxform_single/MTX_motif_top200_validation.json`
  - `processed/99_harness/mtx_full_vdxform_single/MTX_motif_top200_clash_validation.json`
  - `processed/99_harness/mtx_full_vdxform_single/MTX_motif_top200_internal_clash_validation.json`

## 2026-01-09

- Initialized project planning repo with bd issues for vdM × RIFGen integration.
- Vendored upstream repos for code inspection under `external/`:
  - `external/Combs2024` @ `c17476366b07e1861a94e133da5bf3a61ba5d6ee`
  - `external/rifdock` @ `62c3a1a19438f443a7bc391cce5c0574bea8b02b`
- Added initial docs scaffold: `docs/quickstart.md`, `docs/changelog.md`, `docs/takeaway.md`.
- Set up `uv` Python 3.11 env (pins `numpy<2` for RDKit compatibility; optional extra `--extra rdkit`).
- Generated MTX ligand CG mapping: `outputs/01_cgmap/MTX_cg_atommap.json`.
- Updated literature anchors for vdM and RIF (Science 2020; Nature 2022; Nat Biotech 2019) in `docs/quickstart.md` and `docs/design_plan.md`.
- Added polar-site extraction + ligand-site-frame builder (including bb_cnh fallback for uncovered donor atoms) under `scripts/02_polar_sites/`.
- Added vdXform v0 parquet→npz converter + spec: `scripts/03_vdxform/` and `docs/vdxform_spec.md`.
- Added MVP candidate generator that enumerates vdXform placements around ligand site frames with a fast ligand-clash filter and per-site top-K truncation: `scripts/04_candidates/`.
- Added deterministic motif solver (CP-SAT) with coverage + clash constraints and stable tie-breaks: `scripts/05_solver/`.
