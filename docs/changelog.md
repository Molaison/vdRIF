# Changelog

## 2026-02-25

- Hardened candidate chemistry defaults in `scripts/04_candidates/01_generate_candidates.py`:
  - default is now **sidechain-only** satisfaction (backbone H-bonds require explicit `--allow-backbone-hbonds`),
  - default acceptor typing switched to **PLIP-style** (`--acceptor-model plip`),
  - OpenBabel donor-H failure now falls back to RDKit explicit-H coordinates instead of dropping donor-angle geometry entirely.
- Added pocket-completeness scoring during candidate generation:
  - per-candidate sidechain-ligand contact count (`lig_sc_contact_count_u8`) is computed and stored,
  - new controls: `--pocket-contact-cutoff`, `--min-pocket-sidechain-contacts`, `--pocket-contact-weight`,
  - candidate score now includes a weighted contact term to prefer denser, ligand-facing pockets.
- Updated MTX candidate runner scripts (`scripts/04_candidates/01_run_mtx_candidates_debug.sh`, `scripts/04_candidates/02_run_mtx_candidates.sh`) to use strict PLIP + pocket-contact settings by default.
  - Current MTX default run profile explicitly sets `--allow-backbone-hbonds` and `--min-lig-donor-angle-deg 60` to avoid false-negative coverage collapse on ligand donor `N` with the present local vdM library.
- Aligned polar validation defaults with strict typing (`scripts/05_solver/03_validate_motif_polar_satisfaction.py`, solver run wrappers) by using `--acceptor-model plip`.
- Added MTX pocket-knob sweep harness:
  - `scripts/99_harness/05_mtx_pocket_knob_sweep.py`
  - `scripts/99_harness/05_run_mtx_pocket_knob_sweep.sh`
  - Wrapper resolves data and Python from git common root, so worktree runs can consume shared full vdXform assets without `uv sync`.
- Ran full-library local sweeps and recorded reproducible evidence:
  - `processed/99_harness/mtx_pocket_knob_sweep_local_full4_20260225-211454/summary.json`
  - `processed/99_harness/mtx_pocket_knob_sweep_local_ang75_20260225-212315/summary.json`
  - In both sweeps, `--pocket-contact-weight 0.0` produced severe motif-internal overlaps (`worst_internal_overlap=1.709`), while `0.2` produced feasible motifs with full polar satisfaction and clash-clean outputs.
  - `--min-pocket-sidechain-contacts` (`0` vs `1`) had minor impact once contact weighting was enabled; both were feasible.
- Recommended default profile for MTX pocket realism:
  - `--allow-backbone-hbonds`
  - `--acceptor-model plip`
  - `--min-lig-donor-angle-deg 75` (or `60`; both passed in local sweep)
  - `--pocket-contact-weight 0.2`
  - `--min-pocket-sidechain-contacts 1`

## 2026-02-02

- Extended symmetric ligand site-frame handling beyond `coo` to aromatic `ph`/`phenol` (`CD1/CD2` swap frames) in `scripts/02_polar_sites/03_build_ligand_site_frames.py`.
- Added a focused regression harness for symmetry swap frames: `scripts/99_harness/02_run_symmetry_site_frames_regression.sh`.
- Documented quantitative per-voxel storage impact for `Rot12ScoreSat96` vs `Rot10Score6Sat16` in `docs/irot_compression_plan.md`.
- Expanded rifdock integration pointers and a concrete “Tier 1 placement importer” proposal in `docs/rifgen_integration_notes.md`.
- Added MTX rif-export smoke harness: `scripts/99_harness/04_run_mtx_rif_export_smoke.sh`.
- Made `scripts/06_rif_export/01_export_rif_inputs.py` default score mapping rifdock-insert-safe (order-preserving shift to `<=0`).

## 2026-01-13

- Built **full-scale** vdXform libraries for MTX-required CGs under `processed/03_vdxform_full/`:
  - `coo` (348,142), `conh2` (303,932), `ccn` (145,248), `ph` (417,170), `bb_cnh` (314,744) entries.
- Refactored candidate generation to be tractable on 1-core:
  - Prunes by **geometric satisfaction + sidechain-facing** before running expensive full-atom ligand clash checks.
  - Uses a deterministic `heapq` top-K per site (stable tie-breaking by `vdm_id`).
- Added a deterministic `--solver=greedy` option to scale motif selection to ~25k candidates (avoids explicit pairwise conflict graphs); CP-SAT remains default.
- Verified end-to-end determinism on full libs: `processed/99_harness/mtx_det_full_vdxform_det_top200/report.json`.
- Verified determinism at large candidate counts with greedy solver: `processed/99_harness/mtx_det_full_vdxform_det_top2000_greedy/report.json`.
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
