# sh4 PLIP-fill Benchmark (2026-02-25)

## Scope
Issue: `rifvdm-sh4` (Resolve residual MTX PLIP misses after plip-fill).

## Root Cause (evidence-backed)

1. MTX `N5` is typed as a ligand cation.
- Source: `/xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/outputs/02_polar_sites/MTX_polar_sites.json`
- `N5` roles: `["cation"]`

2. Repro case with residual dual misses exists.
- Input motif: `processed/99_harness/mtx_pocketfix_sweep_1772025699/baseline_legacy/MTX_motif.pdb`
- Pre-fix plip-fill report: `processed/99_harness/plipfill_sh4_mtx_baseline_1772034556/report.json`
- Residuals: `['N5', 'O1']`

3. Old plip-fill loop had target-order lock-in.
- Behavior: always picked `target = sorted(missing)[0]` (here `N5`) and aborted on first non-improving target.
- In this repro, `N5` had no improving add-only candidate:
  - 36 candidates covered `N5`
  - 34 rejected by protein clash
  - 2 non-improving
- But `O1` had improving options that were never attempted:
  - 264 candidates covered `O1`
  - 14 were improving (reduce residual set to `['N5']`)

4. PLIP interaction evidence before/after confirms O1 rescue and persistent N5 gap.
- Before XML/PDB:
  - `processed/99_harness/plipfill_sh4_mtx_baseline_1772034556/plip_MTX_motif_plipfill/MTX_motif_plipfill_report.xml`
  - `processed/99_harness/plipfill_sh4_mtx_baseline_1772034556/plip_MTX_motif_plipfill/plipfixed.MTX_motif_plipfill_kvcf3_vp.pdb`
  - ligand hit names: `['N','NA2','NA4','O','O2','OE1','OE2']`
- After XML/PDB:
  - `processed/99_harness/plipfill_sh4_mtx_baseline_after_1772035346/plip_MTX_motif_plipfill/MTX_motif_plipfill_report.xml`
  - `processed/99_harness/plipfill_sh4_mtx_baseline_after_1772035346/plip_MTX_motif_plipfill/plipfixed.MTX_motif_plipfill_2pfy0s3g.pdb`
  - ligand hit names: `['N','NA2','NA4','O','O1','O2','OE1','OE2']`

## Implementation

### A) Global-best rescue in `07_plip_fill_motif.py`
File: `scripts/05_solver/07_plip_fill_motif.py`

Changes:
- Replaced single-target loop with a global proposal pool across all currently-missing ligand atoms each round.
- Deduplicated candidate indices across targets.
- Evaluated proposals globally and selected deterministic best candidate by:
  1. max missing-count reduction
  2. max resolved-target count
  3. deterministic rank/score tie-breaks
- Added richer round diagnostics in report:
  - `missing_before`, `missing_after`, `resolved_targets`
  - `n_proposals`
  - `n_rejected_ligand_clash`, `n_rejected_protein_clash`, `n_rejected_non_improving`

### B) PLIP output variant compatibility in `06_validate_motif_plip.py`
File: `scripts/05_solver/06_validate_motif_plip.py`

Changes:
- Added PLIP output PDB discovery that supports both:
  - `plipfixed.*.pdb`
  - `*_protonated.pdb` (and `*protonated*.pdb`)
- Kept payload key `plipfixed_pdb` stable for compatibility, value now points to the discovered PLIP ligand PDB.

## Benchmarks

### MTX (residual repro case)
- Input motif: `processed/99_harness/mtx_pocketfix_sweep_1772025699/baseline_legacy/MTX_motif.pdb`
- New run report: `processed/99_harness/plipfill_sh4_mtx_baseline_after_1772035346/report.json`
- Result:
  - before: `['N5','O1']`
  - after: `['N5']`
  - at least one residual miss eliminated (`O1`), no clash regression (see below)

### MTX (single residual case)
- Input motif: `processed/99_harness/mtx_pocketfix_sweep_1772029222/pocketfix_default/MTX_motif.pdb`
- New run report: `processed/99_harness/plipfill_sh4_mtx_default_after_1772035376/report.json`
- Result: remains `['N5']` (no false improvement claim)

### Non-MTX (BZA)
- Input motif: `processed/99_harness/mtx_pocketfix_sweep_1772029101_bza/pocketfix_default/BZA_motif.pdb`
- New run report: `processed/99_harness/plipfill_sh4_bza_after_1772035705/report.json`
- Result: `all_satisfied = true`, `steps = []` (no regression)

## Clash checks (post-change outputs)

- MTX improved motif:
  - ligand clash: `processed/99_harness/plipfill_sh4_mtx_baseline_after_1772035346/MTX_motif_plipfill_clash_validation.json` -> `ok=true`, `worst_overlap=0`
  - internal clash: `processed/99_harness/plipfill_sh4_mtx_baseline_after_1772035346/MTX_motif_plipfill_internal_clash_validation.json` -> `ok=true`

- BZA motif:
  - ligand clash: `processed/99_harness/plipfill_sh4_bza_after_1772035705/BZA_motif_plipfill_clash_validation.json` -> `ok=true`, `worst_overlap=0`
  - internal clash: `processed/99_harness/plipfill_sh4_bza_after_1772035705/BZA_motif_plipfill_internal_clash_validation.json` -> `ok=true`

## Verification Commands

- `python3 -m py_compile scripts/05_solver/06_validate_motif_plip.py scripts/05_solver/07_plip_fill_motif.py`
- MTX repro + BZA runs listed above (see report paths).
