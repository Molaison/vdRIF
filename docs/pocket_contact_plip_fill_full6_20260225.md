# MTX Pocket-Contact + PLIP-Fill Full 6-Point Sweep (2026-02-26)

## Goal
Re-run the original 6-point MTX pocket-contact grid with the new PLIP-guided rescue loop and updated candidate retention defaults.

## Command
```bash
python scripts/99_harness/05_benchmark_pocket_contact_grid.py \
  --python-bin /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.venv/bin/python \
  --script-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.worktrees/feat-vdm-pocket-remodel-20260225 \
  --data-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm \
  --vdxform-dir processed/03_vdxform_full \
  --solver greedy \
  --top-per-site 200 \
  --chunk-size 5000 \
  --time-limit-s 120 \
  --clash-tol 0.5 \
  --acceptor-model plip \
  --min-sidechain-contact-count-list 1,2 \
  --sidechain-contact-weight-list 0.05,0.10,0.15 \
  --enable-plip-fill \
  --plip-fill-max-res 15 \
  --plip-fill-top-try-per-atom 120 \
  --tag pocket_contact_grid_plipfill_full6_20260225
```

## Outputs
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_plipfill_full6_20260225/summary.json`
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_plipfill_full6_20260225/summary.md`

## Results
- `plip_all_satisfied_rate = 6/6 = 1.0`
- All runs are feasible:
  - `coverage_complete=true`
  - `ligand_clash_ok=true`
  - `internal_clash_ok=true`
- Across all 6 runs:
  - pre-fill PLIP `n_satisfied=6`, unsatisfied atoms: `N/N5/O`
  - post-fill PLIP `n_satisfied=9`, unsatisfied atoms: none

## Interpretation
- The previous failure mode (`0/6` PLIP all-satisfied) is fixed under the new rescue path.
- The critical change is not only pocket-contact weighting; it is the PLIP-guided post-solver augmentation combined with higher per-atom candidate retention.
