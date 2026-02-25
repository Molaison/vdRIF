# MTX Pocket-Contact Top-Budget Smoke (2026-02-25)

## Goal
Verify whether increasing candidate diversity budget resolves residual `N/N5/O` PLIP misses when using PLIP-guided post-solver filling.

## Command
```bash
python scripts/99_harness/05_benchmark_pocket_contact_grid.py \
  --python-bin /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.venv/bin/python \
  --script-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.worktrees/feat-vdm-pocket-remodel-20260225 \
  --data-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm \
  --vdxform-dir processed/03_vdxform_full \
  --solver greedy \
  --top-per-site-list 200,1000 \
  --top-per-site-per-atom 200 \
  --chunk-size 5000 \
  --time-limit-s 120 \
  --clash-tol 0.5 \
  --acceptor-model plip \
  --min-sidechain-contact-count-list 1 \
  --sidechain-contact-weight-list 0.05 \
  --enable-plip-fill \
  --plip-fill-max-res 15 \
  --plip-fill-top-try-per-atom 120 \
  --max-runs 2 \
  --tag pocket_contact_grid_top_budget_smoke_20260225
```

## Outputs
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_top_budget_smoke_20260225/summary.json`
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_top_budget_smoke_20260225/summary.md`

## Key results
- `run_01_top200_cnt1_w0p05`: feasible `true`, PLIP `9/9` after fill.
- `run_02_top1000_cnt1_w0p05`: feasible `true`, PLIP `9/9` after fill.
- Both runs maintain:
  - `ligand_clash_ok=true`
  - `internal_clash_ok=true`
  - `worst_ligand_overlap=0`
  - `worst_internal_overlap=0`

## Interpretation
- The remaining `N/N5/O` misses are not a hard chemistry impossibility.
- With richer candidate retention (`top_per_site_per_atom=200`) plus PLIP-guided filling, the same pipeline reaches full PLIP satisfaction without introducing clashes.
