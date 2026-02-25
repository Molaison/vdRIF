# MTX Pocket-Contact + PLIP-Fill Smoke (2026-02-25)

## Goal
Verify that the new PLIP-guided post-solver rescue loop is wired end-to-end in the pocket-contact sweep harness.

## Command
```bash
python scripts/99_harness/05_benchmark_pocket_contact_grid.py \
  --python-bin /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.venv/bin/python \
  --script-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.worktrees/feat-vdm-pocket-remodel-20260225 \
  --data-root /xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm \
  --vdxform-dir processed/03_vdxform_full \
  --solver greedy \
  --top-per-site 200 \
  --top-per-site-per-atom 50 \
  --chunk-size 5000 \
  --time-limit-s 120 \
  --clash-tol 0.5 \
  --acceptor-model plip \
  --min-sidechain-contact-count-list 1 \
  --sidechain-contact-weight-list 0.05 \
  --enable-plip-fill \
  --plip-fill-max-res 15 \
  --plip-fill-top-try-per-atom 120 \
  --tag pocket_contact_grid_plipfill_smoke_20260225
```

## Output
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_plipfill_smoke_20260225/summary.json`
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_plipfill_smoke_20260225/summary.md`

## Key result
- Baseline PLIP: `n_satisfied = 6/9`, unsatisfied: `N, N5, O`
- Post-fill PLIP: `n_satisfied = 7/9`, unsatisfied: `N5, O`
- `plip_fill_used_for_metrics = true`
- Ligand clash and motif-internal clash remain clean (`ok=true`, `worst_overlap=0`)

## Interpretation
- The PLIP rescue loop is functional and recovers one missing atom (`N`) deterministically.
- Remaining blockers are concentrated at `N5` and `O`; these likely require improved chemistry typing and/or broader candidate search space rather than contact-weight tuning alone.
