# MTX Pocket-Contact Sweep (2026-02-25)

## Goal
Evaluate whether the new pocket-contact knobs improve **PLIP-grounded** ligand polar satisfaction for MTX on full vdXform libraries.

Swept knobs:
- `MIN_SC_CONTACT_COUNT`: `1, 2`
- `SC_CONTACT_WEIGHT`: `0.05, 0.10, 0.15`

Fixed settings:
- `solver=greedy`
- `top_per_site=200`
- `top_per_site_per_atom=50`
- `acceptor_model=plip`
- `min/max_sidechain_contact_dist=2.8/4.8`
- `vdxform_dir=processed/03_vdxform_full`

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
  --min-sidechain-contact-count-list 1,2 \
  --sidechain-contact-weight-list 0.05,0.10,0.15 \
  --tag pocket_contact_grid_20260225_v2
```

Raw outputs:
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_20260225_v2/summary.json`
- `processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_20260225_v2/summary.md`

## Results
| run_id | cnt | weight | n_candidates | n_selected | simple_sat | plip_sat | clash_ok | internal_ok | total_s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| run_01_cnt1_w0p05 | 1 | 0.05 | 1408 | 8 | 1 | 0 | 1 | 1 | 129.68 |
| run_05_cnt2_w0p1 | 2 | 0.10 | 1408 | 8 | 1 | 0 | 1 | 1 | 130.53 |
| run_06_cnt2_w0p15 | 2 | 0.15 | 1408 | 8 | 1 | 0 | 1 | 1 | 130.67 |
| run_02_cnt1_w0p1 | 1 | 0.10 | 1408 | 8 | 1 | 0 | 1 | 1 | 131.54 |
| run_04_cnt2_w0p05 | 2 | 0.05 | 1408 | 8 | 1 | 0 | 1 | 1 | 132.35 |
| run_03_cnt1_w0p15 | 1 | 0.15 | 1408 | 8 | 1 | 0 | 1 | 1 | 133.07 |

Aggregate:
- `plip_all_satisfied_rate = 0/6 = 0.0`
- `feasible_count = 0`
- PLIP unsatisfied atoms are stable across runs: `N`, `N5`, `O`

## Interpretation
- Pocket-contact knobs changed candidate ranking/metadata but did **not** change final motif quality under this search budget (`top_per_site=200`).
- Current bottleneck is not clash/coverage (all pass), but PLIP-level interaction completeness.
- Next step should target candidate diversity/search space or interaction typing for `N`, `N5`, `O` specifically.
