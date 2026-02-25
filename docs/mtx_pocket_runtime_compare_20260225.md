# MTX Pocket Runtime Compare (local vs dg)

- local summary: `/xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/processed/99_harness/mtx_pocket_contact_grid_pocket_contact_grid_20260225_v2/summary.json`
- dg summary: `/xcfhome/zpzeng/project_run_log/spur_of_moment/001_molten_rifdock_vdm/.worktrees/feat-vdm-pocket-repair-20260225-202728/processed/99_harness/mtx_pocket_knob_sweep_dg_grid6_20260225-223325/summary.json`
- matched runs: 6

| min_sc_contacts | weight | local_total_s | dg_total_s | speedup(local/dg) | local_candidates | dg_candidates |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.05 | 129.68 | 88.22 | 1.470 | 1408 | 1344 |
| 1 | 0.10 | 131.54 | 85.68 | 1.535 | 1408 | 1344 |
| 1 | 0.15 | 133.07 | 85.39 | 1.558 | 1408 | 1344 |
| 2 | 0.05 | 132.35 | 86.07 | 1.538 | 1408 | 1344 |
| 2 | 0.10 | 130.53 | 85.47 | 1.527 | 1408 | 1344 |
| 2 | 0.15 | 130.67 | 85.53 | 1.528 | 1408 | 1344 |

- mean local total_s: 131.31
- mean dg total_s: 86.06
- median local total_s: 131.11
- median dg total_s: 85.60
- mean speedup(local/dg): 1.526x
