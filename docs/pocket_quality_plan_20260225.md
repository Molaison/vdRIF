# Pocket Quality Plan (2026-02-25)

## Problem
`vdM`-sampled placements can satisfy polar constraints but still form loose or incomplete pockets:
- sidechain points roughly toward ligand, but does not make enough sidechain-heavy-atom contacts;
- solver then selects from a candidate pool with weak pocket geometry signal.

## Divergent Ideas

1. Hard geometric gating only
- Increase strict filters (distance/angle/facing) until pockets look dense.
- Risk: candidate collapse and unsatisfied polar atoms.

2. Solver-only objective rewrite
- Keep candidate generation unchanged and add pocket compactness objective in CP-SAT/greedy.
- Risk: low-quality candidate pool remains the bottleneck.

3. Candidate-stage pocket contact signal (selected)
- Add sidechain-ligand contact-shell features and integrate into candidate filtering/ranking.
- Keep deterministic behavior and existing coverage/clash logic.
- Benefit: removes detached placements early and gives solver a stronger pool.

## Selected Strategy
Implement idea 3 with two controls:
- hard gate: minimum number of sidechain atoms in a ligand-contact shell;
- soft score bonus: prefer more sidechain contacts while preserving existing prior/coverage scoring.

## Implemented Changes
- `scripts/04_candidates/01_generate_candidates.py`
  - Added `_sidechain_contact_features(...)`.
  - Added CLI knobs:
    - `--min-sidechain-contact-dist`
    - `--max-sidechain-contact-dist`
    - `--min-sidechain-contact-count`
    - `--sidechain-contact-weight`
  - Integrated contact gating + contact-weighted ranking before top-K retention.
- `scripts/04_candidates/01_run_mtx_candidates_debug.sh`
- `scripts/04_candidates/02_run_mtx_candidates.sh`
  - Added env-tunable defaults for the new knobs.
- `scripts/99_harness/01_mtx_determinism_regression.py`
  - Added harness args and propagated new knobs into candidate generation.

## Suggested Tuning Grid
Use full vdxform libraries and sweep:
- `MIN_SC_CONTACT_COUNT`: `1, 2`
- `SC_CONTACT_WEIGHT`: `0.05, 0.10, 0.15`
- `MIN_SC_CONTACT_DIST`: `2.8`
- `MAX_SC_CONTACT_DIST`: `4.8`

Select settings by:
1. PLIP all-satisfied rate,
2. ligand clash/internal clash pass rate,
3. motif size stability and determinism.
