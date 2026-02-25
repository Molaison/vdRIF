#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="${ROOT}/logs/99_harness/mtx_balanced_solver_sweep.log"
mkdir -p "$(dirname "$LOG")"

CANDIDATES_NPZ="${CANDIDATES_NPZ:-processed/04_candidates/MTX_candidates.npz}"
CANDIDATES_META="${CANDIDATES_META:-processed/04_candidates/MTX_candidates.json}"
LIGAND_PDB="${LIGAND_PDB:-inputs/01_cgmap/MTX.pdb}"

{
  echo "[run] repo: $ROOT"
  echo "[run] candidates_npz: $CANDIDATES_NPZ"
  echo "[run] candidates_meta: $CANDIDATES_META"
  echo "[run] ligand_pdb: $LIGAND_PDB"
  uv run --no-project -p 3.11 --with numpy --with ortools python \
    "${ROOT}/scripts/99_harness/06_mtx_balanced_solver_sweep.py" \
    --repo-root "$ROOT" \
    --candidates-npz "$CANDIDATES_NPZ" \
    --candidates-meta "$CANDIDATES_META" \
    --ligand-pdb "$LIGAND_PDB" \
    "${@}"
} 2>&1 | tee "$LOG"
