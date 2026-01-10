#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIG="${ROOT}/inputs/01_cgmap/MTX.pdb"
CAND_NPZ="${ROOT}/processed/04_candidates/MTX_candidates_debug.npz"
CAND_META="${ROOT}/processed/04_candidates/MTX_candidates_debug.json"

OUT_JSON="${ROOT}/outputs/05_solver/MTX_motif_debug.json"
OUT_PDB="${ROOT}/outputs/05_solver/MTX_motif_debug.pdb"
VAL_JSON="${ROOT}/outputs/05_solver/MTX_motif_debug_validation.json"
CLASH_JSON="${ROOT}/outputs/05_solver/MTX_motif_debug_clash_validation.json"
LOG="${ROOT}/logs/05_solver/MTX_motif_debug.log"

mkdir -p "$(dirname "$OUT_JSON")" "$(dirname "$OUT_PDB")" "$(dirname "$LOG")"

{
  echo "[run] ligand: $LIG"
  echo "[run] candidates: $CAND_NPZ"
  echo "[run] candidates_meta: $CAND_META"
  echo "[run] out_json: $OUT_JSON"
  echo "[run] out_pdb: $OUT_PDB"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/05_solver/01_solve_motif.py" \
    --candidates-npz "$CAND_NPZ" \
    --candidates-meta "$CAND_META" \
    --ligand-pdb "$LIG" \
    --out-json "$OUT_JSON" \
    --out-pdb "$OUT_PDB" \
    --min-res 8 \
    --max-res 15 \
    --time-limit-s 30 \
    --num-workers 1 \
    --grid-size 4.0 \
    --ca-prefilter 8.0 \
    --clash-tol 0.5

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/03_validate_motif_polar_satisfaction.py" \
    --polar-sites "${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json" \
    --motif-pdb "$OUT_PDB" \
    -o "$VAL_JSON"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/04_validate_motif_clashes.py" \
    --motif-pdb "$OUT_PDB" \
    -o "$CLASH_JSON" \
    --tol 0.5 \
    --fail-overlap 0.5
} 2>&1 | tee "$LOG"
