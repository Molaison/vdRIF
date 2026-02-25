#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIG="${ROOT}/inputs/01_cgmap/MTX.pdb"
CAND_NPZ="${ROOT}/processed/04_candidates/MTX_candidates.npz"
CAND_META="${ROOT}/processed/04_candidates/MTX_candidates.json"

OUT_JSON="${ROOT}/outputs/05_solver/MTX_motif.json"
OUT_PDB="${ROOT}/outputs/05_solver/MTX_motif.pdb"
VAL_JSON="${ROOT}/outputs/05_solver/MTX_motif_validation.json"
CLASH_JSON="${ROOT}/outputs/05_solver/MTX_motif_clash_validation.json"
INTERNAL_CLASH_JSON="${ROOT}/outputs/05_solver/MTX_motif_internal_clash_validation.json"
LOG="${ROOT}/logs/05_solver/MTX_motif.log"

mkdir -p "$(dirname "$OUT_JSON")" "$(dirname "$OUT_PDB")" "$(dirname "$LOG")"

TIME_LIMIT_S="${TIME_LIMIT_S:-300}"
CLASH_TOL="${CLASH_TOL:-0.5}"
POCKET_CONTACT_SCORE_WEIGHT="${POCKET_CONTACT_SCORE_WEIGHT:-0.15}"

{
  echo "[run] ligand: $LIG"
  echo "[run] candidates: $CAND_NPZ"
  echo "[run] candidates_meta: $CAND_META"
  echo "[run] out_json: $OUT_JSON"
  echo "[run] out_pdb: $OUT_PDB"
  echo "[run] TIME_LIMIT_S=$TIME_LIMIT_S CLASH_TOL=$CLASH_TOL POCKET_CONTACT_SCORE_WEIGHT=$POCKET_CONTACT_SCORE_WEIGHT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/05_solver/01_solve_motif.py" \
    --candidates-npz "$CAND_NPZ" \
    --candidates-meta "$CAND_META" \
    --ligand-pdb "$LIG" \
    --out-json "$OUT_JSON" \
    --out-pdb "$OUT_PDB" \
    --min-res 8 \
    --max-res 15 \
    --time-limit-s "$TIME_LIMIT_S" \
    --num-workers 1 \
    --grid-size 4.0 \
    --ca-prefilter 8.0 \
    --clash-tol "$CLASH_TOL" \
    --pocket-contact-score-weight "$POCKET_CONTACT_SCORE_WEIGHT"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/03_validate_motif_polar_satisfaction.py" \
    --polar-sites "${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json" \
    --motif-pdb "$OUT_PDB" \
    -o "$VAL_JSON"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/04_validate_motif_clashes.py" \
    --motif-pdb "$OUT_PDB" \
    -o "$CLASH_JSON" \
    --tol "$CLASH_TOL" \
    --fail-overlap 0.5

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/05_validate_motif_internal_clashes.py" \
    --motif-pdb "$OUT_PDB" \
    -o "$INTERNAL_CLASH_JSON" \
    --tol "$CLASH_TOL" \
    --fail-overlap 0.5
} 2>&1 | tee "$LOG"
