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
INTERNAL_CLASH_JSON="${ROOT}/outputs/05_solver/MTX_motif_debug_internal_clash_validation.json"
LOG="${ROOT}/logs/05_solver/MTX_motif_debug.log"

mkdir -p "$(dirname "$OUT_JSON")" "$(dirname "$OUT_PDB")" "$(dirname "$LOG")"

OBJECTIVE_MODE="${OBJECTIVE_MODE:-balanced}"
TARGET_RES="${TARGET_RES:-12}"
TARGET_RES_PENALTY="${TARGET_RES_PENALTY:-1200}"
SITE_DIVERSITY_REWARD="${SITE_DIVERSITY_REWARD:-2000}"
MIN_UNIQUE_SITES="${MIN_UNIQUE_SITES:-6}"
MAX_PER_SITE="${MAX_PER_SITE:-2}"
CA_PREFILTER="${CA_PREFILTER:-12.0}"
CLASH_TOL="${CLASH_TOL:-0.5}"

{
  echo "[run] ligand: $LIG"
  echo "[run] candidates: $CAND_NPZ"
  echo "[run] candidates_meta: $CAND_META"
  echo "[run] out_json: $OUT_JSON"
  echo "[run] out_pdb: $OUT_PDB"
  echo "[run] OBJECTIVE_MODE=$OBJECTIVE_MODE TARGET_RES=$TARGET_RES TARGET_RES_PENALTY=$TARGET_RES_PENALTY"
  echo "[run] SITE_DIVERSITY_REWARD=$SITE_DIVERSITY_REWARD MIN_UNIQUE_SITES=$MIN_UNIQUE_SITES MAX_PER_SITE=$MAX_PER_SITE"
  echo "[run] CA_PREFILTER=$CA_PREFILTER CLASH_TOL=$CLASH_TOL"
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
    --ca-prefilter "$CA_PREFILTER" \
    --clash-tol "$CLASH_TOL" \
    --objective-mode "$OBJECTIVE_MODE" \
    --target-res "$TARGET_RES" \
    --target-res-penalty "$TARGET_RES_PENALTY" \
    --site-diversity-reward "$SITE_DIVERSITY_REWARD" \
    --min-unique-sites "$MIN_UNIQUE_SITES" \
    --max-per-site "$MAX_PER_SITE"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/03_validate_motif_polar_satisfaction.py" \
    --polar-sites "${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json" \
    --motif-pdb "$OUT_PDB" \
    -o "$VAL_JSON"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/04_validate_motif_clashes.py" \
    --motif-pdb "$OUT_PDB" \
    -o "$CLASH_JSON" \
    --tol "$CLASH_TOL" \
    --fail-overlap "$CLASH_TOL"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/05_validate_motif_internal_clashes.py" \
    --motif-pdb "$OUT_PDB" \
    -o "$INTERNAL_CLASH_JSON" \
    --tol "$CLASH_TOL" \
    --fail-overlap "$CLASH_TOL" \
    --ligand-resname MTX
} 2>&1 | tee "$LOG"
