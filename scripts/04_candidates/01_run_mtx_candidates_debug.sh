#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIG="${ROOT}/inputs/01_cgmap/MTX.pdb"
POLAR="${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json"
FRAMES="${ROOT}/outputs/02_polar_sites/MTX_site_frames.json"
VDX="${ROOT}/processed/03_vdxform"

OUT_PREFIX="${ROOT}/processed/04_candidates/MTX_candidates_debug"
LOG="${ROOT}/logs/04_candidates/MTX_candidates_debug.log"

mkdir -p "$(dirname "$OUT_PREFIX")" "$(dirname "$LOG")"

MIN_CA_LIG_DIST="${MIN_CA_LIG_DIST:-4.0}"
MAX_CA_LIG_DIST="${MAX_CA_LIG_DIST:-12.0}"
SCORE_W_FACING="${SCORE_W_FACING:-0.25}"
SCORE_W_CENTROID_FACING="${SCORE_W_CENTROID_FACING:-0.15}"
SCORE_W_CA_SHELL="${SCORE_W_CA_SHELL:-0.20}"
CA_SHELL_CENTER="${CA_SHELL_CENTER:-8.0}"
CA_SHELL_WIDTH="${CA_SHELL_WIDTH:-2.0}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] MIN_CA_LIG_DIST=$MIN_CA_LIG_DIST MAX_CA_LIG_DIST=$MAX_CA_LIG_DIST"
  echo "[run] SCORE_W_FACING=$SCORE_W_FACING SCORE_W_CENTROID_FACING=$SCORE_W_CENTROID_FACING SCORE_W_CA_SHELL=$SCORE_W_CA_SHELL"
  echo "[run] CA_SHELL_CENTER=$CA_SHELL_CENTER CA_SHELL_WIDTH=$CA_SHELL_WIDTH"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/04_candidates/01_generate_candidates.py" \
    --ligand-pdb "$LIG" \
    --polar-sites "$POLAR" \
    --site-frames "$FRAMES" \
    --vdxform-dir "$VDX" \
    --out-prefix "$OUT_PREFIX" \
    --chunk-size 5000 \
    --top-per-site 200 \
    --clash-tol 0.5 \
    --min-ca-lig-centroid-dist "$MIN_CA_LIG_DIST" \
    --max-ca-lig-centroid-dist "$MAX_CA_LIG_DIST" \
    --score-w-facing "$SCORE_W_FACING" \
    --score-w-centroid-facing "$SCORE_W_CENTROID_FACING" \
    --score-w-ca-shell "$SCORE_W_CA_SHELL" \
    --ca-shell-center "$CA_SHELL_CENTER" \
    --ca-shell-width "$CA_SHELL_WIDTH"
} 2>&1 | tee "$LOG"
