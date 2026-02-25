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

SCORE_W_PRIOR="${SCORE_W_PRIOR:-1.0}"
SCORE_W_COVERAGE="${SCORE_W_COVERAGE:-0.05}"
SCORE_W_CONTACT="${SCORE_W_CONTACT:-0.15}"
SCORE_W_SHELL="${SCORE_W_SHELL:-0.3}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] SCORE_W_PRIOR=$SCORE_W_PRIOR SCORE_W_COVERAGE=$SCORE_W_COVERAGE SCORE_W_CONTACT=$SCORE_W_CONTACT SCORE_W_SHELL=$SCORE_W_SHELL"
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
    --score-w-prior "$SCORE_W_PRIOR" \
    --score-w-coverage "$SCORE_W_COVERAGE" \
    --score-w-contact "$SCORE_W_CONTACT" \
    --score-w-shell "$SCORE_W_SHELL"
} 2>&1 | tee "$LOG"
