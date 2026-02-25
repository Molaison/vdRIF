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

TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-32}"
MIN_POCKET_CONTACT_COUNT="${MIN_POCKET_CONTACT_COUNT:-1}"
POCKET_CONTACT_DIST="${POCKET_CONTACT_DIST:-4.5}"
POCKET_CONTACT_SCORE_WEIGHT="${POCKET_CONTACT_SCORE_WEIGHT:-0.15}"
MIN_SIDECHAIN_CENTROID_DOT="${MIN_SIDECHAIN_CENTROID_DOT:-0.05}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] TOP_PER_SITE_PER_ATOM=$TOP_PER_SITE_PER_ATOM MIN_POCKET_CONTACT_COUNT=$MIN_POCKET_CONTACT_COUNT"
  echo "[run] POCKET_CONTACT_DIST=$POCKET_CONTACT_DIST POCKET_CONTACT_SCORE_WEIGHT=$POCKET_CONTACT_SCORE_WEIGHT"
  echo "[run] MIN_SIDECHAIN_CENTROID_DOT=$MIN_SIDECHAIN_CENTROID_DOT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/04_candidates/01_generate_candidates.py" \
    --ligand-pdb "$LIG" \
    --polar-sites "$POLAR" \
    --site-frames "$FRAMES" \
    --vdxform-dir "$VDX" \
    --out-prefix "$OUT_PREFIX" \
    --chunk-size 5000 \
    --top-per-site 200 \
    --top-per-site-per-atom "$TOP_PER_SITE_PER_ATOM" \
    --clash-tol 0.5 \
    --min-sidechain-centroid-dot "$MIN_SIDECHAIN_CENTROID_DOT" \
    --min-pocket-contact-count "$MIN_POCKET_CONTACT_COUNT" \
    --pocket-contact-dist "$POCKET_CONTACT_DIST" \
    --pocket-contact-score-weight "$POCKET_CONTACT_SCORE_WEIGHT"
} 2>&1 | tee "$LOG"
