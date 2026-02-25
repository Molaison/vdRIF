#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIG="${ROOT}/inputs/01_cgmap/MTX.pdb"
POLAR="${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json"
FRAMES="${ROOT}/outputs/02_polar_sites/MTX_site_frames.json"
VDX="${ROOT}/processed/03_vdxform"

OUT_PREFIX="${ROOT}/processed/04_candidates/MTX_candidates"
LOG="${ROOT}/logs/04_candidates/MTX_candidates.log"

mkdir -p "$(dirname "$OUT_PREFIX")" "$(dirname "$LOG")"

TOP_PER_SITE="${TOP_PER_SITE:-2000}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
CLASH_TOL="${CLASH_TOL:-0.5}"
MIN_SC_CONTACT_DIST="${MIN_SC_CONTACT_DIST:-2.8}"
MAX_SC_CONTACT_DIST="${MAX_SC_CONTACT_DIST:-4.8}"
MIN_SC_CONTACT_COUNT="${MIN_SC_CONTACT_COUNT:-1}"
SC_CONTACT_WEIGHT="${SC_CONTACT_WEIGHT:-0.05}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] TOP_PER_SITE=$TOP_PER_SITE CHUNK_SIZE=$CHUNK_SIZE CLASH_TOL=$CLASH_TOL"
  echo "[run] sidechain contact window: [$MIN_SC_CONTACT_DIST, $MAX_SC_CONTACT_DIST]"
  echo "[run] MIN_SC_CONTACT_COUNT=$MIN_SC_CONTACT_COUNT SC_CONTACT_WEIGHT=$SC_CONTACT_WEIGHT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/04_candidates/01_generate_candidates.py" \
    --ligand-pdb "$LIG" \
    --polar-sites "$POLAR" \
    --site-frames "$FRAMES" \
    --vdxform-dir "$VDX" \
    --out-prefix "$OUT_PREFIX" \
    --chunk-size "$CHUNK_SIZE" \
    --top-per-site "$TOP_PER_SITE" \
    --clash-tol "$CLASH_TOL" \
    --min-sidechain-contact-dist "$MIN_SC_CONTACT_DIST" \
    --max-sidechain-contact-dist "$MAX_SC_CONTACT_DIST" \
    --min-sidechain-contact-count "$MIN_SC_CONTACT_COUNT" \
    --sidechain-contact-weight "$SC_CONTACT_WEIGHT" \
    --require-full-coverage
} 2>&1 | tee "$LOG"
