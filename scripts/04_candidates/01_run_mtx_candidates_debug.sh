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

SC_CONTACT_CUTOFF="${SC_CONTACT_CUTOFF:-4.2}"
MIN_SC_LIG_CONTACTS="${MIN_SC_LIG_CONTACTS:-2}"
MAX_SC_MIN_DIST="${MAX_SC_MIN_DIST:-4.5}"
SC_CONTACT_SCORE_W="${SC_CONTACT_SCORE_W:-0.05}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] SC_CONTACT_CUTOFF=$SC_CONTACT_CUTOFF MIN_SC_LIG_CONTACTS=$MIN_SC_LIG_CONTACTS"
  echo "[run] MAX_SC_MIN_DIST=$MAX_SC_MIN_DIST SC_CONTACT_SCORE_W=$SC_CONTACT_SCORE_W"
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
    --sidechain-contact-cutoff "$SC_CONTACT_CUTOFF" \
    --min-sidechain-ligand-contacts "$MIN_SC_LIG_CONTACTS" \
    --max-sidechain-min-dist "$MAX_SC_MIN_DIST" \
    --score-weight-contact-count "$SC_CONTACT_SCORE_W"
} 2>&1 | tee "$LOG"
