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
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-200}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
CLASH_TOL="${CLASH_TOL:-0.5}"
POCKET_CONTACT_CUTOFF="${POCKET_CONTACT_CUTOFF:-4.5}"
MIN_POCKET_CONTACTS="${MIN_POCKET_CONTACTS:-0}"
POCKET_CONTACT_WEIGHT="${POCKET_CONTACT_WEIGHT:-0.2}"
MIN_LIG_DONOR_ANGLE_DEG="${MIN_LIG_DONOR_ANGLE_DEG:-60}"

{
  echo "[run] ligand: $LIG"
  echo "[run] polar: $POLAR"
  echo "[run] frames: $FRAMES"
  echo "[run] vdx: $VDX"
  echo "[run] out_prefix: $OUT_PREFIX"
  echo "[run] TOP_PER_SITE=$TOP_PER_SITE TOP_PER_SITE_PER_ATOM=$TOP_PER_SITE_PER_ATOM CHUNK_SIZE=$CHUNK_SIZE CLASH_TOL=$CLASH_TOL"
  echo "[run] POCKET_CONTACT_CUTOFF=$POCKET_CONTACT_CUTOFF MIN_POCKET_CONTACTS=$MIN_POCKET_CONTACTS POCKET_CONTACT_WEIGHT=$POCKET_CONTACT_WEIGHT"
  echo "[run] MIN_LIG_DONOR_ANGLE_DEG=$MIN_LIG_DONOR_ANGLE_DEG"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/04_candidates/01_generate_candidates.py" \
    --ligand-pdb "$LIG" \
    --polar-sites "$POLAR" \
    --site-frames "$FRAMES" \
    --vdxform-dir "$VDX" \
    --out-prefix "$OUT_PREFIX" \
    --chunk-size "$CHUNK_SIZE" \
    --top-per-site "$TOP_PER_SITE" \
    --top-per-site-per-atom "$TOP_PER_SITE_PER_ATOM" \
    --allow-backbone-hbonds \
    --acceptor-model plip \
    --min-lig-donor-angle-deg "$MIN_LIG_DONOR_ANGLE_DEG" \
    --pocket-contact-cutoff "$POCKET_CONTACT_CUTOFF" \
    --min-pocket-sidechain-contacts "$MIN_POCKET_CONTACTS" \
    --pocket-contact-weight "$POCKET_CONTACT_WEIGHT" \
    --clash-tol "$CLASH_TOL" \
    --require-full-coverage
} 2>&1 | tee "$LOG"
