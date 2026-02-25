#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT" rev-parse --path-format=absolute --git-common-dir)"
COMMON_ROOT="$(dirname "$COMMON_GIT_DIR")"

TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform_full}"
LIGAND_PDB="${LIGAND_PDB:-${COMMON_ROOT}/inputs/01_cgmap/MTX.pdb}"
POLAR_SITES="${POLAR_SITES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_polar_sites.json}"
SITE_FRAMES="${SITE_FRAMES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_site_frames.json}"
PYTHON_BIN="${PYTHON_BIN:-${COMMON_ROOT}/.venv/bin/python}"
SOLVER="${SOLVER:-greedy}"
TOP_PER_SITE="${TOP_PER_SITE:-200}"
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-50}"
ANGLE_LIST="${ANGLE_LIST:-60,75}"
WEIGHT_LIST="${WEIGHT_LIST:-0.0,0.2}"
MIN_SC_LIST="${MIN_SC_LIST:-0,1}"
TIME_LIMIT_S="${TIME_LIMIT_S:-120}"
MAX_RUNS="${MAX_RUNS:-0}"

LOG="${ROOT}/logs/99_harness/mtx_pocket_knob_sweep_${TAG}.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  echo "[run] tag: $TAG"
  echo "[run] vdxform: $VDXFORM_DIR"
  echo "[run] ligand_pdb: $LIGAND_PDB"
  echo "[run] polar_sites: $POLAR_SITES"
  echo "[run] site_frames: $SITE_FRAMES"
  echo "[run] python_bin: $PYTHON_BIN"
  echo "[run] solver: $SOLVER"
  echo "[run] top_per_site: $TOP_PER_SITE"
  echo "[run] top_per_site_per_atom: $TOP_PER_SITE_PER_ATOM"
  echo "[run] angle_list: $ANGLE_LIST"
  echo "[run] weight_list: $WEIGHT_LIST"
  echo "[run] min_sc_list: $MIN_SC_LIST"
  echo "[run] time_limit_s: $TIME_LIMIT_S"
  echo "[run] max_runs: $MAX_RUNS"

  "$PYTHON_BIN" "${ROOT}/scripts/99_harness/05_mtx_pocket_knob_sweep.py" \
    --repo-root "$ROOT" \
    --tag "$TAG" \
    --python-bin "$PYTHON_BIN" \
    --ligand-pdb "$LIGAND_PDB" \
    --polar-sites "$POLAR_SITES" \
    --site-frames "$SITE_FRAMES" \
    --vdxform-dir "$VDXFORM_DIR" \
    --solver "$SOLVER" \
    --top-per-site "$TOP_PER_SITE" \
    --top-per-site-per-atom "$TOP_PER_SITE_PER_ATOM" \
    --time-limit-s "$TIME_LIMIT_S" \
    --min-lig-donor-angle-list "$ANGLE_LIST" \
    --pocket-contact-weight-list "$WEIGHT_LIST" \
    --min-pocket-sidechain-contacts-list "$MIN_SC_LIST" \
    --max-runs "$MAX_RUNS"
} 2>&1 | tee "$LOG"
