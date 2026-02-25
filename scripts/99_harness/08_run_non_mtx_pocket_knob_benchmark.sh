#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT" rev-parse --path-format=absolute --git-common-dir)"
COMMON_ROOT="$(dirname "$COMMON_GIT_DIR")"

TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
PYTHON_BIN="${PYTHON_BIN:-${COMMON_ROOT}/.venv/bin/python}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform_full}"
PREPARED_ROOT="${PREPARED_ROOT:-${ROOT}/processed/99_harness/non_mtx_cation_validate_20260226-cation-nonmtx-v2/ligands}"
LIGANDS="${LIGANDS:-GDM,BZA,BTM}"
SOLVER="${SOLVER:-greedy}"
TOP_PER_SITE="${TOP_PER_SITE:-200}"
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-50}"
TIME_LIMIT_S="${TIME_LIMIT_S:-120}"
MIN_RES="${MIN_RES:-3}"
MAX_RES="${MAX_RES:-12}"
CA_PREFILTER="${CA_PREFILTER:-12.0}"
ANGLE_LIST="${ANGLE_LIST:-60}"
WEIGHT_LIST="${WEIGHT_LIST:-0.0,0.2}"
MIN_SC_LIST="${MIN_SC_LIST:-0,1}"
MAX_RUNS="${MAX_RUNS:-0}"

LOG="${ROOT}/logs/99_harness/non_mtx_pocket_knob_benchmark_${TAG}.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  echo "[run] tag: $TAG"
  echo "[run] python_bin: $PYTHON_BIN"
  echo "[run] prepared_root: $PREPARED_ROOT"
  echo "[run] vdxform_dir: $VDXFORM_DIR"
  echo "[run] ligands: $LIGANDS"
  echo "[run] solver: $SOLVER"
  echo "[run] knobs: angle=$ANGLE_LIST weight=$WEIGHT_LIST min_sc=$MIN_SC_LIST"
  echo

  "$PYTHON_BIN" "$ROOT/scripts/99_harness/08_non_mtx_pocket_knob_benchmark.py" \
    --repo-root "$ROOT" \
    --tag "$TAG" \
    --python-bin "$PYTHON_BIN" \
    --prepared-root "$PREPARED_ROOT" \
    --ligands "$LIGANDS" \
    --vdxform-dir "$VDXFORM_DIR" \
    --solver "$SOLVER" \
    --top-per-site "$TOP_PER_SITE" \
    --top-per-site-per-atom "$TOP_PER_SITE_PER_ATOM" \
    --time-limit-s "$TIME_LIMIT_S" \
    --min-res "$MIN_RES" \
    --max-res "$MAX_RES" \
    --ca-prefilter "$CA_PREFILTER" \
    --angle-list "$ANGLE_LIST" \
    --weight-list "$WEIGHT_LIST" \
    --min-sc-list "$MIN_SC_LIST" \
    --max-runs "$MAX_RUNS"
} 2>&1 | tee "$LOG"

