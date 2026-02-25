#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_ROOT="${COMMON_ROOT:-$(git -C "$ROOT" worktree list --porcelain | awk '/^worktree / {print $2; exit}')}"

TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
PYTHON_RUN="${PYTHON_RUN:-${COMMON_ROOT}/.venv/bin/python}"
PYTHON_PREP="${PYTHON_PREP:-$PYTHON_RUN}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform_full}"

LIGANDS="${LIGANDS:-PAB:Nc1ccc(cc1)C(=O)O,PBA:[NH3+]c1ccc(cc1)C(=O)[O-]}"
POLAR_BACKEND="${POLAR_BACKEND:-auto}"
SOLVER="${SOLVER:-cp_sat}"
TOP_PER_SITE="${TOP_PER_SITE:-120}"
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-40}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
TIME_LIMIT_S="${TIME_LIMIT_S:-120}"
CA_PREFILTER="${CA_PREFILTER:-14.0}"
ACCEPTOR_MODEL="${ACCEPTOR_MODEL:-plip}"
MAX_RUNS="${MAX_RUNS:-8}"

SCORE_W_COVERAGE_LIST="${SCORE_W_COVERAGE_LIST:-0.03,0.05}"
SCORE_W_CONTACT_LIST="${SCORE_W_CONTACT_LIST:-0.1,0.2}"
SCORE_W_SHELL_LIST="${SCORE_W_SHELL_LIST:-0.2}"
TARGET_RES_LIST="${TARGET_RES_LIST:-10,12}"
MIN_COVER_PER_POLAR_LIST="${MIN_COVER_PER_POLAR_LIST:-1}"

RUN_PLIP="${RUN_PLIP:-1}"
REQUIRE_PLIP_SUCCESS="${REQUIRE_PLIP_SUCCESS:-1}"
PLIP_BIN="${PLIP_BIN:-$(command -v plip 2>/dev/null || echo plip)}"

LOG="${ROOT}/logs/99_harness/non_mtx_transfer_benchmark_${TAG}.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] root: $ROOT"
  echo "[run] common_root: $COMMON_ROOT"
  echo "[run] tag: $TAG"
  echo "[run] python_run: $PYTHON_RUN"
  echo "[run] python_prep: $PYTHON_PREP"
  echo "[run] vdxform_dir: $VDXFORM_DIR"
  echo "[run] ligands: $LIGANDS"
  echo "[run] backend=$POLAR_BACKEND solver=$SOLVER top_per_site=$TOP_PER_SITE top_per_site_per_atom=$TOP_PER_SITE_PER_ATOM chunk_size=$CHUNK_SIZE time_limit_s=$TIME_LIMIT_S ca_prefilter=$CA_PREFILTER"
  echo "[run] score_w_coverage_list=$SCORE_W_COVERAGE_LIST score_w_contact_list=$SCORE_W_CONTACT_LIST score_w_shell_list=$SCORE_W_SHELL_LIST"
  echo "[run] target_res_list=$TARGET_RES_LIST min_cover_per_polar_list=$MIN_COVER_PER_POLAR_LIST max_runs=$MAX_RUNS"
  echo "[run] run_plip=$RUN_PLIP require_plip_success=$REQUIRE_PLIP_SUCCESS plip_bin=$PLIP_BIN"

  EXTRA_FLAGS=()
  if [[ "$RUN_PLIP" == "1" ]]; then
    EXTRA_FLAGS+=(--run-plip --plip-bin "$PLIP_BIN")
  fi
  if [[ "$REQUIRE_PLIP_SUCCESS" == "1" ]]; then
    EXTRA_FLAGS+=(--require-plip-success)
  fi

  "$PYTHON_RUN" "$ROOT/scripts/99_harness/06_non_mtx_transfer_benchmark.py" \
    --repo-root "$ROOT" \
    --tag "$TAG" \
    --python-run "$PYTHON_RUN" \
    --python-prep "$PYTHON_PREP" \
    --vdxform-dir "$VDXFORM_DIR" \
    --ligands "$LIGANDS" \
    --polar-backend "$POLAR_BACKEND" \
    --solver "$SOLVER" \
    --top-per-site "$TOP_PER_SITE" \
    --top-per-site-per-atom "$TOP_PER_SITE_PER_ATOM" \
    --chunk-size "$CHUNK_SIZE" \
    --time-limit-s "$TIME_LIMIT_S" \
    --ca-prefilter "$CA_PREFILTER" \
    --acceptor-model "$ACCEPTOR_MODEL" \
    --score-w-coverage-list "$SCORE_W_COVERAGE_LIST" \
    --score-w-contact-list "$SCORE_W_CONTACT_LIST" \
    --score-w-shell-list "$SCORE_W_SHELL_LIST" \
    --target-res-list "$TARGET_RES_LIST" \
    --min-cover-per-polar-list "$MIN_COVER_PER_POLAR_LIST" \
    --max-runs "$MAX_RUNS" \
    "${EXTRA_FLAGS[@]}"
} 2>&1 | tee "$LOG"
