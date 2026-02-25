#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_ROOT="${COMMON_ROOT:-$(git -C "$ROOT" worktree list --porcelain | awk '/^worktree / {print $2; exit}')}"

TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
PYTHON_BIN="${PYTHON_BIN:-${COMMON_ROOT}/.venv/bin/python}"
LIGAND_PDB="${LIGAND_PDB:-${COMMON_ROOT}/inputs/01_cgmap/MTX.pdb}"
POLAR_SITES="${POLAR_SITES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_polar_sites.json}"
SITE_FRAMES="${SITE_FRAMES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_site_frames.json}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform_full}"

SOLVER="${SOLVER:-cp_sat}"
TOP_PER_SITE="${TOP_PER_SITE:-120}"
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-40}"
TIME_LIMIT_S="${TIME_LIMIT_S:-120}"
CA_PREFILTER="${CA_PREFILTER:-14.0}"
MAX_RUNS="${MAX_RUNS:-4}"

SCORE_W_COVERAGE_LIST="${SCORE_W_COVERAGE_LIST:-0.03,0.05}"
SCORE_W_CONTACT_LIST="${SCORE_W_CONTACT_LIST:-0.1,0.2}"
SCORE_W_SHELL_LIST="${SCORE_W_SHELL_LIST:-0.2}"
TARGET_RES_LIST="${TARGET_RES_LIST:-10,12}"
MIN_COVER_PER_POLAR_LIST="${MIN_COVER_PER_POLAR_LIST:-1}"

LOG="${ROOT}/logs/99_harness/mtx_pocket_quality_sweep_${TAG}.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] root: $ROOT"
  echo "[run] common_root: $COMMON_ROOT"
  echo "[run] tag: $TAG"
  echo "[run] python_bin: $PYTHON_BIN"
  echo "[run] ligand: $LIGAND_PDB"
  echo "[run] polar_sites: $POLAR_SITES"
  echo "[run] site_frames: $SITE_FRAMES"
  echo "[run] vdxform_dir: $VDXFORM_DIR"
  echo "[run] solver=$SOLVER top_per_site=$TOP_PER_SITE top_per_site_per_atom=$TOP_PER_SITE_PER_ATOM time_limit_s=$TIME_LIMIT_S ca_prefilter=$CA_PREFILTER max_runs=$MAX_RUNS"
  echo "[run] score_w_coverage_list=$SCORE_W_COVERAGE_LIST score_w_contact_list=$SCORE_W_CONTACT_LIST score_w_shell_list=$SCORE_W_SHELL_LIST"
  echo "[run] target_res_list=$TARGET_RES_LIST min_cover_per_polar_list=$MIN_COVER_PER_POLAR_LIST"
  "$PYTHON_BIN" "$ROOT/scripts/99_harness/05_mtx_pocket_quality_sweep.py" \
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
    --ca-prefilter "$CA_PREFILTER" \
    --score-w-coverage-list "$SCORE_W_COVERAGE_LIST" \
    --score-w-contact-list "$SCORE_W_CONTACT_LIST" \
    --score-w-shell-list "$SCORE_W_SHELL_LIST" \
    --target-res-list "$TARGET_RES_LIST" \
    --min-cover-per-polar-list "$MIN_COVER_PER_POLAR_LIST" \
    --run-plip \
    --max-runs "$MAX_RUNS"
} 2>&1 | tee "$LOG"
