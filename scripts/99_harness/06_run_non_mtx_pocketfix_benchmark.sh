#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_ROOT="${COMMON_ROOT:-$(git -C "$ROOT" worktree list --porcelain | awk '/^worktree / {print $2; exit}')}"
LOG="${ROOT}/logs/99_harness/non_mtx_pocketfix_benchmark.log"
mkdir -p "$(dirname "$LOG")"

PYTHON_RUN="${PYTHON_RUN:-${COMMON_ROOT}/.venv/bin/python}"
PYTHON_PREP="${PYTHON_PREP:-python3}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform}"

{
  echo "[run] repo: $ROOT"
  echo "[run] common_root: $COMMON_ROOT"
  echo "[run] python_run: $PYTHON_RUN"
  echo "[run] python_prep: $PYTHON_PREP"
  "$PYTHON_RUN" "${ROOT}/scripts/99_harness/06_non_mtx_pocketfix_benchmark.py" \
    --repo-root "$ROOT" \
    --python-run "$PYTHON_RUN" \
    --python-prep "$PYTHON_PREP" \
    --vdxform-dir "$VDXFORM_DIR" \
    --configs "${CONFIGS:-baseline_legacy,pocketfix_default,pocketfix_weighted,pocketfix_strict}" \
    --ligands "${LIGANDS:-NIC:NC(=O)c1ccncc1,BZA:O=C(O)c1ccccc1}" \
    --top-per-site "${TOP_PER_SITE:-200}" \
    --solver "${SOLVER:-cp_sat}" \
    --time-limit-s "${TIME_LIMIT_S:-30}" \
    --max-internal-overlap-rank "${MAX_INTERNAL_OVERLAP_RANK:-0.25}" \
    "${@}"
} 2>&1 | tee "$LOG"
