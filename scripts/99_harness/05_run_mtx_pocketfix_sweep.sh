#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_ROOT="${COMMON_ROOT:-$(git -C "$ROOT" worktree list --porcelain | awk '/^worktree / {print $2; exit}')}"
LOG="${ROOT}/logs/99_harness/mtx_pocketfix_sweep.log"
mkdir -p "$(dirname "$LOG")"

PYTHON_BIN="${PYTHON_BIN:-${COMMON_ROOT}/.venv/bin/python}"
LIGAND_PDB="${LIGAND_PDB:-${COMMON_ROOT}/inputs/01_cgmap/MTX.pdb}"
POLAR_SITES="${POLAR_SITES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_polar_sites.json}"
SITE_FRAMES="${SITE_FRAMES:-${COMMON_ROOT}/outputs/02_polar_sites/MTX_site_frames.json}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi

{
  echo "[run] repo: $ROOT"
  echo "[run] common_root: $COMMON_ROOT"
  echo "[run] python: $PYTHON_BIN"
  "$PYTHON_BIN" "${ROOT}/scripts/99_harness/05_mtx_pocketfix_sweep.py" \
    --repo-root "$ROOT" \
    --ligand-pdb "$LIGAND_PDB" \
    --polar-sites "$POLAR_SITES" \
    --site-frames "$SITE_FRAMES" \
    --vdxform-dir "$VDXFORM_DIR" \
    --top-per-site "${TOP_PER_SITE:-200}" \
    --solver "${SOLVER:-cp_sat}" \
    --time-limit-s "${TIME_LIMIT_S:-30}" \
    --python "$PYTHON_BIN" \
    "${@}"
} 2>&1 | tee "$LOG"
