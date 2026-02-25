#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT" rev-parse --path-format=absolute --git-common-dir)"
COMMON_ROOT="$(dirname "$COMMON_GIT_DIR")"

TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
PY311="${PY311:-${COMMON_ROOT}/.venv/bin/python}"
VDXFORM_DIR="${VDXFORM_DIR:-${COMMON_ROOT}/processed/03_vdxform_full}"
CGMAP_SCRIPT="${CGMAP_SCRIPT:-${HOME}/practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py}"
VDM_DB="${VDM_DB:-${HOME}/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs}"

LOG="${ROOT}/logs/99_harness/non_mtx_cation_validate_${TAG}.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  echo "[run] tag: $TAG"
  echo "[run] py311: $PY311"
  echo "[run] vdxform_dir: $VDXFORM_DIR"
  echo "[run] cgmap_script: $CGMAP_SCRIPT"
  echo "[run] vdm_db: $VDM_DB"
  echo

  python3 "$ROOT/scripts/99_harness/07_validate_cation_typing_non_mtx.py" \
    --repo-root "$ROOT" \
    --tag "$TAG" \
    --python311 "$PY311" \
    --vdxform-dir "$VDXFORM_DIR" \
    --cgmap-script "$CGMAP_SCRIPT" \
    --vdm-db "$VDM_DB"
} 2>&1 | tee "$LOG"

