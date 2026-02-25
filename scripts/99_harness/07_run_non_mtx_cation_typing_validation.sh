#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="${ROOT}/logs/99_harness/non_mtx_cation_typing_validation.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/99_harness/07_non_mtx_cation_typing_validation.py" \
    --repo-root "$ROOT" \
    --tag "$(date +%Y%m%d-%H%M%S)"
} 2>&1 | tee "$LOG"

