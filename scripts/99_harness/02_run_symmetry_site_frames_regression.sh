#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="${ROOT}/logs/99_harness/symmetry_site_frames_regression.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/99_harness/02_symmetry_site_frames_regression.py" --repo-root "$ROOT"
} 2>&1 | tee "$LOG"

