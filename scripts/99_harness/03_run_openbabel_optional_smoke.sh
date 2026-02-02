#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="${ROOT}/logs/99_harness/openbabel_optional_smoke.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/99_harness/03_openbabel_optional_smoke.py" --repo-root "$ROOT"
} 2>&1 | tee "$LOG"

