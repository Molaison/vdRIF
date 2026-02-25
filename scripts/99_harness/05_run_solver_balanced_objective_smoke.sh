#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG="${ROOT}/logs/99_harness/solver_balanced_objective_smoke.log"
mkdir -p "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  uv run --no-project -p 3.11 --with numpy --with ortools python \
    "${ROOT}/scripts/99_harness/05_solver_balanced_objective_smoke.py" \
    --repo-root "$ROOT"
} 2>&1 | tee "$LOG"
