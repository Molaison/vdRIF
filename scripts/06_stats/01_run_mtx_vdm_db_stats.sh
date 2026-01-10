#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

DB="${HOME}/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs"
OUT="${ROOT}/processed/06_stats/mtx_vdm_db_stats.json"
LOG="${ROOT}/logs/06_stats/mtx_vdm_db_stats.log"

mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"

{
  echo "[run] db: $DB"
  echo "[run] out: $OUT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/06_stats/01_vdm_db_stats.py" \
    --vdm-db "$DB" \
    --cgs coo conh2 ccn ph bb_cnh \
    -o "$OUT"
} 2>&1 | tee "$LOG"

echo "Wrote: $OUT"

