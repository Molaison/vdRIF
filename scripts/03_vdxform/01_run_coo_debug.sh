#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

DB="${HOME}/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs"
OUT="${ROOT}/processed/03_vdxform"
LOG="${ROOT}/logs/03_vdxform/coo_debug.log"

mkdir -p "$(dirname "$LOG")"

{
  echo "[run] db: $DB"
  echo "[run] out: $OUT"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/03_vdxform/01_convert_parquet_to_vdxform_parts.py" \
    --vdm-db "$DB" \
    --cg coo \
    --outdir "$OUT" \
    --frame-defs "${ROOT}/configs/cg_frame_defs.json" \
    --workers 1 \
    --max-instances-per-file 500 \
    --resume

  uv run -p 3.11 python "${ROOT}/scripts/03_vdxform/02_pack_vdxform_npz.py" \
    --parts-dir "${OUT}/coo/parts" \
    -o "${OUT}/coo/vdxform_coo.npz" \
    --resume
} 2>&1 | tee "$LOG"
