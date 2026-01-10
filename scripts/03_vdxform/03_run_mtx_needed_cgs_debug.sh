#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

DB="${HOME}/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs"
OUT="${ROOT}/processed/03_vdxform"
LOG="${ROOT}/logs/03_vdxform/mtx_needed_cgs_debug.log"

CGS=(coo conh2 ccn ph bb_cnh)

mkdir -p "$(dirname "$LOG")"

{
  echo "[run] db: $DB"
  echo "[run] out: $OUT"
  echo "[run] cgs: ${CGS[*]}"
  uv sync -p 3.11 --extra rdkit

  for cg in "${CGS[@]}"; do
    echo "[convert-debug] $cg"
    uv run -p 3.11 python "${ROOT}/scripts/03_vdxform/01_convert_parquet_to_vdxform_parts.py" \
      --vdm-db "$DB" \
      --cg "$cg" \
      --outdir "$OUT" \
      --frame-defs "${ROOT}/configs/cg_frame_defs.json" \
      --workers 1 \
      --max-instances-per-file 500 \
      --resume

    uv run -p 3.11 python "${ROOT}/scripts/03_vdxform/02_pack_vdxform_npz.py" \
      --parts-dir "${OUT}/${cg}/parts" \
      -o "${OUT}/${cg}/vdxform_${cg}.npz" \
      --resume
  done
} 2>&1 | tee "$LOG"

