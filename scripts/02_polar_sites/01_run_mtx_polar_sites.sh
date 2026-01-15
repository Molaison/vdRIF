#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIG="${ROOT}/inputs/01_cgmap/MTX.pdb"
OUT_SITES="${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json"
OUT_COV="${ROOT}/outputs/02_polar_sites/MTX_polar_coverage_vs_cgmap.json"
OUT_FRAMES="${ROOT}/outputs/02_polar_sites/MTX_site_frames.json"
LOG="${ROOT}/logs/02_polar_sites/MTX_polar_sites.log"

mkdir -p "$(dirname "$OUT_SITES")" "$(dirname "$OUT_COV")" "$(dirname "$LOG")"

{
  echo "[run] ligand: $LIG"
  echo "[run] out_sites: $OUT_SITES"
  echo "[run] out_coverage: $OUT_COV"
  echo "[run] out_frames: $OUT_FRAMES"
  uv sync -p 3.11 --extra rdkit
  uv run -p 3.11 python "${ROOT}/scripts/02_polar_sites/01_extract_polar_sites.py" "$LIG" -o "$OUT_SITES" --backend openbabel
  uv run -p 3.11 python "${ROOT}/scripts/02_polar_sites/02_check_polar_coverage_vs_cgmap.py" \
    --polar-sites "$OUT_SITES" \
    --cg-atommap "${ROOT}/outputs/01_cgmap/MTX_cg_atommap.json" \
    -o "$OUT_COV"
  uv run -p 3.11 python "${ROOT}/scripts/02_polar_sites/03_build_ligand_site_frames.py" \
    "$LIG" \
    --cg-atommap "${ROOT}/outputs/01_cgmap/MTX_cg_atommap.json" \
    --polar-sites "$OUT_SITES" \
    --cg-frame-defs "${ROOT}/configs/cg_frame_defs.json" \
    --add-bb-cnh-for-uncovered-donors \
    --add-ccn-for-cations \
    -o "$OUT_FRAMES"
} 2>&1 | tee "$LOG"
