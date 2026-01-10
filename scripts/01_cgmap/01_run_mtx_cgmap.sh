#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LIGAND_IN="${ROOT}/inputs/01_cgmap/MTX.pdb"
OUT_JSON="${ROOT}/outputs/01_cgmap/MTX_cg_atommap.json"
LOG="${ROOT}/logs/01_cgmap/MTX_cg_atommap.log"

VDM_DB="${HOME}/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs"
CGMAP="${HOME}/practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py"

mkdir -p "$(dirname "${OUT_JSON}")" "$(dirname "${LOG}")"

{
  date
  echo "ligand=${LIGAND_IN}"
  echo "vdm_db=${VDM_DB}"
  echo "script=${CGMAP}"
  echo
  uv run -p 3.11 python "${CGMAP}" "${LIGAND_IN}" \
    --format json \
    --vdm-database "${VDM_DB}" \
    -o "${OUT_JSON}"
} 2>&1 | tee "${LOG}"

echo "Wrote: ${OUT_JSON}"
