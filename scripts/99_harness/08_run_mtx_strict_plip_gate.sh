#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TAG="${1:-$(date +%Y%m%d-%H%M%S)}"
OUTDIR="${ROOT}/processed/99_harness/mtx_strict_plip_gate_${TAG}"
LOG="${ROOT}/logs/99_harness/mtx_strict_plip_gate.log"
PLIP_BIN="${PLIP_BIN:-/xcfhome/zpzeng/miniconda/bin/plip}"

mkdir -p "$OUTDIR" "$(dirname "$LOG")"

{
  echo "[run] repo: $ROOT"
  echo "[run] tag: $TAG"
  echo "[run] outdir: $OUTDIR"
  echo "[run] plip_bin: $PLIP_BIN"

  # Hard cutover: always regenerate polar sites from current extraction logic.
  POLAR_SITES_BACKEND="${POLAR_SITES_BACKEND:-rdkit}" bash "${ROOT}/scripts/02_polar_sites/01_run_mtx_polar_sites.sh"
  bash "${ROOT}/scripts/04_candidates/01_run_mtx_candidates_debug.sh"
  bash "${ROOT}/scripts/05_solver/01_run_mtx_solver_debug.sh"

  uv run -p 3.11 python "${ROOT}/scripts/05_solver/06_validate_motif_plip.py" \
    --motif-pdb "${ROOT}/outputs/05_solver/MTX_motif_debug.pdb" \
    --polar-sites "${ROOT}/outputs/02_polar_sites/MTX_polar_sites.json" \
    --plip-bin "$PLIP_BIN" \
    -o "${OUTDIR}/MTX_motif_debug_plip_validation.json"

  python3 - <<'PY' "$ROOT" "$OUTDIR" "$TAG"
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
outdir = Path(sys.argv[2])
tag = str(sys.argv[3])

polar_path = root / "outputs/02_polar_sites/MTX_polar_sites.json"
plip_path = outdir / "MTX_motif_debug_plip_validation.json"
sat_path = root / "outputs/05_solver/MTX_motif_debug_validation.json"
clash_path = root / "outputs/05_solver/MTX_motif_debug_clash_validation.json"
motif_pdb = root / "outputs/05_solver/MTX_motif_debug.pdb"
cand_npz = root / "processed/04_candidates/MTX_candidates_debug.npz"
cand_meta = root / "processed/04_candidates/MTX_candidates_debug.json"

polar = json.loads(polar_path.read_text(encoding="utf-8"))
plip = json.loads(plip_path.read_text(encoding="utf-8"))
sat = json.loads(sat_path.read_text(encoding="utf-8"))
clash = json.loads(clash_path.read_text(encoding="utf-8"))

n5_roles = []
for s in polar.get("sites", []):
    if str(s.get("atom_name")) == "N5":
        n5_roles = list(s.get("roles") or [])
        break

report = {
    "tag": tag,
    "inputs": {
        "polar_sites": str(polar_path),
        "motif_pdb": str(motif_pdb),
        "candidates_npz": str(cand_npz),
        "candidates_meta": str(cand_meta),
    },
    "metrics": {
        "n_polar_atoms": int(plip.get("n_polar_atoms", 0)),
        "n_plip_satisfied": int(plip.get("n_satisfied", 0)),
        "plip_all_satisfied": bool(plip.get("all_satisfied")),
        "sat_all_satisfied": bool(sat.get("all_satisfied")),
        "clash_ok": bool(clash.get("ok")),
    },
    "n5": {
        "roles": n5_roles,
        "has_cation_role": bool("cation" in set(str(x) for x in n5_roles)),
    },
    "unsatisfied_polar_atoms": list(plip.get("unsatisfied_polar_atoms") or []),
}

(outdir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

if not report["metrics"]["plip_all_satisfied"]:
    raise SystemExit("Strict PLIP gate failed: motif is not fully PLIP-satisfied.")
if report["n5"]["has_cation_role"]:
    raise SystemExit("Strict PLIP gate failed: stale N5 cation target still present in MTX polar sites.")
if not report["metrics"]["clash_ok"]:
    raise SystemExit("Strict PLIP gate failed: motif has ligand clashes.")

print(f"[done] strict report: {outdir / 'report.json'}")
PY
} 2>&1 | tee "$LOG"

