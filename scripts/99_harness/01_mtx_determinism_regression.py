#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _ensure_file(path: Path, hint_cmd: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}\nHint: {hint_cmd}")


@dataclass(frozen=True)
class RunOut:
    candidates_npz: Path
    candidates_json: Path
    motif_json: Path
    motif_pdb: Path
    motif_validation_json: Path
    motif_clash_validation_json: Path
    motif_internal_clash_validation_json: Path


def main() -> None:
    ap = argparse.ArgumentParser(description="Determinism regression for MTX: candidates + solver must be byte-identical.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--time-limit-s", type=float, default=30.0)
    ap.add_argument(
        "--solver",
        type=str,
        default="cp_sat",
        choices=["cp_sat", "greedy"],
        help="Which solver to use for motif selection (default cp_sat).",
    )
    ap.add_argument(
        "--vdxform-dir",
        type=Path,
        default=Path("processed/03_vdxform"),
        help="Directory containing per-cg vdxform_<cg>.npz files.",
    )
    ap.add_argument(
        "--allow-backbone-only",
        action="store_true",
        help=(
            "Allow motif PDBs to be written as backbone-only stubs (N/CA/C/CB) when candidates lack "
            "`center_atom_xyz_stub_f32` or candidates meta lacks `atom_order`."
        ),
    )
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"mtx_det_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand = root / "inputs/01_cgmap/MTX.pdb"
    polar = root / "outputs/02_polar_sites/MTX_polar_sites.json"
    frames = root / "outputs/02_polar_sites/MTX_site_frames.json"

    vdx = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()
    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        _ensure_file(
            vdx / cg / f"vdxform_{cg}.npz",
            hint_cmd="bash scripts/03_vdxform/03_run_mtx_needed_cgs_debug.sh (debug) or 02_run_mtx_needed_cgs.sh (full)",
        )

    _ensure_file(ligand, hint_cmd="Place MTX at inputs/01_cgmap/MTX.pdb")
    _ensure_file(polar, hint_cmd="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")
    _ensure_file(frames, hint_cmd="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")

    candidate_py = root / "scripts/04_candidates/01_generate_candidates.py"
    solver_py = root / "scripts/05_solver/01_solve_motif.py"
    validate_sat_py = root / "scripts/05_solver/03_validate_motif_polar_satisfaction.py"
    validate_clash_py = root / "scripts/05_solver/04_validate_motif_clashes.py"
    validate_internal_clash_py = root / "scripts/05_solver/05_validate_motif_internal_clashes.py"

    def run_once(run_name: str) -> RunOut:
        run_dir = outdir / run_name
        run_dir.mkdir(parents=True, exist_ok=True)
        cand_prefix = run_dir / "MTX_candidates"
        motif_json = run_dir / "MTX_motif.json"
        motif_pdb = run_dir / "MTX_motif.pdb"
        motif_val_json = run_dir / "MTX_motif_validation.json"
        motif_clash_json = run_dir / "MTX_motif_clash_validation.json"
        motif_internal_clash_json = run_dir / "MTX_motif_internal_clash_validation.json"

        t0 = time.time()
        _run(
            [
                "uv",
                "run",
                "-p",
                "3.11",
                "python",
                str(candidate_py),
                "--ligand-pdb",
                str(ligand),
                "--polar-sites",
                str(polar),
                "--site-frames",
                str(frames),
                "--vdxform-dir",
                str(vdx),
                "--out-prefix",
                str(cand_prefix),
                "--chunk-size",
                "5000",
                "--top-per-site",
                str(args.top_per_site),
                "--clash-tol",
                "0.5",
                "--exclude-aa3",
                "PRO,CYS",
                "--require-sidechain-facing",
                "--min-sidechain-facing-dot",
                "0.2",
                "--min-sidechain-centroid-dot",
                "0.0",
                "--require-full-coverage",
            ]
        )

        # Guardrail: if candidates lack full-atom coordinates, the solver will fall back to
        # backbone-only stub residues (N/CA/C/CB), which can look like every residue is "C-C(-N)-C".
        cand_meta = json.loads(cand_prefix.with_suffix(".json").read_text(encoding="utf-8"))
        has_atom_order = bool(cand_meta.get("atom_order"))
        z = np.load(cand_prefix.with_suffix(".npz"))
        has_center = "center_atom_xyz_stub_f32" in z.files
        if not args.allow_backbone_only and (not has_atom_order or not has_center):
            raise RuntimeError(
                "Candidates are missing full-atom residue coordinates required for full-atom motif PDB output. "
                f"has_atom_order={has_atom_order} has_center_atom_xyz_stub_f32={has_center}. "
                "Re-run vdXform conversion + candidate generation with the current pipeline."
            )

        t1 = time.time()
        _run(
            [
                "uv",
                "run",
                "-p",
                "3.11",
                "python",
                str(solver_py),
                "--candidates-npz",
                str(cand_prefix.with_suffix(".npz")),
                "--candidates-meta",
                str(cand_prefix.with_suffix(".json")),
                "--ligand-pdb",
                str(ligand),
                "--out-json",
                str(motif_json),
                "--out-pdb",
                str(motif_pdb),
                "--solver",
                str(args.solver),
                "--min-res",
                "8",
                "--max-res",
                "15",
                "--time-limit-s",
                str(args.time_limit_s),
                "--num-workers",
                "1",
                "--grid-size",
                "4.0",
                "--ca-prefilter",
                "12.0",
                "--clash-tol",
                "0.5",
            ]
        )
        t2 = time.time()

        _run(
            [
                "uv",
                "run",
                "-p",
                "3.11",
                "python",
                str(validate_sat_py),
                "--polar-sites",
                str(polar),
                "--motif-pdb",
                str(motif_pdb),
                "-o",
                str(motif_val_json),
            ]
        )
        _run(
            [
                "uv",
                "run",
                "-p",
                "3.11",
                "python",
                str(validate_clash_py),
                "--motif-pdb",
                str(motif_pdb),
                "--ligand-resname",
                "MTX",
                "-o",
                str(motif_clash_json),
            ]
        )
        _run(
            [
                "uv",
                "run",
                "-p",
                "3.11",
                "python",
                str(validate_internal_clash_py),
                "--motif-pdb",
                str(motif_pdb),
                "--ligand-resname",
                "MTX",
                "-o",
                str(motif_internal_clash_json),
            ]
        )

        val = json.loads(motif_val_json.read_text(encoding="utf-8"))
        if not bool(val.get("all_satisfied")):
            raise RuntimeError(f"Polar satisfaction validation FAILED for run {run_name}: {motif_val_json}")

        clash = json.loads(motif_clash_json.read_text(encoding="utf-8"))
        if not bool(clash.get("ok")):
            raise RuntimeError(f"Clash validation FAILED for run {run_name}: {motif_clash_json}")

        internal = json.loads(motif_internal_clash_json.read_text(encoding="utf-8"))
        if not bool(internal.get("ok")):
            raise RuntimeError(f"Internal clash validation FAILED for run {run_name}: {motif_internal_clash_json}")

        return RunOut(
            candidates_npz=cand_prefix.with_suffix(".npz"),
            candidates_json=cand_prefix.with_suffix(".json"),
            motif_json=motif_json,
            motif_pdb=motif_pdb,
            motif_validation_json=motif_val_json,
            motif_clash_validation_json=motif_clash_json,
            motif_internal_clash_validation_json=motif_internal_clash_json,
        ), {"candidates_s": t1 - t0, "solve_s": t2 - t1, "total_s": t2 - t0}

    A, statsA = run_once("A")
    B, statsB = run_once("B")

    files = [
        ("candidates.npz", A.candidates_npz, B.candidates_npz),
        ("candidates.json", A.candidates_json, B.candidates_json),
        ("motif.json", A.motif_json, B.motif_json),
        ("motif.pdb", A.motif_pdb, B.motif_pdb),
    ]
    checks = []
    ok = True
    for label, pa, pb in files:
        ha = _sha256(pa)
        hb = _sha256(pb)
        same = ha == hb
        ok = ok and same
        checks.append({"label": label, "A": str(pa), "B": str(pb), "sha256_A": ha, "sha256_B": hb, "same": same})

    report: dict[str, Any] = {
        "tag": args.tag,
        "params": {"top_per_site": args.top_per_site, "time_limit_s": args.time_limit_s},
        "stats": {"A": statsA, "B": statsB},
        "checks": checks,
        "deterministic": bool(ok),
        "outdir": str(outdir),
    }
    (outdir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if not ok:
        raise SystemExit(f"Determinism check FAILED; see {outdir/'report.json'}")


if __name__ == "__main__":
    main()
