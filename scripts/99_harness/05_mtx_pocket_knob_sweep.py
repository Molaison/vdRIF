#!/usr/bin/env python
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd))


def _ensure_file(path: Path, hint_cmd: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}\nHint: {hint_cmd}")


def _parse_float_list(csv: str) -> list[float]:
    vals = [float(x.strip()) for x in csv.split(",") if x.strip()]
    if not vals:
        raise ValueError(f"Expected non-empty float list, got: {csv!r}")
    return vals


def _parse_int_list(csv: str) -> list[int]:
    vals = [int(x.strip()) for x in csv.split(",") if x.strip()]
    if not vals:
        raise ValueError(f"Expected non-empty int list, got: {csv!r}")
    return vals


def _safe_num(v: float | int) -> str:
    if isinstance(v, int):
        return str(v)
    s = f"{float(v):g}"
    return s.replace("-", "m").replace(".", "p")


def _read_json_if_exists(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _py_cmd(python_bin: str, script: Path, *args: str) -> list[str]:
    return [python_bin, str(script), *args]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Sweep MTX pocket-quality knobs on vdM candidate generation + motif solve, "
            "and summarize feasibility/quality/runtime."
        )
    )
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--python-bin", type=Path, default=None, help="Python interpreter used to run pipeline scripts.")
    ap.add_argument("--ligand-pdb", type=Path, default=Path("inputs/01_cgmap/MTX.pdb"))
    ap.add_argument("--polar-sites", type=Path, default=Path("outputs/02_polar_sites/MTX_polar_sites.json"))
    ap.add_argument("--site-frames", type=Path, default=Path("outputs/02_polar_sites/MTX_site_frames.json"))
    ap.add_argument(
        "--vdxform-dir",
        type=Path,
        default=Path("processed/03_vdxform_full"),
        help="Directory containing per-cg vdxform_<cg>.npz files.",
    )
    ap.add_argument("--solver", type=str, default="greedy", choices=["cp_sat", "greedy"])
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--top-per-site-per-atom", type=int, default=50)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument("--allow-backbone-hbonds", action="store_true", default=True)
    ap.add_argument("--acceptor-model", type=str, default="plip", choices=["legacy", "plip"])
    ap.add_argument("--pocket-contact-cutoff", type=float, default=4.5)
    ap.add_argument(
        "--min-lig-donor-angle-list",
        type=str,
        default="60,75",
        help="CSV float list for --min-lig-donor-angle-deg.",
    )
    ap.add_argument(
        "--pocket-contact-weight-list",
        type=str,
        default="0.0,0.2",
        help="CSV float list for --pocket-contact-weight.",
    )
    ap.add_argument(
        "--min-pocket-sidechain-contacts-list",
        type=str,
        default="0,1",
        help="CSV int list for --min-pocket-sidechain-contacts.",
    )
    ap.add_argument("--max-runs", type=int, default=0, help="If >0, only run the first N combinations.")
    ap.add_argument("--stop-on-error", action="store_true")
    args = ap.parse_args()

    root = args.repo_root.resolve()
    if args.python_bin is not None:
        python_bin_path = args.python_bin if args.python_bin.is_absolute() else (root / args.python_bin)
        python_bin = str(python_bin_path)
    else:
        default_py = root / ".venv" / "bin" / "python"
        python_bin = str(default_py if default_py.exists() else Path("python3"))
    outdir = root / "processed/99_harness" / f"mtx_pocket_knob_sweep_{args.tag}"
    runs_dir = outdir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    ligand = (root / args.ligand_pdb).resolve() if not args.ligand_pdb.is_absolute() else args.ligand_pdb.resolve()
    polar = (root / args.polar_sites).resolve() if not args.polar_sites.is_absolute() else args.polar_sites.resolve()
    frames = (root / args.site_frames).resolve() if not args.site_frames.is_absolute() else args.site_frames.resolve()

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

    min_lig_donor_angles = _parse_float_list(args.min_lig_donor_angle_list)
    pocket_contact_weights = _parse_float_list(args.pocket_contact_weight_list)
    min_pocket_sidechain_contacts = _parse_int_list(args.min_pocket_sidechain_contacts_list)

    combos = list(
        itertools.product(
            min_lig_donor_angles,
            pocket_contact_weights,
            min_pocket_sidechain_contacts,
        )
    )
    if int(args.max_runs) > 0:
        combos = combos[: int(args.max_runs)]
    if not combos:
        raise ValueError("No combinations to run.")

    runs: list[dict[str, Any]] = []
    for i, (min_lig_ang, pocket_w, min_sc_contacts) in enumerate(combos, start=1):
        run_id = (
            f"run_{i:02d}_ang{_safe_num(min_lig_ang)}_w{_safe_num(pocket_w)}_"
            f"minsc{_safe_num(min_sc_contacts)}"
        )
        run_dir = runs_dir / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        cand_prefix = run_dir / "MTX_candidates"
        motif_json = run_dir / "MTX_motif.json"
        motif_pdb = run_dir / "MTX_motif.pdb"
        motif_val_json = run_dir / "MTX_motif_validation.json"
        motif_clash_json = run_dir / "MTX_motif_clash_validation.json"
        motif_internal_clash_json = run_dir / "MTX_motif_internal_clash_validation.json"

        run_record: dict[str, Any] = {
            "run_id": run_id,
            "params": {
                "min_lig_donor_angle_deg": float(min_lig_ang),
                "pocket_contact_weight": float(pocket_w),
                "min_pocket_sidechain_contacts": int(min_sc_contacts),
                "acceptor_model": str(args.acceptor_model),
                "allow_backbone_hbonds": bool(args.allow_backbone_hbonds),
                "pocket_contact_cutoff": float(args.pocket_contact_cutoff),
                "top_per_site": int(args.top_per_site),
                "top_per_site_per_atom": int(args.top_per_site_per_atom),
                "chunk_size": int(args.chunk_size),
                "solver": str(args.solver),
                "time_limit_s": float(args.time_limit_s),
            },
            "paths": {
                "run_dir": str(run_dir),
                "candidates_npz": str(cand_prefix.with_suffix(".npz")),
                "candidates_json": str(cand_prefix.with_suffix(".json")),
                "motif_json": str(motif_json),
                "motif_pdb": str(motif_pdb),
                "motif_validation_json": str(motif_val_json),
                "motif_clash_validation_json": str(motif_clash_json),
                "motif_internal_clash_validation_json": str(motif_internal_clash_json),
            },
            "status": "ok",
            "timing_s": {},
            "error": None,
        }

        t0 = time.time()
        step = "candidates"
        try:
            cand_cmd = _py_cmd(
                python_bin,
                candidate_py,
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
                str(int(args.chunk_size)),
                "--top-per-site",
                str(int(args.top_per_site)),
                "--top-per-site-per-atom",
                str(int(args.top_per_site_per_atom)),
                "--acceptor-model",
                str(args.acceptor_model),
                "--min-lig-donor-angle-deg",
                str(float(min_lig_ang)),
                "--pocket-contact-cutoff",
                str(float(args.pocket_contact_cutoff)),
                "--min-pocket-sidechain-contacts",
                str(int(min_sc_contacts)),
                "--pocket-contact-weight",
                str(float(pocket_w)),
                "--clash-tol",
                str(float(args.clash_tol)),
                "--exclude-aa3",
                "PRO,CYS",
                "--require-sidechain-facing",
                "--min-sidechain-facing-dot",
                "0.2",
                "--min-sidechain-centroid-dot",
                "0.0",
                "--require-full-coverage",
            )
            if bool(args.allow_backbone_hbonds):
                cand_cmd.append("--allow-backbone-hbonds")
            _run(cand_cmd, cwd=root)
            t1 = time.time()

            step = "solver"
            _run(
                _py_cmd(
                    python_bin,
                    solver_py,
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
                    str(float(args.time_limit_s)),
                    "--num-workers",
                    "1",
                    "--grid-size",
                    "4.0",
                    "--ca-prefilter",
                    "8.0",
                    "--clash-tol",
                    str(float(args.clash_tol)),
                ),
                cwd=root,
            )
            t2 = time.time()

            step = "validate_sat"
            _run(
                _py_cmd(
                    python_bin,
                    validate_sat_py,
                    "--polar-sites",
                    str(polar),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--acceptor-model",
                    str(args.acceptor_model),
                    "-o",
                    str(motif_val_json),
                ),
                cwd=root,
            )
            step = "validate_ligand_clash"
            _run(
                _py_cmd(
                    python_bin,
                    validate_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    "MTX",
                    "-o",
                    str(motif_clash_json),
                ),
                cwd=root,
            )
            step = "validate_internal_clash"
            _run(
                _py_cmd(
                    python_bin,
                    validate_internal_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    "MTX",
                    "-o",
                    str(motif_internal_clash_json),
                ),
                cwd=root,
            )
            t3 = time.time()

            run_record["timing_s"] = {
                "candidates_s": float(t1 - t0),
                "solve_s": float(t2 - t1),
                "validate_s": float(t3 - t2),
                "total_s": float(t3 - t0),
            }
        except subprocess.CalledProcessError as e:
            t_err = time.time()
            run_record["status"] = "error"
            run_record["error"] = {
                "step": step,
                "returncode": int(e.returncode),
                "cmd": [str(x) for x in e.cmd] if isinstance(e.cmd, (list, tuple)) else str(e.cmd),
            }
            run_record["timing_s"] = {"total_s": float(t_err - t0)}
            if args.stop_on_error:
                raise

        cand_meta = _read_json_if_exists(cand_prefix.with_suffix(".json")) or {}
        motif_report = _read_json_if_exists(motif_json) or {}
        sat_report = _read_json_if_exists(motif_val_json) or {}
        clash_report = _read_json_if_exists(motif_clash_json) or {}
        internal_clash_report = _read_json_if_exists(motif_internal_clash_json) or {}

        contact_stats = cand_meta.get("pocket_contact_stats", {})
        if (cand_prefix.with_suffix(".npz")).exists():
            with np.load(cand_prefix.with_suffix(".npz"), allow_pickle=True) as z:
                if "lig_sc_contact_count_u8" in z.files:
                    c = z["lig_sc_contact_count_u8"].astype(np.float64)
                    if c.size:
                        contact_stats = {
                            "mean": float(c.mean()),
                            "max": int(c.max()),
                            "min": int(c.min()),
                        }

        coverage_complete = bool(motif_report.get("coverage_complete", False))
        all_satisfied = bool(sat_report.get("all_satisfied", False))
        clash_ok = bool(clash_report.get("ok", False))
        internal_clash_ok = bool(internal_clash_report.get("ok", False))
        n_selected = motif_report.get("n_selected")
        total_s = float(run_record["timing_s"].get("total_s", 0.0))
        contact_mean = float(contact_stats.get("mean", 0.0))
        feasible = bool(run_record["status"] == "ok" and coverage_complete and all_satisfied and clash_ok and internal_clash_ok)

        quality_score = 0.0
        quality_score += 10000.0 if feasible else 0.0
        quality_score += 1000.0 if coverage_complete else 0.0
        quality_score += 500.0 if all_satisfied else 0.0
        quality_score += 200.0 if clash_ok else 0.0
        quality_score += 200.0 if internal_clash_ok else 0.0
        quality_score += 5.0 * contact_mean
        if isinstance(n_selected, int):
            quality_score -= 10.0 * float(n_selected)
        quality_score -= 0.01 * total_s

        run_record["metrics"] = {
            "n_candidates": int(cand_meta.get("n_candidates", 0)),
            "frame_coverage_union_bitmask": int(cand_meta.get("frame_coverage_union_bitmask", 0)),
            "satisfaction_union_bitmask": int(cand_meta.get("satisfaction_union_bitmask", 0)),
            "n_selected": int(n_selected) if isinstance(n_selected, int) else None,
            "coverage_complete": coverage_complete,
            "all_satisfied": all_satisfied,
            "ligand_clash_ok": clash_ok,
            "internal_clash_ok": internal_clash_ok,
            "solver_status": motif_report.get("solver", {}).get("status"),
            "pocket_contact_stats": contact_stats,
            "worst_ligand_overlap": clash_report.get("worst_overlap"),
            "worst_internal_overlap": internal_clash_report.get("worst_overlap"),
        }
        run_record["feasible"] = feasible
        run_record["quality_score"] = float(quality_score)

        runs.append(run_record)

    runs_sorted = sorted(
        runs,
        key=lambda r: (
            int(bool(r.get("feasible"))),
            float(r.get("quality_score", 0.0)),
        ),
        reverse=True,
    )
    best = runs_sorted[0] if runs_sorted else None
    best_feasible = next((r for r in runs_sorted if bool(r.get("feasible"))), None)
    recommended = best_feasible or best

    summary = {
        "tag": args.tag,
        "outdir": str(outdir),
        "inputs": {
            "ligand": str(ligand),
            "polar": str(polar),
            "site_frames": str(frames),
            "vdxform_dir": str(vdx),
        },
        "sweep": {
            "min_lig_donor_angle_list": min_lig_donor_angles,
            "pocket_contact_weight_list": pocket_contact_weights,
            "min_pocket_sidechain_contacts_list": min_pocket_sidechain_contacts,
            "n_runs": len(runs),
        },
        "recommended_run_id": recommended.get("run_id") if recommended else None,
        "recommended_feasible": bool(recommended.get("feasible")) if recommended else False,
        "runs": runs_sorted,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        f"# MTX Pocket Knob Sweep ({args.tag})",
        "",
        f"- outdir: `{outdir}`",
        f"- recommended: `{summary['recommended_run_id']}` (feasible={summary['recommended_feasible']})",
        "",
        "| run_id | feasible | n_candidates | n_selected | all_satisfied | clash_ok | internal_ok | contact_mean | total_s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in runs_sorted:
        m = r.get("metrics", {})
        t = r.get("timing_s", {})
        contact_mean = (m.get("pocket_contact_stats") or {}).get("mean", 0.0)
        md_lines.append(
            "| "
            + f"{r.get('run_id')} | {int(bool(r.get('feasible')))} | {m.get('n_candidates', 0)} | "
            + f"{m.get('n_selected')} | {int(bool(m.get('all_satisfied')))} | {int(bool(m.get('ligand_clash_ok')))} | "
            + f"{int(bool(m.get('internal_clash_ok')))} | {float(contact_mean):.3f} | {float(t.get('total_s', 0.0)):.2f} |"
        )
    (outdir / "summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    if not any(bool(r.get("feasible")) for r in runs):
        raise SystemExit(f"No feasible run found. See {outdir/'summary.json'}")


if __name__ == "__main__":
    main()
