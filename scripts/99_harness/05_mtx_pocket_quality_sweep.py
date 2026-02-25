#!/usr/bin/env python
from __future__ import annotations

import argparse
import itertools
import json
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any


def _run(cmd: list[str], *, cwd: Path, log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"$ {' '.join(shlex.quote(x) for x in cmd)}\n")
        log.flush()
        subprocess.run(cmd, check=True, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT, env=env)


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
    return f"{float(v):g}".replace("-", "m").replace(".", "p")


def _safe_token(v: str) -> str:
    t = "".join(ch if ch.isalnum() else "_" for ch in v.strip())
    return t or "MTX"


def _read_json_if_exists(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _py_cmd(python_bin: Path, script: Path, *args: str) -> list[str]:
    return [str(python_bin), str(script), *args]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Sweep pocket-quality knobs for candidate generation and solver, "
            "then summarize feasibility/quality/runtime."
        )
    )
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--python-bin", type=Path, default=None, help="Python interpreter used to run pipeline scripts.")
    ap.add_argument("--ligand-tag", type=str, default="MTX", help="Prefix token for output filenames and outdir names.")
    ap.add_argument("--ligand-resname", type=str, default="MTX", help="Residue name used for clash validators.")
    ap.add_argument("--ligand-pdb", type=Path, default=Path("inputs/01_cgmap/MTX.pdb"))
    ap.add_argument("--polar-sites", type=Path, default=Path("outputs/02_polar_sites/MTX_polar_sites.json"))
    ap.add_argument("--site-frames", type=Path, default=Path("outputs/02_polar_sites/MTX_site_frames.json"))
    ap.add_argument("--vdxform-dir", type=Path, default=Path("processed/03_vdxform_full"))
    ap.add_argument("--solver", type=str, default="cp_sat", choices=["cp_sat", "greedy"])
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--top-per-site-per-atom", type=int, default=50)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument("--ca-prefilter", type=float, default=14.0, help="CA prefilter radius used for solver clash graph.")
    ap.add_argument("--acceptor-model", type=str, default="plip", choices=["legacy", "plip"])
    ap.add_argument("--score-w-prior", type=float, default=1.0)
    ap.add_argument("--score-w-coverage-list", type=str, default="0.03,0.05")
    ap.add_argument("--score-w-contact-list", type=str, default="0.1,0.2")
    ap.add_argument("--score-w-shell-list", type=str, default="0.2,0.4")
    ap.add_argument("--target-res-list", type=str, default="10,12")
    ap.add_argument("--min-cover-per-polar-list", type=str, default="1")
    ap.add_argument("--contact-dist-min", type=float, default=2.6)
    ap.add_argument("--contact-dist-max", type=float, default=4.6)
    ap.add_argument("--shell-center-dist", type=float, default=3.4)
    ap.add_argument("--shell-half-width", type=float, default=1.8)
    ap.add_argument("--min-sidechain-facing-dot", type=float, default=0.2)
    ap.add_argument("--min-sidechain-centroid-dot", type=float, default=0.0)
    ap.add_argument("--run-plip", action="store_true")
    ap.add_argument("--require-plip-success", action="store_true")
    ap.add_argument(
        "--plip-bin",
        type=str,
        default="plip",
        help="PLIP executable path/name passed to validate_motif_plip.py.",
    )
    ap.add_argument("--max-runs", type=int, default=0, help="If >0, run only first N combinations.")
    ap.add_argument("--stop-on-error", action="store_true")
    args = ap.parse_args()

    root = args.repo_root.resolve()
    if args.python_bin is None:
        default_py = root / ".venv" / "bin" / "python"
        python_bin = default_py if default_py.exists() else Path("python3")
    else:
        python_bin = args.python_bin if args.python_bin.is_absolute() else (root / args.python_bin)
    if not python_bin.exists():
        raise FileNotFoundError(f"Python binary not found: {python_bin}")

    run_env = os.environ.copy()
    py_bin_dir = str(python_bin.parent.resolve())
    run_env["PATH"] = py_bin_dir + os.pathsep + run_env.get("PATH", "")

    ligand_tag = _safe_token(str(args.ligand_tag))
    ligand_resname = str(args.ligand_resname).strip() or "MTX"
    outdir = root / "processed/99_harness" / f"{ligand_tag.lower()}_pocket_quality_sweep_{args.tag}"
    runs_dir = outdir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    ligand = (root / args.ligand_pdb).resolve() if not args.ligand_pdb.is_absolute() else args.ligand_pdb.resolve()
    polar = (root / args.polar_sites).resolve() if not args.polar_sites.is_absolute() else args.polar_sites.resolve()
    frames = (root / args.site_frames).resolve() if not args.site_frames.is_absolute() else args.site_frames.resolve()
    vdx = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()

    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        _ensure_file(
            vdx / cg / f"vdxform_{cg}.npz",
            hint_cmd="bash scripts/03_vdxform/02_run_mtx_needed_cgs.sh (full) or 03_run_mtx_needed_cgs_debug.sh (debug)",
        )
    _ensure_file(ligand, hint_cmd="Provide --ligand-pdb, e.g. inputs/01_cgmap/<LIG>.pdb")
    _ensure_file(polar, hint_cmd="Provide --polar-sites, e.g. outputs/02_polar_sites/<LIG>_polar_sites.json")
    _ensure_file(frames, hint_cmd="Provide --site-frames, e.g. outputs/02_polar_sites/<LIG>_site_frames.json")

    candidate_py = root / "scripts/04_candidates/01_generate_candidates.py"
    solver_py = root / "scripts/05_solver/01_solve_motif.py"
    validate_sat_py = root / "scripts/05_solver/03_validate_motif_polar_satisfaction.py"
    validate_clash_py = root / "scripts/05_solver/04_validate_motif_clashes.py"
    validate_internal_clash_py = root / "scripts/05_solver/05_validate_motif_internal_clashes.py"
    validate_plip_py = root / "scripts/05_solver/06_validate_motif_plip.py"

    for path in [candidate_py, solver_py, validate_sat_py, validate_clash_py, validate_internal_clash_py]:
        _ensure_file(path, hint_cmd="sync repository scripts")
    if args.run_plip:
        _ensure_file(validate_plip_py, hint_cmd="sync repository scripts")

    score_w_coverage_list = _parse_float_list(args.score_w_coverage_list)
    score_w_contact_list = _parse_float_list(args.score_w_contact_list)
    score_w_shell_list = _parse_float_list(args.score_w_shell_list)
    target_res_list = _parse_int_list(args.target_res_list)
    min_cover_list = _parse_int_list(args.min_cover_per_polar_list)

    combos = list(itertools.product(score_w_coverage_list, score_w_contact_list, score_w_shell_list, target_res_list, min_cover_list))
    if int(args.max_runs) > 0:
        combos = combos[: int(args.max_runs)]
    if not combos:
        raise ValueError("No combinations to run.")

    runs: list[dict[str, Any]] = []
    for i, (w_cov, w_contact, w_shell, target_res, min_cover) in enumerate(combos, start=1):
        run_id = (
            f"run_{i:02d}_wcov{_safe_num(w_cov)}_wct{_safe_num(w_contact)}_"
            f"wsh{_safe_num(w_shell)}_tr{_safe_num(target_res)}_mc{_safe_num(min_cover)}"
        )
        run_dir = runs_dir / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        log_path = run_dir / "pipeline.log"

        cand_prefix = run_dir / f"{ligand_tag}_candidates"
        motif_json = run_dir / f"{ligand_tag}_motif.json"
        motif_pdb = run_dir / f"{ligand_tag}_motif.pdb"
        motif_val_json = run_dir / f"{ligand_tag}_motif_validation.json"
        motif_clash_json = run_dir / f"{ligand_tag}_motif_clash_validation.json"
        motif_internal_clash_json = run_dir / f"{ligand_tag}_motif_internal_clash_validation.json"
        motif_plip_json = run_dir / f"{ligand_tag}_motif_plip_validation.json"

        run_record: dict[str, Any] = {
            "run_id": run_id,
            "status": "ok",
            "error": None,
            "plip_error": None,
            "params": {
                "score_w_prior": float(args.score_w_prior),
                "score_w_coverage": float(w_cov),
                "score_w_contact": float(w_contact),
                "score_w_shell": float(w_shell),
                "target_res": int(target_res),
                "min_cover_per_polar": int(min_cover),
                "solver": str(args.solver),
                "top_per_site": int(args.top_per_site),
                "top_per_site_per_atom": int(args.top_per_site_per_atom),
                "chunk_size": int(args.chunk_size),
                "time_limit_s": float(args.time_limit_s),
                "ca_prefilter": float(args.ca_prefilter),
                "acceptor_model": str(args.acceptor_model),
                "run_plip": bool(args.run_plip),
                "plip_bin": str(args.plip_bin),
            },
            "paths": {
                "run_dir": str(run_dir),
                "log": str(log_path),
                "candidates_npz": str(cand_prefix.with_suffix(".npz")),
                "candidates_json": str(cand_prefix.with_suffix(".json")),
                "motif_json": str(motif_json),
                "motif_pdb": str(motif_pdb),
                "sat_json": str(motif_val_json),
                "clash_json": str(motif_clash_json),
                "internal_clash_json": str(motif_internal_clash_json),
                "plip_json": str(motif_plip_json),
            },
            "timing_s": {},
        }

        step = "candidates"
        t0 = time.time()
        try:
            _run(
                _py_cmd(
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
                    "--clash-tol",
                    str(float(args.clash_tol)),
                    "--exclude-aa3",
                    "PRO,CYS",
                    "--acceptor-model",
                    str(args.acceptor_model),
                    "--require-sidechain-facing",
                    "--min-sidechain-facing-dot",
                    str(float(args.min_sidechain_facing_dot)),
                    "--min-sidechain-centroid-dot",
                    str(float(args.min_sidechain_centroid_dot)),
                    "--score-w-prior",
                    str(float(args.score_w_prior)),
                    "--score-w-coverage",
                    str(float(w_cov)),
                    "--score-w-contact",
                    str(float(w_contact)),
                    "--score-w-shell",
                    str(float(w_shell)),
                    "--contact-dist-min",
                    str(float(args.contact_dist_min)),
                    "--contact-dist-max",
                    str(float(args.contact_dist_max)),
                    "--shell-center-dist",
                    str(float(args.shell_center_dist)),
                    "--shell-half-width",
                    str(float(args.shell_half_width)),
                    "--require-full-coverage",
                ),
                cwd=root,
                log_path=log_path,
                env=run_env,
            )
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
                    "--target-res",
                    str(int(target_res)),
                    "--min-cover-per-polar",
                    str(int(min_cover)),
                    "--time-limit-s",
                    str(float(args.time_limit_s)),
                    "--num-workers",
                    "1",
                    "--grid-size",
                    "4.0",
                    "--ca-prefilter",
                    str(float(args.ca_prefilter)),
                    "--clash-tol",
                    str(float(args.clash_tol)),
                ),
                cwd=root,
                log_path=log_path,
                env=run_env,
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
                    "--ligand-resname",
                    str(ligand_resname),
                    "--acceptor-model",
                    str(args.acceptor_model),
                    "-o",
                    str(motif_val_json),
                ),
                cwd=root,
                log_path=log_path,
                env=run_env,
            )
            step = "validate_ligand_clash"
            _run(
                _py_cmd(
                    python_bin,
                    validate_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    str(ligand_resname),
                    "-o",
                    str(motif_clash_json),
                    "--tol",
                    "0.5",
                    "--fail-overlap",
                    "0.5",
                ),
                cwd=root,
                log_path=log_path,
                env=run_env,
            )
            step = "validate_internal_clash"
            _run(
                _py_cmd(
                    python_bin,
                    validate_internal_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    str(ligand_resname),
                    "-o",
                    str(motif_internal_clash_json),
                    "--tol",
                    "0.5",
                    "--fail-overlap",
                    "0.5",
                ),
                cwd=root,
                log_path=log_path,
                env=run_env,
            )
            if args.run_plip:
                step = "validate_plip"
                try:
                    _run(
                        _py_cmd(
                            python_bin,
                            validate_plip_py,
                            "--motif-pdb",
                            str(motif_pdb),
                            "--polar-sites",
                            str(polar),
                            "-o",
                            str(motif_plip_json),
                            "--timeout-s",
                            "180",
                            "--plip-bin",
                            str(args.plip_bin),
                        ),
                        cwd=root,
                        log_path=log_path,
                        env=run_env,
                    )
                except subprocess.CalledProcessError as e:
                    run_record["plip_error"] = {
                        "step": step,
                        "returncode": int(e.returncode),
                        "cmd": [str(x) for x in e.cmd] if isinstance(e.cmd, (list, tuple)) else str(e.cmd),
                    }
                    if args.require_plip_success:
                        raise
                    if args.stop_on_error:
                        raise
                    with log_path.open("a", encoding="utf-8") as log:
                        log.write("[warn] Optional PLIP validation failed; continuing because --require-plip-success is not set.\n")
            t3 = time.time()
            run_record["timing_s"] = {
                "candidates_s": float(t1 - t0),
                "solve_s": float(t2 - t1),
                "validate_s": float(t3 - t2),
                "total_s": float(t3 - t0),
            }
        except subprocess.CalledProcessError as e:
            run_record["status"] = "error"
            run_record["error"] = {
                "step": step,
                "returncode": int(e.returncode),
                "cmd": [str(x) for x in e.cmd] if isinstance(e.cmd, (list, tuple)) else str(e.cmd),
            }
            run_record["timing_s"] = {"total_s": float(time.time() - t0)}
            if args.stop_on_error:
                raise

        cand_meta = _read_json_if_exists(cand_prefix.with_suffix(".json")) or {}
        motif_report = _read_json_if_exists(motif_json) or {}
        sat_report = _read_json_if_exists(motif_val_json) or {}
        clash_report = _read_json_if_exists(motif_clash_json) or {}
        internal_clash_report = _read_json_if_exists(motif_internal_clash_json) or {}
        plip_report = _read_json_if_exists(motif_plip_json) or {}

        coverage_complete = bool(motif_report.get("coverage_complete", False))
        all_satisfied = bool(sat_report.get("all_satisfied", False))
        ligand_clash_ok = bool(clash_report.get("ok", False))
        internal_clash_ok = bool(internal_clash_report.get("ok", False))
        plip_all_satisfied = None
        if args.run_plip and "all_satisfied" in plip_report:
            plip_all_satisfied = bool(plip_report.get("all_satisfied"))
        n_selected = motif_report.get("n_selected")
        n_candidates = int(cand_meta.get("n_candidates", 0))
        total_s = float(run_record["timing_s"].get("total_s", 0.0))

        feasible_core = bool(run_record["status"] == "ok" and coverage_complete and all_satisfied and ligand_clash_ok and internal_clash_ok)
        feasible = feasible_core
        if args.run_plip and args.require_plip_success:
            feasible = feasible and bool(plip_all_satisfied)

        quality_score = 0.0
        quality_score += 10000.0 if feasible else 0.0
        quality_score += 2000.0 if feasible_core else 0.0
        quality_score += 1000.0 if coverage_complete else 0.0
        quality_score += 500.0 if all_satisfied else 0.0
        quality_score += 300.0 if ligand_clash_ok else 0.0
        quality_score += 300.0 if internal_clash_ok else 0.0
        if args.run_plip:
            quality_score += 400.0 if bool(plip_all_satisfied) else 0.0
            quality_score += 5.0 * float(plip_report.get("n_satisfied", 0))
        if isinstance(n_selected, int):
            quality_score -= 20.0 * abs(int(n_selected) - int(target_res))
        quality_score -= 0.05 * total_s

        run_record["metrics"] = {
            "n_candidates": n_candidates,
            "n_selected": int(n_selected) if isinstance(n_selected, int) else None,
            "coverage_complete": coverage_complete,
            "all_satisfied": all_satisfied,
            "ligand_clash_ok": ligand_clash_ok,
            "internal_clash_ok": internal_clash_ok,
            "plip_all_satisfied": plip_all_satisfied,
            "plip_n_satisfied": plip_report.get("n_satisfied") if args.run_plip and "n_satisfied" in plip_report else None,
            "plip_n_polar_atoms": plip_report.get("n_polar_atoms") if args.run_plip and "n_polar_atoms" in plip_report else None,
            "worst_ligand_overlap": clash_report.get("worst_overlap"),
            "worst_internal_overlap": internal_clash_report.get("worst_overlap"),
            "solver_status": motif_report.get("solver", {}).get("status"),
            "size_deviation": motif_report.get("solver", {}).get("size_deviation"),
        }
        run_record["feasible"] = bool(feasible)
        run_record["quality_score"] = float(quality_score)
        runs.append(run_record)

    runs_sorted = sorted(runs, key=lambda r: (int(bool(r.get("feasible"))), float(r.get("quality_score", 0.0))), reverse=True)
    best = runs_sorted[0] if runs_sorted else None
    best_feasible = next((r for r in runs_sorted if bool(r.get("feasible"))), None)
    recommended = best_feasible or best

    summary = {
        "tag": args.tag,
        "outdir": str(outdir),
        "inputs": {
            "ligand_tag": str(ligand_tag),
            "ligand_resname": str(ligand_resname),
            "ligand_pdb": str(ligand),
            "polar_sites": str(polar),
            "site_frames": str(frames),
            "vdxform_dir": str(vdx),
        },
        "sweep": {
            "score_w_coverage_list": score_w_coverage_list,
            "score_w_contact_list": score_w_contact_list,
            "score_w_shell_list": score_w_shell_list,
            "target_res_list": target_res_list,
            "min_cover_per_polar_list": min_cover_list,
            "n_runs": len(runs),
        },
        "run_plip": bool(args.run_plip),
        "require_plip_success": bool(args.require_plip_success),
        "recommended_run_id": recommended.get("run_id") if recommended else None,
        "recommended_feasible": bool(recommended.get("feasible")) if recommended else False,
        "runs": runs_sorted,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        f"# {ligand_tag} Pocket Quality Sweep ({args.tag})",
        "",
        f"- outdir: `{outdir}`",
        f"- recommended: `{summary['recommended_run_id']}` (feasible={summary['recommended_feasible']})",
        "",
        "| run_id | feasible | n_candidates | n_selected | sat | plip | lig_clash | int_clash | total_s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in runs_sorted:
        m = r.get("metrics", {})
        t = r.get("timing_s", {})
        plip_ok = m.get("plip_all_satisfied")
        plip_val = "-" if plip_ok is None else str(int(bool(plip_ok)))
        md_lines.append(
            "| "
            + f"{r.get('run_id')} | {int(bool(r.get('feasible')))} | {m.get('n_candidates', 0)} | "
            + f"{m.get('n_selected')} | {int(bool(m.get('all_satisfied')))} | {plip_val} | "
            + f"{int(bool(m.get('ligand_clash_ok')))} | {int(bool(m.get('internal_clash_ok')))} | "
            + f"{float(t.get('total_s', 0.0)):.2f} |"
        )
    (outdir / "summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    if not any(bool(r.get("feasible")) for r in runs):
        raise SystemExit(f"No feasible run found. See {outdir / 'summary.json'}")


if __name__ == "__main__":
    main()
