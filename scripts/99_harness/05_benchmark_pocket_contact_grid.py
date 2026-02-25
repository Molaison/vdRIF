#!/usr/bin/env python
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import time
from pathlib import Path
from typing import Any


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd))


def _py_cmd(python_bin: str, script: Path, *args: str) -> list[str]:
    return [python_bin, str(script), *args]


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


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Sweep MTX pocket-contact knobs with full vdxform and report PLIP/clash/coverage/runtime metrics."
        )
    )
    default_script_root = Path(__file__).resolve().parents[2]
    ap.add_argument("--script-root", type=Path, default=default_script_root, help="Repo root that provides scripts/.")
    ap.add_argument(
        "--data-root",
        type=Path,
        default=default_script_root,
        help="Repo root that provides inputs/, outputs/, processed/. Can differ from --script-root.",
    )
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--outdir", type=Path, default=None, help="Optional explicit output directory.")
    ap.add_argument(
        "--python-bin",
        type=str,
        default="python",
        help="Python executable used to run candidate/solver/validator scripts.",
    )
    ap.add_argument(
        "--vdxform-dir",
        type=Path,
        default=Path("processed/03_vdxform_full"),
        help="Directory containing per-cg vdxform_<cg>.npz files (absolute or relative to --data-root).",
    )

    ap.add_argument("--solver", type=str, default="greedy", choices=["cp_sat", "greedy"])
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument(
        "--top-per-site-list",
        type=str,
        default="",
        help="Optional CSV int list to sweep top_per_site. If empty, uses --top-per-site only.",
    )
    ap.add_argument("--top-per-site-per-atom", type=int, default=50)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument("--acceptor-model", type=str, default="plip", choices=["legacy", "plip"])
    ap.add_argument("--min-sidechain-contact-dist", type=float, default=2.8)
    ap.add_argument("--max-sidechain-contact-dist", type=float, default=4.8)

    ap.add_argument(
        "--min-sidechain-contact-count-list",
        type=str,
        default="1,2",
        help="CSV int list for --min-sidechain-contact-count.",
    )
    ap.add_argument(
        "--sidechain-contact-weight-list",
        type=str,
        default="0.05,0.10,0.15",
        help="CSV float list for --sidechain-contact-weight.",
    )
    ap.add_argument("--max-runs", type=int, default=0, help="If >0, only run first N combinations.")
    ap.add_argument("--stop-on-error", action="store_true")
    ap.add_argument("--keep-plip-outdir", action="store_true")
    ap.add_argument(
        "--enable-plip-fill",
        action="store_true",
        help="Run PLIP-guided post-solver residue filling when baseline PLIP coverage is incomplete.",
    )
    ap.add_argument("--plip-fill-max-res", type=int, default=15)
    ap.add_argument("--plip-fill-top-try-per-atom", type=int, default=120)
    args = ap.parse_args()

    if float(args.max_sidechain_contact_dist) <= float(args.min_sidechain_contact_dist):
        raise ValueError(
            "--max-sidechain-contact-dist must be > --min-sidechain-contact-dist, "
            f"got min={args.min_sidechain_contact_dist} max={args.max_sidechain_contact_dist}"
        )

    script_root = args.script_root.resolve()
    data_root = args.data_root.resolve()

    outdir = args.outdir.resolve() if args.outdir else (data_root / "processed/99_harness" / f"mtx_pocket_contact_grid_{args.tag}")
    runs_dir = outdir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    ligand = data_root / "inputs/01_cgmap/MTX.pdb"
    polar = data_root / "outputs/02_polar_sites/MTX_polar_sites.json"
    frames = data_root / "outputs/02_polar_sites/MTX_site_frames.json"

    vdx = args.vdxform_dir.resolve() if args.vdxform_dir.is_absolute() else (data_root / args.vdxform_dir).resolve()
    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        _ensure_file(
            vdx / cg / f"vdxform_{cg}.npz",
            hint_cmd="bash scripts/03_vdxform/03_run_mtx_needed_cgs_debug.sh (debug) or 02_run_mtx_needed_cgs.sh (full)",
        )
    _ensure_file(ligand, hint_cmd="Place MTX at inputs/01_cgmap/MTX.pdb")
    _ensure_file(polar, hint_cmd="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")
    _ensure_file(frames, hint_cmd="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")

    candidate_py = script_root / "scripts/04_candidates/01_generate_candidates.py"
    solver_py = script_root / "scripts/05_solver/01_solve_motif.py"
    validate_sat_py = script_root / "scripts/05_solver/03_validate_motif_polar_satisfaction.py"
    validate_clash_py = script_root / "scripts/05_solver/04_validate_motif_clashes.py"
    validate_internal_clash_py = script_root / "scripts/05_solver/05_validate_motif_internal_clashes.py"
    validate_plip_py = script_root / "scripts/05_solver/06_validate_motif_plip.py"
    plip_fill_py = script_root / "scripts/05_solver/07_plip_fill_motif.py"

    for p in [
        candidate_py,
        solver_py,
        validate_sat_py,
        validate_clash_py,
        validate_internal_clash_py,
        validate_plip_py,
        plip_fill_py,
    ]:
        _ensure_file(p, hint_cmd="Ensure scripts are present in --script-root")

    min_contact_counts = _parse_int_list(args.min_sidechain_contact_count_list)
    contact_weights = _parse_float_list(args.sidechain_contact_weight_list)

    top_per_site_list = (
        _parse_int_list(args.top_per_site_list) if str(args.top_per_site_list).strip() else [int(args.top_per_site)]
    )
    combos = list(itertools.product(top_per_site_list, min_contact_counts, contact_weights))
    if int(args.max_runs) > 0:
        combos = combos[: int(args.max_runs)]
    if not combos:
        raise ValueError("No parameter combinations to run.")

    runs: list[dict[str, Any]] = []
    for i, (top_per_site_run, min_sc_contacts, sc_weight) in enumerate(combos, start=1):
        run_id = (
            f"run_{i:02d}_top{_safe_num(top_per_site_run)}_"
            f"cnt{_safe_num(min_sc_contacts)}_w{_safe_num(sc_weight)}"
        )
        run_dir = runs_dir / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        cand_prefix = run_dir / "MTX_candidates"
        motif_json = run_dir / "MTX_motif.json"
        motif_pdb = run_dir / "MTX_motif.pdb"
        motif_val_json = run_dir / "MTX_motif_validation.json"
        motif_clash_json = run_dir / "MTX_motif_clash_validation.json"
        motif_internal_clash_json = run_dir / "MTX_motif_internal_clash_validation.json"
        motif_plip_json = run_dir / "MTX_motif_plip_validation.json"
        motif_fill_pdb = run_dir / "MTX_motif_plip_fill.pdb"
        motif_fill_report_json = run_dir / "MTX_motif_plip_fill_report.json"
        motif_fill_val_json = run_dir / "MTX_motif_plip_fill_validation.json"
        motif_fill_clash_json = run_dir / "MTX_motif_plip_fill_clash_validation.json"
        motif_fill_internal_clash_json = run_dir / "MTX_motif_plip_fill_internal_clash_validation.json"
        motif_fill_plip_json = run_dir / "MTX_motif_plip_fill_plip_validation.json"

        # Metric defaults; may switch to PLIP-fill outputs if improved.
        final_motif_pdb = motif_pdb
        final_val_json = motif_val_json
        final_clash_json = motif_clash_json
        final_internal_clash_json = motif_internal_clash_json
        final_plip_json = motif_plip_json
        plip_fill_used = False

        run_record: dict[str, Any] = {
            "run_id": run_id,
            "params": {
                "min_sidechain_contact_count": int(min_sc_contacts),
                "sidechain_contact_weight": float(sc_weight),
                "min_sidechain_contact_dist": float(args.min_sidechain_contact_dist),
                "max_sidechain_contact_dist": float(args.max_sidechain_contact_dist),
                "acceptor_model": str(args.acceptor_model),
                "top_per_site": int(top_per_site_run),
                "top_per_site_per_atom": int(args.top_per_site_per_atom),
                "chunk_size": int(args.chunk_size),
                "solver": str(args.solver),
                "time_limit_s": float(args.time_limit_s),
                "clash_tol": float(args.clash_tol),
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
                "motif_plip_validation_json": str(motif_plip_json),
                "motif_plip_fill_pdb": str(motif_fill_pdb),
                "motif_plip_fill_report_json": str(motif_fill_report_json),
                "motif_plip_fill_validation_json": str(motif_fill_val_json),
                "motif_plip_fill_clash_validation_json": str(motif_fill_clash_json),
                "motif_plip_fill_internal_clash_validation_json": str(motif_fill_internal_clash_json),
                "motif_plip_fill_plip_validation_json": str(motif_fill_plip_json),
            },
            "status": "ok",
            "timing_s": {},
            "error": None,
        }

        t0 = time.time()
        step = "candidates"
        try:
            _run(
                _py_cmd(
                    str(args.python_bin),
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
                    str(int(top_per_site_run)),
                    "--top-per-site-per-atom",
                    str(int(args.top_per_site_per_atom)),
                    "--acceptor-model",
                    str(args.acceptor_model),
                    "--min-sidechain-contact-dist",
                    str(float(args.min_sidechain_contact_dist)),
                    "--max-sidechain-contact-dist",
                    str(float(args.max_sidechain_contact_dist)),
                    "--min-sidechain-contact-count",
                    str(int(min_sc_contacts)),
                    "--sidechain-contact-weight",
                    str(float(sc_weight)),
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
                ),
                cwd=data_root,
            )
            t1 = time.time()

            step = "solver"
            _run(
                _py_cmd(
                    str(args.python_bin),
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
                cwd=data_root,
            )
            t2 = time.time()

            step = "validate_sat"
            _run(
                _py_cmd(
                    str(args.python_bin),
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
                cwd=data_root,
            )

            step = "validate_ligand_clash"
            _run(
                _py_cmd(
                    str(args.python_bin),
                    validate_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    "MTX",
                    "-o",
                    str(motif_clash_json),
                ),
                cwd=data_root,
            )

            step = "validate_internal_clash"
            _run(
                _py_cmd(
                    str(args.python_bin),
                    validate_internal_clash_py,
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    "MTX",
                    "-o",
                    str(motif_internal_clash_json),
                ),
                cwd=data_root,
            )

            step = "validate_plip"
            plip_cmd = _py_cmd(
                str(args.python_bin),
                validate_plip_py,
                "--motif-pdb",
                str(motif_pdb),
                "--polar-sites",
                str(polar),
                "-o",
                str(motif_plip_json),
            )
            if args.keep_plip_outdir:
                plip_cmd.append("--keep-plip-outdir")
            _run(plip_cmd, cwd=data_root)

            plip_before = _read_json_if_exists(motif_plip_json) or {}
            if bool(args.enable_plip_fill) and not bool(plip_before.get("all_satisfied", False)):
                step = "plip_fill"
                _run(
                    _py_cmd(
                        str(args.python_bin),
                        plip_fill_py,
                        "--motif-pdb",
                        str(motif_pdb),
                        "--candidates-npz",
                        str(cand_prefix.with_suffix(".npz")),
                        "--candidates-meta",
                        str(cand_prefix.with_suffix(".json")),
                        "--polar-sites",
                        str(polar),
                        "--ligand-pdb",
                        str(ligand),
                        "--out-pdb",
                        str(motif_fill_pdb),
                        "--out-report-json",
                        str(motif_fill_report_json),
                        "--max-res",
                        str(int(args.plip_fill_max_res)),
                        "--clash-tol",
                        str(float(args.clash_tol)),
                        "--top-try-per-atom",
                        str(int(args.plip_fill_top_try_per_atom)),
                    ),
                    cwd=data_root,
                )

                step = "validate_sat_post_fill"
                _run(
                    _py_cmd(
                        str(args.python_bin),
                        validate_sat_py,
                        "--polar-sites",
                        str(polar),
                        "--motif-pdb",
                        str(motif_fill_pdb),
                        "--acceptor-model",
                        str(args.acceptor_model),
                        "-o",
                        str(motif_fill_val_json),
                    ),
                    cwd=data_root,
                )

                step = "validate_ligand_clash_post_fill"
                _run(
                    _py_cmd(
                        str(args.python_bin),
                        validate_clash_py,
                        "--motif-pdb",
                        str(motif_fill_pdb),
                        "--ligand-resname",
                        "MTX",
                        "-o",
                        str(motif_fill_clash_json),
                    ),
                    cwd=data_root,
                )

                step = "validate_internal_clash_post_fill"
                _run(
                    _py_cmd(
                        str(args.python_bin),
                        validate_internal_clash_py,
                        "--motif-pdb",
                        str(motif_fill_pdb),
                        "--ligand-resname",
                        "MTX",
                        "-o",
                        str(motif_fill_internal_clash_json),
                    ),
                    cwd=data_root,
                )

                step = "validate_plip_post_fill"
                fill_plip_cmd = _py_cmd(
                    str(args.python_bin),
                    validate_plip_py,
                    "--motif-pdb",
                    str(motif_fill_pdb),
                    "--polar-sites",
                    str(polar),
                    "-o",
                    str(motif_fill_plip_json),
                )
                if args.keep_plip_outdir:
                    fill_plip_cmd.append("--keep-plip-outdir")
                _run(fill_plip_cmd, cwd=data_root)

                plip_after = _read_json_if_exists(motif_fill_plip_json) or {}
                n_sat_before = int(plip_before.get("n_satisfied", -1))
                n_sat_after = int(plip_after.get("n_satisfied", -1))
                plip_fill_used = n_sat_after >= n_sat_before
                run_record["plip_fill"] = {
                    "attempted": True,
                    "used_post_fill": bool(plip_fill_used),
                    "all_satisfied_pre_fill": bool(plip_before.get("all_satisfied", False)),
                    "all_satisfied_post_fill": bool(plip_after.get("all_satisfied", False)),
                    "n_satisfied_pre_fill": n_sat_before,
                    "n_satisfied_post_fill": n_sat_after,
                    "unsatisfied_pre_fill": plip_before.get("unsatisfied_polar_atoms", []),
                    "unsatisfied_post_fill": plip_after.get("unsatisfied_polar_atoms", []),
                }
                if plip_fill_used:
                    final_motif_pdb = motif_fill_pdb
                    final_val_json = motif_fill_val_json
                    final_clash_json = motif_fill_clash_json
                    final_internal_clash_json = motif_fill_internal_clash_json
                    final_plip_json = motif_fill_plip_json

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

        run_record["paths"]["final_motif_pdb_for_metrics"] = str(final_motif_pdb)
        run_record["paths"]["final_motif_validation_json_for_metrics"] = str(final_val_json)
        run_record["paths"]["final_motif_clash_json_for_metrics"] = str(final_clash_json)
        run_record["paths"]["final_motif_internal_clash_json_for_metrics"] = str(final_internal_clash_json)
        run_record["paths"]["final_motif_plip_json_for_metrics"] = str(final_plip_json)

        cand_meta = _read_json_if_exists(cand_prefix.with_suffix(".json")) or {}
        motif_report = _read_json_if_exists(motif_json) or {}
        sat_report = _read_json_if_exists(final_val_json) or {}
        clash_report = _read_json_if_exists(final_clash_json) or {}
        internal_clash_report = _read_json_if_exists(final_internal_clash_json) or {}
        plip_report = _read_json_if_exists(final_plip_json) or {}
        plip_report_pre = _read_json_if_exists(motif_plip_json) or {}

        coverage_complete = bool(motif_report.get("coverage_complete", False))
        all_satisfied_simple = bool(sat_report.get("all_satisfied", False))
        plip_all_satisfied = bool(plip_report.get("all_satisfied", False))
        clash_ok = bool(clash_report.get("ok", False))
        internal_clash_ok = bool(internal_clash_report.get("ok", False))
        n_selected = motif_report.get("n_selected")
        total_s = float(run_record["timing_s"].get("total_s", 0.0))

        feasible = bool(
            run_record["status"] == "ok"
            and coverage_complete
            and clash_ok
            and internal_clash_ok
            and plip_all_satisfied
        )

        quality_score = 0.0
        quality_score += 10000.0 if feasible else 0.0
        quality_score += 2000.0 if plip_all_satisfied else 0.0
        quality_score += 1000.0 if coverage_complete else 0.0
        quality_score += 500.0 if clash_ok else 0.0
        quality_score += 500.0 if internal_clash_ok else 0.0
        if isinstance(n_selected, int):
            quality_score -= 5.0 * float(n_selected)
        quality_score -= 0.01 * total_s

        run_record["metrics"] = {
            "n_candidates": int(cand_meta.get("n_candidates", 0)),
            "frame_coverage_union_bitmask": int(cand_meta.get("frame_coverage_union_bitmask", 0)),
            "satisfaction_union_bitmask": int(cand_meta.get("satisfaction_union_bitmask", 0)),
            "n_selected": int(n_selected) if isinstance(n_selected, int) else None,
            "coverage_complete": coverage_complete,
            "all_satisfied_simple": all_satisfied_simple,
            "all_satisfied_plip": plip_all_satisfied,
            "all_satisfied_plip_pre_fill": bool(plip_report_pre.get("all_satisfied", False)),
            "plip_fill_used_for_metrics": bool(plip_fill_used),
            "n_satisfied_plip": plip_report.get("n_satisfied"),
            "n_satisfied_plip_pre_fill": plip_report_pre.get("n_satisfied"),
            "n_polar_atoms_plip": plip_report.get("n_polar_atoms"),
            "unsatisfied_polar_atoms_plip": plip_report.get("unsatisfied_polar_atoms", []),
            "unsatisfied_polar_atoms_plip_pre_fill": plip_report_pre.get("unsatisfied_polar_atoms", []),
            "ligand_clash_ok": clash_ok,
            "internal_clash_ok": internal_clash_ok,
            "worst_ligand_overlap": clash_report.get("worst_overlap"),
            "worst_internal_overlap": internal_clash_report.get("worst_overlap"),
            "solver_status": motif_report.get("solver", {}).get("status"),
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

    plip_sat_count = sum(int(bool((r.get("metrics") or {}).get("all_satisfied_plip"))) for r in runs)
    plip_sat_rate = float(plip_sat_count / len(runs)) if runs else 0.0

    summary: dict[str, Any] = {
        "tag": args.tag,
        "script_root": str(script_root),
        "data_root": str(data_root),
        "outdir": str(outdir),
        "inputs": {
            "ligand": str(ligand),
            "polar": str(polar),
            "site_frames": str(frames),
            "vdxform_dir": str(vdx),
        },
        "sweep": {
            "top_per_site_list": top_per_site_list,
            "min_sidechain_contact_count_list": min_contact_counts,
            "sidechain_contact_weight_list": contact_weights,
            "enable_plip_fill": bool(args.enable_plip_fill),
            "plip_fill_max_res": int(args.plip_fill_max_res),
            "plip_fill_top_try_per_atom": int(args.plip_fill_top_try_per_atom),
            "n_runs": len(runs),
        },
        "aggregate": {
            "plip_all_satisfied_count": int(plip_sat_count),
            "plip_all_satisfied_rate": float(plip_sat_rate),
            "feasible_count": int(sum(int(bool(r.get("feasible"))) for r in runs)),
        },
        "recommended_run_id": recommended.get("run_id") if recommended else None,
        "recommended_feasible": bool(recommended.get("feasible")) if recommended else False,
        "runs": runs_sorted,
    }

    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        f"# MTX Pocket Contact Sweep ({args.tag})",
        "",
        f"- outdir: `{outdir}`",
        f"- recommended: `{summary['recommended_run_id']}` (feasible={summary['recommended_feasible']})",
        f"- plip all-satisfied rate: {summary['aggregate']['plip_all_satisfied_rate']:.3f} ({summary['aggregate']['plip_all_satisfied_count']}/{len(runs)})",
        "",
        "| run_id | feasible | top | cnt | weight | n_candidates | n_selected | simple_sat | plip_sat | clash_ok | internal_ok | total_s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in runs_sorted:
        m = r.get("metrics", {})
        t = r.get("timing_s", {})
        p = r.get("params", {})
        md_lines.append(
            "| "
            + f"{r.get('run_id')} | {int(bool(r.get('feasible')))} | {p.get('top_per_site')} | {p.get('min_sidechain_contact_count')} | "
            + f"{float(p.get('sidechain_contact_weight', 0.0)):.2f} | {m.get('n_candidates', 0)} | {m.get('n_selected')} | "
            + f"{int(bool(m.get('all_satisfied_simple')))} | {int(bool(m.get('all_satisfied_plip')))} | "
            + f"{int(bool(m.get('ligand_clash_ok')))} | {int(bool(m.get('internal_clash_ok')))} | "
            + f"{float(t.get('total_s', 0.0)):.2f} |"
        )
    (outdir / "summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
