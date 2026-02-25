#!/usr/bin/env python
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd))


def _run_allow_failure(cmd: list[str], cwd: Path) -> int:
    return int(subprocess.run(cmd, check=False, cwd=str(cwd)).returncode)


def _parse_int_list(s: str) -> list[int]:
    vals = [x.strip() for x in str(s).split(",") if x.strip()]
    if not vals:
        raise ValueError(f"Expected non-empty int list, got: {s!r}")
    out: list[int] = []
    for v in vals:
        out.append(int(v))
    return out


def _bit_count_u16(x: int) -> int:
    return int(bin(int(x) & 0xFFFF).count("1"))


def _run_dir_name(
    run_id: int,
    target_res: int,
    target_res_penalty: int,
    site_diversity_reward: int,
    min_unique_sites: int,
    max_per_site: int,
) -> str:
    return (
        f"run_{run_id:03d}_tr{target_res}_tp{target_res_penalty}"
        f"_sd{site_diversity_reward}_mus{min_unique_sites}_mps{max_per_site}"
    )


def _evenly_spaced_indices(n_items: int, n_pick: int) -> list[int]:
    if n_items <= 0:
        return []
    if n_pick <= 1:
        return [0]
    if n_pick >= n_items:
        return list(range(n_items))

    step = float(n_items - 1) / float(n_pick - 1)
    out: list[int] = []
    last = -1
    for i in range(n_pick):
        idx = int(round(float(i) * step))
        if idx <= last:
            idx = last + 1
        if idx >= n_items:
            idx = n_items - 1
        out.append(idx)
        last = idx
    return out


def _to_md(summary: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# MTX balanced solver sweep")
    lines.append("")
    lines.append(f"- tag: `{summary['tag']}`")
    lines.append(f"- total runs: `{summary['sweep']['n_runs']}`")
    rec = summary.get("recommended")
    if rec is None:
        lines.append("- recommended: none (no successful run)")
    else:
        lines.append(f"- recommended run_id: `{rec['run_id']}`")
        lines.append(
            "- recommended metrics: "
            f"coverage_complete={rec['coverage_complete']}, "
            f"n_unique_sites_selected={rec['n_unique_sites_selected']}, "
            f"n_selected={rec['n_selected']}, "
            f"target_dev={rec['target_res_dev']}"
        )
    lines.append("")
    lines.append("## Top runs")
    lines.append("")
    lines.append(
        "| rank | run_id | feasible | internal_ok | worst_overlap | coverage_bits | unique_sites | n_selected | target_dev | avg_score | params |"
    )
    lines.append("|---:|---:|:---:|:---:|---:|---:|---:|---:|---:|---:|---|")
    for rank, run in enumerate(summary["runs"][:10], start=1):
        params = run["params"]
        worst_overlap = run.get("internal_clash_worst_overlap")
        worst_overlap_str = "-" if worst_overlap is None else f"{float(worst_overlap):.3f}"
        lines.append(
            f"| {rank} | {run['run_id']} | {str(run['feasible']).lower()} | "
            f"{str(bool(run.get('internal_clash_ok', False))).lower()} | {worst_overlap_str} | "
            f"{run['coverage_bits']} | {run['n_unique_sites_selected']} | {run['n_selected']} | "
            f"{run['target_res_dev']} | {run['avg_selected_score']:.3f} | "
            f"`tr={params['target_res']},tp={params['target_res_penalty']},"
            f"sd={params['site_diversity_reward']},mus={params['min_unique_sites']},"
            f"mps={params['max_per_site']}` |"
        )
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Sweep balanced solver parameters on a fixed MTX candidates set.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--python", type=str, default=sys.executable, help="Python executable for solver script.")

    ap.add_argument("--candidates-npz", type=Path, default=Path("processed/04_candidates/MTX_candidates.npz"))
    ap.add_argument("--candidates-meta", type=Path, default=Path("processed/04_candidates/MTX_candidates.json"))
    ap.add_argument("--ligand-pdb", type=Path, default=Path("inputs/01_cgmap/MTX.pdb"))

    # Balanced objective is implemented for CP-SAT only in the solver.
    ap.add_argument("--solver", type=str, default="cp_sat", choices=["cp_sat"])
    ap.add_argument("--min-res", type=int, default=8)
    ap.add_argument("--max-res", type=int, default=15)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--num-workers", type=int, default=1)
    ap.add_argument("--grid-size", type=float, default=4.0)
    ap.add_argument("--ca-prefilter", type=float, default=12.0)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument("--ligand-resname", type=str, default="MTX")
    ap.add_argument("--internal-clash-fail-overlap", type=float, default=0.5)
    ap.add_argument("--internal-clash-tol", type=float, default=0.5)
    ap.add_argument(
        "--no-enforce-internal-clash",
        action="store_false",
        dest="enforce_internal_clash",
        help="Do not require internal clash validator to pass for feasibility.",
    )
    ap.set_defaults(enforce_internal_clash=True)

    ap.add_argument("--target-res-list", type=str, default="10,11,12")
    ap.add_argument("--target-res-penalty-list", type=str, default="1200,2000,3000")
    ap.add_argument("--site-diversity-reward-list", type=str, default="800,1200,2000")
    ap.add_argument("--min-unique-sites-list", type=str, default="0,6")
    ap.add_argument("--max-per-site-list", type=str, default="0,2")
    ap.add_argument("--max-runs", type=int, default=24)
    args = ap.parse_args()

    root = args.repo_root.resolve()
    candidates_npz = (root / args.candidates_npz).resolve()
    candidates_meta = (root / args.candidates_meta).resolve()
    ligand_pdb = (root / args.ligand_pdb).resolve()

    if not candidates_npz.exists():
        raise FileNotFoundError(f"Missing candidates npz: {candidates_npz}")
    if not candidates_meta.exists():
        raise FileNotFoundError(f"Missing candidates meta: {candidates_meta}")
    if not ligand_pdb.exists():
        raise FileNotFoundError(f"Missing ligand PDB: {ligand_pdb}")

    target_res_list = _parse_int_list(args.target_res_list)
    target_res_penalty_list = _parse_int_list(args.target_res_penalty_list)
    site_diversity_reward_list = _parse_int_list(args.site_diversity_reward_list)
    min_unique_sites_list = _parse_int_list(args.min_unique_sites_list)
    max_per_site_list = _parse_int_list(args.max_per_site_list)

    available_site_count: int | None = None
    with np.load(candidates_npz, allow_pickle=False) as z:
        if "site_index_u16" in z.files:
            available_site_count = int(np.unique(z["site_index_u16"].astype(np.uint16)).size)

    all_combos = list(
        itertools.product(
            target_res_list,
            target_res_penalty_list,
            site_diversity_reward_list,
            min_unique_sites_list,
            max_per_site_list,
        )
    )
    max_runs = max(int(args.max_runs), 1)
    combo_indices = _evenly_spaced_indices(len(all_combos), max_runs)
    combos = [all_combos[i] for i in combo_indices]

    outdir = root / "processed/99_harness" / f"mtx_balanced_solver_sweep_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)
    solver_py = root / "scripts/05_solver/01_solve_motif.py"
    internal_validator_py = root / "scripts/05_solver/05_validate_motif_internal_clashes.py"

    runs: list[dict[str, Any]] = []
    for i, (target_res, target_res_penalty, site_diversity_reward, min_unique_sites, max_per_site) in enumerate(
        combos, start=1
    ):
        run_name = _run_dir_name(
            run_id=i,
            target_res=target_res,
            target_res_penalty=target_res_penalty,
            site_diversity_reward=site_diversity_reward,
            min_unique_sites=min_unique_sites,
            max_per_site=max_per_site,
        )
        run_dir = outdir / run_name
        run_dir.mkdir(parents=True, exist_ok=True)
        out_json = run_dir / "motif.json"
        out_pdb = run_dir / "motif.pdb"
        internal_clash_json = run_dir / "internal_clash.json"

        params = {
            "target_res": int(target_res),
            "target_res_penalty": int(target_res_penalty),
            "site_diversity_reward": int(site_diversity_reward),
            "min_unique_sites": int(min_unique_sites),
            "max_per_site": int(max_per_site),
        }

        cmd = [
            str(args.python),
            str(solver_py),
            "--candidates-npz",
            str(candidates_npz),
            "--candidates-meta",
            str(candidates_meta),
            "--ligand-pdb",
            str(ligand_pdb),
            "--out-json",
            str(out_json),
            "--out-pdb",
            str(out_pdb),
            "--solver",
            str(args.solver),
            "--min-res",
            str(int(args.min_res)),
            "--max-res",
            str(int(args.max_res)),
            "--time-limit-s",
            str(float(args.time_limit_s)),
            "--num-workers",
            str(int(args.num_workers)),
            "--grid-size",
            str(float(args.grid_size)),
            "--ca-prefilter",
            str(float(args.ca_prefilter)),
            "--clash-tol",
            str(float(args.clash_tol)),
            "--objective-mode",
            "balanced",
            "--target-res",
            str(int(target_res)),
            "--target-res-penalty",
            str(int(target_res_penalty)),
            "--site-diversity-reward",
            str(int(site_diversity_reward)),
            "--min-unique-sites",
            str(int(min_unique_sites)),
            "--max-per-site",
            str(int(max_per_site)),
        ]

        t0 = time.time()
        run_rec: dict[str, Any] = {
            "run_id": i,
            "run_name": run_name,
            "params": params,
            "motif_json": str(out_json),
            "motif_pdb": str(out_pdb),
            "internal_clash_json": str(internal_clash_json),
        }
        if available_site_count is not None and int(min_unique_sites) > int(available_site_count):
            run_rec.update(
                {
                    "ok": False,
                    "feasible": False,
                    "coverage_complete": False,
                    "coverage_bitmask": 0,
                    "coverage_bits": 0,
                    "n_selected": 0,
                    "n_unique_sites_selected": 0,
                    "selected_site_indices": [],
                    "target_res_dev": 10**9,
                    "avg_selected_score": -1e9,
                    "internal_clash_ok": False,
                    "internal_clash_ok_raw": False,
                    "internal_clash_worst_overlap": None,
                    "internal_clash_validation_returncode": None,
                    "error": (
                        f"Skipped: min_unique_sites={int(min_unique_sites)} exceeds "
                        f"available_site_count={int(available_site_count)}"
                    ),
                    "skipped_precheck": True,
                }
            )
            run_rec["elapsed_s"] = 0.0
            runs.append(run_rec)
            continue

        try:
            _run(cmd, cwd=root)
            report = _load_json(out_json)
            selected = report.get("selected", []) or []
            avg_score = 0.0
            if selected:
                avg_score = float(sum(float(x.get("score", 0.0)) for x in selected) / len(selected))
            n_selected = int(report.get("n_selected", 0))
            target_dev = abs(n_selected - int(target_res))
            coverage_bitmask = int(report.get("coverage_bitmask", 0))
            coverage_bits = _bit_count_u16(coverage_bitmask)
            coverage_complete = bool(report.get("coverage_complete", False))
            solver_status = str((report.get("solver") or {}).get("status", ""))
            internal_cmd = [
                str(args.python),
                str(internal_validator_py),
                "--motif-pdb",
                str(out_pdb),
                "-o",
                str(internal_clash_json),
                "--ligand-resname",
                str(args.ligand_resname),
                "--tol",
                str(float(args.internal_clash_tol)),
                "--fail-overlap",
                str(float(args.internal_clash_fail_overlap)),
            ]
            internal_rc = _run_allow_failure(internal_cmd, cwd=root)
            internal_ok_raw = False
            internal_worst_overlap: float | None = None
            if internal_clash_json.exists():
                internal_report = _load_json(internal_clash_json)
                internal_ok_raw = bool(internal_report.get("ok", False))
                worst = internal_report.get("worst_overlap")
                if worst is not None:
                    internal_worst_overlap = float(worst)
            internal_ok = bool(internal_ok_raw) if bool(args.enforce_internal_clash) else True
            feasible = coverage_complete and (solver_status in {"OPTIMAL", "FEASIBLE", "GREEDY"}) and internal_ok

            run_rec.update(
                {
                    "ok": True,
                    "solver_status": solver_status,
                    "coverage_complete": coverage_complete,
                    "coverage_bitmask": coverage_bitmask,
                    "coverage_bits": coverage_bits,
                    "n_selected": n_selected,
                    "n_unique_sites_selected": int(report.get("n_unique_sites_selected", 0)),
                    "selected_site_indices": report.get("selected_site_indices", []),
                    "target_res_dev": int(target_dev),
                    "avg_selected_score": float(avg_score),
                    "feasible": bool(feasible),
                    "internal_clash_ok": bool(internal_ok),
                    "internal_clash_ok_raw": bool(internal_ok_raw),
                    "internal_clash_worst_overlap": internal_worst_overlap,
                    "internal_clash_validation_returncode": int(internal_rc),
                    "objective_value": (report.get("solver") or {}).get("objective_value"),
                }
            )
        except Exception as e:
            run_rec.update(
                {
                    "ok": False,
                    "feasible": False,
                    "coverage_complete": False,
                    "coverage_bitmask": 0,
                    "coverage_bits": 0,
                    "n_selected": 0,
                    "n_unique_sites_selected": 0,
                    "selected_site_indices": [],
                    "target_res_dev": 10**9,
                    "avg_selected_score": -1e9,
                    "internal_clash_ok": False,
                    "internal_clash_ok_raw": False,
                    "internal_clash_worst_overlap": None,
                    "internal_clash_validation_returncode": None,
                    "error": str(e),
                    "skipped_precheck": False,
                }
            )
        finally:
            run_rec["elapsed_s"] = float(time.time() - t0)
            runs.append(run_rec)

    runs_sorted = sorted(
        runs,
        key=lambda r: (
            int(bool(r.get("feasible", False))),
            int(r.get("coverage_bits", 0)),
            int(r.get("n_unique_sites_selected", 0)),
            -float(r.get("internal_clash_worst_overlap", 1e9)),
            -int(r.get("target_res_dev", 10**9)),
            int(r.get("n_selected", 0)),
            float(r.get("avg_selected_score", -1e9)),
            -int(r.get("run_id", 10**9)),
        ),
        reverse=True,
    )
    best = next((r for r in runs_sorted if bool(r.get("feasible", False))), None)

    n_ok = int(sum(1 for r in runs if bool(r.get("ok", False))))
    n_feasible = int(sum(1 for r in runs if bool(r.get("feasible", False))))
    n_internal_ok = int(sum(1 for r in runs if bool(r.get("internal_clash_ok_raw", False))))
    n_skipped_precheck = int(sum(1 for r in runs if bool(r.get("skipped_precheck", False))))

    summary: dict[str, Any] = {
        "tag": str(args.tag),
        "outdir": str(outdir),
        "inputs": {
            "repo_root": str(root),
            "python": str(args.python),
            "candidates_npz": str(candidates_npz),
            "candidates_meta": str(candidates_meta),
            "ligand_pdb": str(ligand_pdb),
            "available_site_count": available_site_count,
        },
        "solver_args": {
            "solver": str(args.solver),
            "min_res": int(args.min_res),
            "max_res": int(args.max_res),
            "time_limit_s": float(args.time_limit_s),
            "num_workers": int(args.num_workers),
            "grid_size": float(args.grid_size),
            "ca_prefilter": float(args.ca_prefilter),
            "clash_tol": float(args.clash_tol),
            "ligand_resname": str(args.ligand_resname),
            "enforce_internal_clash": bool(args.enforce_internal_clash),
            "internal_clash_fail_overlap": float(args.internal_clash_fail_overlap),
            "internal_clash_tol": float(args.internal_clash_tol),
        },
        "sweep": {
            "n_all_combos": len(all_combos),
            "target_res_list": target_res_list,
            "target_res_penalty_list": target_res_penalty_list,
            "site_diversity_reward_list": site_diversity_reward_list,
            "min_unique_sites_list": min_unique_sites_list,
            "max_per_site_list": max_per_site_list,
            "n_runs": len(runs),
            "combo_indices": combo_indices,
            "n_ok": n_ok,
            "n_feasible": n_feasible,
            "n_internal_ok": n_internal_ok,
            "n_skipped_precheck": n_skipped_precheck,
        },
        "recommended": best,
        "runs": runs_sorted,
    }

    report_json = outdir / "report.json"
    report_md = outdir / "report.md"
    _write_json(report_json, summary)
    report_md.write_text(_to_md(summary), encoding="utf-8")
    print(f"[done] report: {report_json}")
    if best is None:
        print("[done] no feasible run found")
    else:
        print(
            "[done] recommended run: "
            f"id={best['run_id']} unique_sites={best.get('n_unique_sites_selected')} "
            f"n_selected={best.get('n_selected')} target_dev={best.get('target_res_dev')}"
        )


if __name__ == "__main__":
    main()
