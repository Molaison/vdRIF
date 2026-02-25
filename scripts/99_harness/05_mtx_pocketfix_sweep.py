#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class SweepConfig:
    name: str
    top_per_site_per_atom: int
    min_sidechain_centroid_dot: float
    min_pocket_contact_count: int
    pocket_contact_dist: float
    pocket_contact_score_weight: float
    solver_pocket_contact_score_weight: float


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _ensure_file(path: Path, hint: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}\nHint: {hint}")


def _default_configs() -> list[SweepConfig]:
    return [
        SweepConfig(
            name="baseline_legacy",
            top_per_site_per_atom=0,
            min_sidechain_centroid_dot=0.0,
            min_pocket_contact_count=0,
            pocket_contact_dist=4.5,
            pocket_contact_score_weight=0.0,
            solver_pocket_contact_score_weight=0.0,
        ),
        SweepConfig(
            name="pocketfix_default",
            top_per_site_per_atom=32,
            min_sidechain_centroid_dot=0.05,
            min_pocket_contact_count=1,
            pocket_contact_dist=4.5,
            pocket_contact_score_weight=0.15,
            solver_pocket_contact_score_weight=0.15,
        ),
        SweepConfig(
            name="pocketfix_weighted",
            top_per_site_per_atom=32,
            min_sidechain_centroid_dot=0.05,
            min_pocket_contact_count=1,
            pocket_contact_dist=4.5,
            pocket_contact_score_weight=0.25,
            solver_pocket_contact_score_weight=0.25,
        ),
        SweepConfig(
            name="pocketfix_strict",
            top_per_site_per_atom=64,
            min_sidechain_centroid_dot=0.10,
            min_pocket_contact_count=2,
            pocket_contact_dist=4.5,
            pocket_contact_score_weight=0.25,
            solver_pocket_contact_score_weight=0.25,
        ),
    ]


def _config_map(configs: list[SweepConfig]) -> dict[str, SweepConfig]:
    return {c.name: c for c in configs}


def _to_md(summary: dict[str, Any]) -> str:
    lines = []
    lines.append(f"# MTX Pocketfix Sweep ({summary['tag']})")
    lines.append("")
    lines.append("## Inputs")
    lines.append(f"- repo_root: `{summary['inputs']['repo_root']}`")
    lines.append(f"- python: `{summary['inputs']['python']}`")
    lines.append(f"- top_per_site: `{summary['inputs']['top_per_site']}`")
    lines.append(f"- solver: `{summary['inputs']['solver']}`")
    lines.append("")
    best = summary.get("best")
    if best is not None:
        lines.append(f"## Best Config: `{best['name']}`")
        lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append(
        "| config | status | cand_n | selected_n | polar | ligand_clash | internal_clash | selected_pocket_mean | internal_worst |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in summary["results"]:
        lines.append(
            "| "
            + " | ".join(
                [
                    r["name"],
                    r.get("status", "error"),
                    str(r.get("n_candidates", "")),
                    str(r.get("n_selected", "")),
                    str(r.get("polar_all_satisfied", "")),
                    str(r.get("ligand_clash_ok", "")),
                    str(r.get("internal_clash_ok", "")),
                    f"{r.get('selected_pocket_score_mean', 0.0):.3f}" if isinstance(r.get("selected_pocket_score_mean"), (int, float)) else "",
                    f"{r.get('internal_worst_overlap', 0.0):.3f}" if isinstance(r.get("internal_worst_overlap"), (int, float)) else "",
                ]
            )
            + " |"
        )
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Sweep MTX pocket-quality parameters and compare motif quality.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--solver", type=str, default="cp_sat", choices=["cp_sat", "greedy"])
    ap.add_argument("--time-limit-s", type=float, default=30.0)
    ap.add_argument("--python", type=str, default=sys.executable, help="Python executable to run pipeline scripts.")
    ap.add_argument("--ligand-pdb", type=Path, default=Path("inputs/01_cgmap/MTX.pdb"))
    ap.add_argument("--polar-sites", type=Path, default=Path("outputs/02_polar_sites/MTX_polar_sites.json"))
    ap.add_argument("--site-frames", type=Path, default=Path("outputs/02_polar_sites/MTX_site_frames.json"))
    ap.add_argument(
        "--configs",
        type=str,
        default="baseline_legacy,pocketfix_default,pocketfix_weighted,pocketfix_strict",
        help="Comma-separated preset names to run.",
    )
    ap.add_argument(
        "--vdxform-dir",
        type=Path,
        default=Path("processed/03_vdxform"),
        help="Directory containing per-cg vdxform_<cg>.npz files.",
    )
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"mtx_pocketfix_sweep_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand = (root / args.ligand_pdb).resolve() if not args.ligand_pdb.is_absolute() else args.ligand_pdb.resolve()
    polar = (root / args.polar_sites).resolve() if not args.polar_sites.is_absolute() else args.polar_sites.resolve()
    frames = (root / args.site_frames).resolve() if not args.site_frames.is_absolute() else args.site_frames.resolve()
    vdx = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()

    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        _ensure_file(
            vdx / cg / f"vdxform_{cg}.npz",
            hint="Generate vdxform first, e.g. bash scripts/03_vdxform/03_run_mtx_needed_cgs_debug.sh",
        )
    _ensure_file(ligand, hint="Place MTX at inputs/01_cgmap/MTX.pdb")
    _ensure_file(polar, hint="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")
    _ensure_file(frames, hint="bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh")

    candidate_py = root / "scripts/04_candidates/01_generate_candidates.py"
    solver_py = root / "scripts/05_solver/01_solve_motif.py"
    validate_sat_py = root / "scripts/05_solver/03_validate_motif_polar_satisfaction.py"
    validate_clash_py = root / "scripts/05_solver/04_validate_motif_clashes.py"
    validate_internal_clash_py = root / "scripts/05_solver/05_validate_motif_internal_clashes.py"

    all_cfg = _config_map(_default_configs())
    chosen_names = [x.strip() for x in str(args.configs).split(",") if x.strip()]
    configs = [all_cfg[n] for n in chosen_names]

    results: list[dict[str, Any]] = []
    for cfg in configs:
        run_dir = outdir / cfg.name
        run_dir.mkdir(parents=True, exist_ok=True)
        cand_prefix = run_dir / "MTX_candidates"
        motif_json = run_dir / "MTX_motif.json"
        motif_pdb = run_dir / "MTX_motif.pdb"
        val_json = run_dir / "MTX_motif_validation.json"
        clash_json = run_dir / "MTX_motif_clash_validation.json"
        internal_json = run_dir / "MTX_motif_internal_clash_validation.json"

        t0 = time.time()
        row: dict[str, Any] = {"name": cfg.name, "params": asdict(cfg)}
        try:
            _run(
                [
                    args.python,
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
                    str(int(args.top_per_site)),
                    "--top-per-site-per-atom",
                    str(cfg.top_per_site_per_atom),
                    "--clash-tol",
                    "0.5",
                    "--exclude-aa3",
                    "PRO,CYS",
                    "--require-sidechain-facing",
                    "--min-sidechain-facing-dot",
                    "0.2",
                    "--min-sidechain-centroid-dot",
                    str(cfg.min_sidechain_centroid_dot),
                    "--min-pocket-contact-count",
                    str(cfg.min_pocket_contact_count),
                    "--pocket-contact-dist",
                    str(cfg.pocket_contact_dist),
                    "--pocket-contact-score-weight",
                    str(cfg.pocket_contact_score_weight),
                    "--require-full-coverage",
                ]
            )
            t1 = time.time()

            _run(
                [
                    args.python,
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
                    str(float(args.time_limit_s)),
                    "--num-workers",
                    "1",
                    "--grid-size",
                    "4.0",
                    "--ca-prefilter",
                    "8.0",
                    "--clash-tol",
                    "0.5",
                    "--pocket-contact-score-weight",
                    str(cfg.solver_pocket_contact_score_weight),
                ]
            )
            t2 = time.time()

            _run([args.python, str(validate_sat_py), "--polar-sites", str(polar), "--motif-pdb", str(motif_pdb), "-o", str(val_json)])
            _run([args.python, str(validate_clash_py), "--motif-pdb", str(motif_pdb), "--ligand-resname", "MTX", "-o", str(clash_json)])
            _run(
                [
                    args.python,
                    str(validate_internal_clash_py),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    "MTX",
                    "-o",
                    str(internal_json),
                ]
            )

            cand_meta = _load_json(cand_prefix.with_suffix(".json"))
            solve = _load_json(motif_json)
            val = _load_json(val_json)
            clash = _load_json(clash_json)
            internal = _load_json(internal_json)

            selected = solve.get("selected", [])
            p_scores = [float(x.get("pocket_contact_score", 0.0)) for x in selected]
            p_counts = [float(x.get("pocket_contact_count", 0.0)) for x in selected]

            row.update(
                {
                    "status": "ok",
                    "candidate_s": t1 - t0,
                    "solve_s": t2 - t1,
                    "n_candidates": int(cand_meta.get("n_candidates", 0)),
                    "empty_candidates": bool(cand_meta.get("empty_candidates", False)),
                    "n_selected": int(solve.get("n_selected", 0)),
                    "coverage_complete": bool(solve.get("coverage_complete", False)),
                    "polar_all_satisfied": bool(val.get("all_satisfied", False)),
                    "ligand_clash_ok": bool(clash.get("ok", False)),
                    "ligand_worst_overlap": float(clash.get("worst_overlap", 0.0)),
                    "internal_clash_ok": bool(internal.get("ok", False)),
                    "internal_worst_overlap": float(internal.get("worst_overlap", 0.0)),
                    "selected_pocket_score_mean": float(np.mean(p_scores)) if p_scores else 0.0,
                    "selected_pocket_count_mean": float(np.mean(p_counts)) if p_counts else 0.0,
                    "run_dir": str(run_dir),
                }
            )
        except Exception as e:
            row.update({"status": "error", "error": repr(e), "run_dir": str(run_dir)})
        results.append(row)

    valid = [
        r
        for r in results
        if r.get("status") == "ok"
        and r.get("coverage_complete")
        and r.get("polar_all_satisfied")
        and r.get("ligand_clash_ok")
        and r.get("internal_clash_ok")
        and not r.get("empty_candidates")
    ]
    valid_sorted = sorted(
        valid,
        key=lambda r: (
            -float(r.get("selected_pocket_score_mean", 0.0)),
            -float(r.get("selected_pocket_count_mean", 0.0)),
            float(r.get("internal_worst_overlap", 999.0)),
            float(r.get("ligand_worst_overlap", 999.0)),
            int(r.get("n_selected", 999)),
        ),
    )
    best = valid_sorted[0] if valid_sorted else None

    summary = {
        "tag": str(args.tag),
        "inputs": {
            "repo_root": str(root),
            "python": str(args.python),
            "ligand_pdb": str(ligand),
            "polar_sites": str(polar),
            "site_frames": str(frames),
            "top_per_site": int(args.top_per_site),
            "solver": str(args.solver),
            "time_limit_s": float(args.time_limit_s),
            "vdxform_dir": str(vdx),
            "configs": chosen_names,
        },
        "results": results,
        "best": best,
    }
    report_json = outdir / "report.json"
    report_md = outdir / "report.md"
    _write_json(report_json, summary)
    report_md.write_text(_to_md(summary), encoding="utf-8")

    print(f"[done] report: {report_json}")
    if best is None:
        print("[done] no valid config found")
    else:
        print(f"[done] best: {best['name']}")


if __name__ == "__main__":
    main()
