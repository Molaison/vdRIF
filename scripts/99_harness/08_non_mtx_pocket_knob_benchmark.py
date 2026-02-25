#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class LigandAssets:
    code: str
    ligand_pdb: Path
    polar_sites: Path
    site_frames: Path


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd))


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_codes(csv: str) -> list[str]:
    out = [x.strip().upper() for x in str(csv).split(",") if x.strip()]
    if not out:
        raise ValueError("Expected non-empty ligand code list.")
    for c in out:
        if len(c) > 3:
            raise ValueError(f"Ligand code must be <=3 chars for PDB residue name, got: {c}")
    return out


def _resolve_ligand_assets(prepared_root: Path, code: str) -> LigandAssets:
    cdir = prepared_root / code
    ligand_pdb = cdir / f"{code}.pdb"
    polar_sites = cdir / f"{code}_polar_new.json"
    site_frames_candidates = [
        cdir / "new_pipeline" / f"{code}_site_frames.json",
        cdir / f"{code}_site_frames.json",
    ]
    site_frames = next((p for p in site_frames_candidates if p.exists()), None)
    if site_frames is None:
        raise FileNotFoundError(
            f"Missing site-frames for {code}. Tried: " + ", ".join(str(p) for p in site_frames_candidates)
        )
    for need in [ligand_pdb, polar_sites]:
        if not need.exists():
            raise FileNotFoundError(f"Missing prepared ligand asset for {code}: {need}")
    return LigandAssets(code=code, ligand_pdb=ligand_pdb, polar_sites=polar_sites, site_frames=site_frames)


def _run_single_ligand(
    *,
    repo_root: Path,
    python_bin: Path,
    sweep_py: Path,
    vdxform_dir: Path,
    lig: LigandAssets,
    tag: str,
    solver: str,
    top_per_site: int,
    top_per_site_per_atom: int,
    time_limit_s: float,
    min_res: int,
    max_res: int,
    ca_prefilter: float,
    angle_list: str,
    weight_list: str,
    min_sc_list: str,
    max_runs: int,
) -> dict[str, Any]:
    sweep_tag = f"{tag}_{lig.code.lower()}"
    _run(
        [
            str(python_bin),
            str(sweep_py),
            "--repo-root",
            str(repo_root),
            "--tag",
            sweep_tag,
            "--python-bin",
            str(python_bin),
            "--ligand-tag",
            str(lig.code),
            "--ligand-resname",
            str(lig.code),
            "--ligand-pdb",
            str(lig.ligand_pdb),
            "--polar-sites",
            str(lig.polar_sites),
            "--site-frames",
            str(lig.site_frames),
            "--vdxform-dir",
            str(vdxform_dir),
            "--solver",
            str(solver),
            "--top-per-site",
            str(int(top_per_site)),
            "--top-per-site-per-atom",
            str(int(top_per_site_per_atom)),
            "--time-limit-s",
            str(float(time_limit_s)),
            "--min-res",
            str(int(min_res)),
            "--max-res",
            str(int(max_res)),
            "--ca-prefilter",
            str(float(ca_prefilter)),
            "--min-lig-donor-angle-list",
            str(angle_list),
            "--pocket-contact-weight-list",
            str(weight_list),
            "--min-pocket-sidechain-contacts-list",
            str(min_sc_list),
            "--max-runs",
            str(int(max_runs)),
        ],
        cwd=repo_root,
    )
    summary_path = repo_root / "processed/99_harness" / f"mtx_pocket_knob_sweep_{sweep_tag}" / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Expected sweep summary not found: {summary_path}")
    summary = _read_json(summary_path)
    runs = list(summary.get("runs", []))
    feasible_runs = [r for r in runs if bool(r.get("feasible", False))]
    rec_id = summary.get("recommended_run_id")
    rec = next((r for r in runs if str(r.get("run_id")) == str(rec_id)), None)
    return {
        "status": "ok",
        "summary_path": str(summary_path),
        "n_runs": len(runs),
        "n_feasible_runs": len(feasible_runs),
        "recommended_run_id": rec_id,
        "recommended_feasible": bool(summary.get("recommended_feasible", False)),
        "recommended": rec,
        "summary": summary,
    }


def _config_key(run: dict[str, Any]) -> tuple[float, float, int]:
    p = run.get("params") or {}
    return (
        float(p.get("min_lig_donor_angle_deg", 0.0)),
        float(p.get("pocket_contact_weight", 0.0)),
        int(p.get("min_pocket_sidechain_contacts", 0)),
    )


def _safe_mean(vals: list[float], default: float) -> float:
    return float(statistics.mean(vals)) if vals else float(default)


def _build_config_summary(per_ligand: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    ok_rows = [r for r in per_ligand if r.get("status") == "ok"]
    n_ligands = len(ok_rows)
    all_keys: set[tuple[float, float, int]] = set()
    for row in ok_rows:
        for run in ((row.get("sweep") or {}).get("runs") or []):
            all_keys.add(_config_key(run))

    rows: list[dict[str, Any]] = []
    for key in sorted(all_keys):
        pass_count = 0
        qualities: list[float] = []
        contact_means: list[float] = []
        total_s_list: list[float] = []
        for row in ok_rows:
            runs = ((row.get("sweep") or {}).get("runs") or [])
            run = next((x for x in runs if _config_key(x) == key), None)
            if run is None:
                continue
            if bool(run.get("feasible", False)):
                pass_count += 1
                qualities.append(float(run.get("quality_score", 0.0)))
                contact_means.append(float(((run.get("metrics") or {}).get("pocket_contact_stats") or {}).get("mean", 0.0)))
                total_s_list.append(float((run.get("timing_s") or {}).get("total_s", 0.0)))

        row = {
            "min_lig_donor_angle_deg": float(key[0]),
            "pocket_contact_weight": float(key[1]),
            "min_pocket_sidechain_contacts": int(key[2]),
            "pass_count": int(pass_count),
            "n_ligands": int(n_ligands),
            "pass_all": bool(n_ligands > 0 and pass_count == n_ligands),
            "avg_quality_score_feasible": _safe_mean(qualities, default=0.0),
            "avg_contact_mean_feasible": _safe_mean(contact_means, default=0.0),
            "avg_total_s_feasible": _safe_mean(total_s_list, default=0.0),
        }
        rows.append(row)

    rows_sorted = sorted(
        rows,
        key=lambda r: (
            int(r["pass_count"]),
            float(r["avg_quality_score_feasible"]),
            float(r["avg_contact_mean_feasible"]),
            -float(r["avg_total_s_feasible"]),
        ),
        reverse=True,
    )
    best = rows_sorted[0] if rows_sorted else None
    return rows_sorted, best


def _write_md(path: Path, summary: dict[str, Any]) -> None:
    lines: list[str] = [
        f"# Non-MTX Pocket Knob Benchmark ({summary['tag']})",
        "",
        "## Inputs",
        f"- repo_root: `{summary['inputs']['repo_root']}`",
        f"- python_bin: `{summary['inputs']['python_bin']}`",
        f"- prepared_root: `{summary['inputs']['prepared_root']}`",
        f"- vdxform_dir: `{summary['inputs']['vdxform_dir']}`",
        f"- ligands: `{','.join(summary['inputs']['ligands'])}`",
        f"- sweep knobs: angle=`{summary['inputs']['angle_list']}`, weight=`{summary['inputs']['weight_list']}`, min_sc=`{summary['inputs']['min_sc_list']}`",
        "",
    ]

    best = summary.get("best_config")
    if best:
        lines.append(
            "## Best Shared Knob Set: "
            + f"`angle={best['min_lig_donor_angle_deg']}, weight={best['pocket_contact_weight']}, min_sc={best['min_pocket_sidechain_contacts']}`"
        )
        lines.append("")

    lines.extend(
        [
            "## Per Ligand",
            "",
            "| ligand | status | n_runs | n_feasible | recommended_run | recommended_feasible | best_quality | best_contact_mean | best_total_s |",
            "|---|---|---:|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary.get("per_ligand", []):
        rec = row.get("recommended") or {}
        m = rec.get("metrics") or {}
        contact_mean = float((m.get("pocket_contact_stats") or {}).get("mean", 0.0))
        total_s = float((rec.get("timing_s") or {}).get("total_s", 0.0))
        lines.append(
            "| "
            + f"{row.get('code')} | {row.get('status')} | {row.get('n_runs', 0)} | {row.get('n_feasible_runs', 0)} | "
            + f"{row.get('recommended_run_id', '-')} | {int(bool(row.get('recommended_feasible', False)))} | "
            + f"{float(rec.get('quality_score', 0.0)):.2f} | {contact_mean:.3f} | {total_s:.2f} |"
        )

    lines.extend(
        [
            "",
            "## Shared Knob Aggregate",
            "",
            "| angle | weight | min_sc | pass_count | n_ligands | pass_all | avg_quality(feasible) | avg_contact(feasible) | avg_total_s(feasible) |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.get("config_summary", []):
        lines.append(
            "| "
            + f"{row.get('min_lig_donor_angle_deg', 0.0):.1f} | {row.get('pocket_contact_weight', 0.0):.2f} | "
            + f"{row.get('min_pocket_sidechain_contacts', 0)} | {row.get('pass_count', 0)} | {row.get('n_ligands', 0)} | "
            + f"{int(bool(row.get('pass_all', False)))} | {float(row.get('avg_quality_score_feasible', 0.0)):.2f} | "
            + f"{float(row.get('avg_contact_mean_feasible', 0.0)):.3f} | {float(row.get('avg_total_s_feasible', 0.0)):.2f} |"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Benchmark pocket-quality knob transferability on prepared non-MTX ligands.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=time.strftime("%Y%m%d-%H%M%S"))
    ap.add_argument("--python-bin", type=Path, required=True)
    ap.add_argument("--prepared-root", type=Path, default=Path("processed/99_harness/non_mtx_cation_validate_20260226-cation-nonmtx-v2/ligands"))
    ap.add_argument("--ligands", type=str, default="GDM,BZA,BTM")
    ap.add_argument("--vdxform-dir", type=Path, required=True)
    ap.add_argument("--solver", type=str, default="greedy", choices=["cp_sat", "greedy"])
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--top-per-site-per-atom", type=int, default=50)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--min-res", type=int, default=3)
    ap.add_argument("--max-res", type=int, default=12)
    ap.add_argument("--ca-prefilter", type=float, default=12.0)
    ap.add_argument("--angle-list", type=str, default="60")
    ap.add_argument("--weight-list", type=str, default="0.0,0.2")
    ap.add_argument("--min-sc-list", type=str, default="0,1")
    ap.add_argument("--max-runs", type=int, default=0)
    args = ap.parse_args()

    root = args.repo_root.resolve()
    python_bin = args.python_bin if args.python_bin.is_absolute() else (root / args.python_bin)
    prepared_root = args.prepared_root if args.prepared_root.is_absolute() else (root / args.prepared_root)
    vdxform_dir = args.vdxform_dir if args.vdxform_dir.is_absolute() else (root / args.vdxform_dir)
    sweep_py = root / "scripts/99_harness/05_mtx_pocket_knob_sweep.py"
    if not sweep_py.exists():
        raise FileNotFoundError(f"Missing sweep script: {sweep_py}")
    if not python_bin.exists():
        raise FileNotFoundError(f"Python binary not found: {python_bin}")
    if not prepared_root.exists():
        raise FileNotFoundError(f"Prepared ligand root not found: {prepared_root}")

    outdir = root / "processed/99_harness" / f"non_mtx_pocket_knob_benchmark_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand_codes = _parse_codes(args.ligands)
    per_ligand: list[dict[str, Any]] = []
    for code in ligand_codes:
        row: dict[str, Any] = {"code": code}
        try:
            assets = _resolve_ligand_assets(prepared_root, code)
            row["assets"] = {
                "ligand_pdb": str(assets.ligand_pdb),
                "polar_sites": str(assets.polar_sites),
                "site_frames": str(assets.site_frames),
            }
            result = _run_single_ligand(
                repo_root=root,
                python_bin=python_bin,
                sweep_py=sweep_py,
                vdxform_dir=vdxform_dir,
                lig=assets,
                tag=str(args.tag),
                solver=str(args.solver),
                top_per_site=int(args.top_per_site),
                top_per_site_per_atom=int(args.top_per_site_per_atom),
                time_limit_s=float(args.time_limit_s),
                min_res=int(args.min_res),
                max_res=int(args.max_res),
                ca_prefilter=float(args.ca_prefilter),
                angle_list=str(args.angle_list),
                weight_list=str(args.weight_list),
                min_sc_list=str(args.min_sc_list),
                max_runs=int(args.max_runs),
            )
            row.update({k: v for k, v in result.items() if k != "summary"})
            row["sweep"] = result["summary"]
        except Exception as e:
            row["status"] = "error"
            row["error"] = repr(e)
        per_ligand.append(row)

    config_summary, best_cfg = _build_config_summary(per_ligand)
    summary = {
        "tag": str(args.tag),
        "inputs": {
            "repo_root": str(root),
            "python_bin": str(python_bin),
            "prepared_root": str(prepared_root),
            "vdxform_dir": str(vdxform_dir),
            "ligands": ligand_codes,
            "solver": str(args.solver),
            "top_per_site": int(args.top_per_site),
            "top_per_site_per_atom": int(args.top_per_site_per_atom),
            "time_limit_s": float(args.time_limit_s),
            "min_res": int(args.min_res),
            "max_res": int(args.max_res),
            "ca_prefilter": float(args.ca_prefilter),
            "angle_list": str(args.angle_list),
            "weight_list": str(args.weight_list),
            "min_sc_list": str(args.min_sc_list),
            "max_runs": int(args.max_runs),
        },
        "per_ligand": per_ligand,
        "config_summary": config_summary,
        "best_config": best_cfg,
    }

    out_json = outdir / "summary.json"
    out_md = outdir / "summary.md"
    _write_json(out_json, summary)
    _write_md(out_md, summary)
    print(str(out_json))
    print(str(out_md))


if __name__ == "__main__":
    main()
