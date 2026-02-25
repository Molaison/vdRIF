#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass(frozen=True)
class LigandSpec:
    code: str
    smiles: str


def _run(cmd: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd is not None else None)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_md(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        f"# Non-MTX Transfer Benchmark ({summary['tag']})",
        "",
        f"- repo_root: `{summary['inputs']['repo_root']}`",
        f"- vdxform_dir: `{summary['inputs']['vdxform_dir']}`",
        f"- ligands: `{','.join(summary['inputs']['ligands'])}`",
        f"- sweep_grid: `wcov={summary['inputs']['score_w_coverage_list']} wct={summary['inputs']['score_w_contact_list']} wsh={summary['inputs']['score_w_shell_list']} target={summary['inputs']['target_res_list']} min_cover={summary['inputs']['min_cover_per_polar_list']}`",
        "",
    ]
    best = summary.get("best_param_set")
    if isinstance(best, dict):
        lines.append(f"- best_param_set: `{json.dumps(best, sort_keys=True)}`")
        lines.append("")
    lines.extend(
        [
            "| ligand | status | feasible | recommended_run | plip | sweep_summary |",
            "|---|---|---:|---|---:|---|",
        ]
    )
    for row in summary.get("per_ligand", []):
        best_row = row.get("recommended_metrics") or {}
        plip = best_row.get("plip_all_satisfied")
        plip_txt = "-" if plip is None else str(int(bool(plip)))
        lines.append(
            "| "
            + f"{row.get('ligand')} | {row.get('status')} | {int(bool(row.get('recommended_feasible')))} | "
            + f"{row.get('recommended_run_id')} | {plip_txt} | {row.get('sweep_summary', '-')}"
            + " |"
        )
    lines.append("")
    lines.extend(
        [
            "| params | pass_count | n_ligands | pass_all | avg_quality | avg_total_s |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.get("param_summary", []):
        lines.append(
            "| "
            + f"{json.dumps(row.get('params', {}), sort_keys=True)} | {row.get('pass_count')} | {row.get('n_ligands')} | "
            + f"{int(bool(row.get('pass_all')))} | {float(row.get('avg_quality_score', 0.0)):.2f} | "
            + f"{float(row.get('avg_total_s', 0.0)):.2f} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_ligands(raw: str) -> list[LigandSpec]:
    out: list[LigandSpec] = []
    for token in [x.strip() for x in str(raw).split(",") if x.strip()]:
        if ":" not in token:
            raise ValueError(f"Invalid ligand token (expect CODE:SMILES): {token}")
        code, smiles = token.split(":", 1)
        code = code.strip().upper()
        if not code:
            raise ValueError(f"Empty ligand code in token: {token}")
        if len(code) > 3:
            raise ValueError(f"Ligand code must be <=3 chars for PDB residue name, got: {code}")
        if not smiles.strip():
            raise ValueError(f"Empty SMILES in token: {token}")
        out.append(LigandSpec(code=code, smiles=smiles.strip()))
    if not out:
        raise ValueError("No ligands provided.")
    return out


def _build_ligand_pdb(smiles: str, resname: str, out_pdb: Path) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise RuntimeError(f"RDKit embedding failed for SMILES: {smiles}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    for i, atom in enumerate(mol.GetAtoms(), start=1):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{i:>3}"[:4])
        info.SetResidueName(resname)
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(out_pdb))


def _safe_summary_name(code: str, tag: str) -> str:
    return f"{code.lower()}_pocket_quality_sweep_{tag}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Benchmark pocket-quality objective transferability on non-MTX ligands.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=time.strftime("%Y%m%d-%H%M%S"))
    ap.add_argument("--python-run", type=str, default=sys.executable, help="Python used for solver/sweep scripts.")
    ap.add_argument("--python-prep", type=str, default="", help="Python used for ligand prep + cgmap (default: --python-run).")
    ap.add_argument(
        "--ligands",
        type=str,
        default="PAB:Nc1ccc(cc1)C(=O)O,PBA:[NH3+]c1ccc(cc1)C(=O)[O-]",
        help="Comma-separated CODE:SMILES tokens for transfer benchmark.",
    )
    ap.add_argument(
        "--cgmap-script",
        type=Path,
        default=Path.home() / "practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py",
    )
    ap.add_argument(
        "--vdm-database",
        type=Path,
        default=Path.home() / "practice_arena/050_test_COMBS/Combs2024/database/database/vdMs",
    )
    ap.add_argument("--vdxform-dir", type=Path, default=Path("processed/03_vdxform_full"))
    ap.add_argument("--polar-backend", type=str, default="auto", choices=["auto", "openbabel", "rdkit"])
    ap.add_argument("--solver", type=str, default="cp_sat", choices=["cp_sat", "greedy"])
    ap.add_argument("--top-per-site", type=int, default=120)
    ap.add_argument("--top-per-site-per-atom", type=int, default=40)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--time-limit-s", type=float, default=120.0)
    ap.add_argument("--ca-prefilter", type=float, default=14.0)
    ap.add_argument("--acceptor-model", type=str, default="plip", choices=["legacy", "plip"])
    ap.add_argument("--score-w-coverage-list", type=str, default="0.03,0.05")
    ap.add_argument("--score-w-contact-list", type=str, default="0.1,0.2")
    ap.add_argument("--score-w-shell-list", type=str, default="0.2")
    ap.add_argument("--target-res-list", type=str, default="10,12")
    ap.add_argument("--min-cover-per-polar-list", type=str, default="1")
    ap.add_argument("--max-runs", type=int, default=8)
    ap.add_argument("--run-plip", action="store_true")
    ap.add_argument("--require-plip-success", action="store_true")
    ap.add_argument("--plip-bin", type=str, default="plip")
    ap.add_argument("--stop-on-error", action="store_true")
    args = ap.parse_args()

    root = args.repo_root.resolve()
    python_run = str(args.python_run)
    python_prep = str(args.python_prep).strip() or python_run
    vdxform_dir = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()
    ligands = _parse_ligands(args.ligands)

    sweep_py = root / "scripts/99_harness/05_mtx_pocket_quality_sweep.py"
    extract_py = root / "scripts/02_polar_sites/01_extract_polar_sites.py"
    cov_py = root / "scripts/02_polar_sites/02_check_polar_coverage_vs_cgmap.py"
    frames_py = root / "scripts/02_polar_sites/03_build_ligand_site_frames.py"
    cg_frame_defs = root / "configs/cg_frame_defs.json"

    for need in [sweep_py, extract_py, cov_py, frames_py, cg_frame_defs, args.cgmap_script, args.vdm_database, vdxform_dir]:
        if not Path(need).exists():
            raise FileNotFoundError(f"Missing required path: {need}")

    outdir = root / "processed/99_harness" / f"non_mtx_transfer_benchmark_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    per_ligand: list[dict[str, Any]] = []
    sweep_summaries_by_ligand: dict[str, dict[str, Any]] = {}
    for lig in ligands:
        row: dict[str, Any] = {"ligand": lig.code, "smiles": lig.smiles, "status": "error"}
        lig_dir = outdir / lig.code
        lig_dir.mkdir(parents=True, exist_ok=True)
        ligand_pdb = lig_dir / f"{lig.code}.pdb"
        cg_atommap = lig_dir / f"{lig.code}_cg_atommap.json"
        polar_sites = lig_dir / f"{lig.code}_polar_sites.json"
        polar_cov = lig_dir / f"{lig.code}_polar_coverage.json"
        site_frames = lig_dir / f"{lig.code}_site_frames.json"

        try:
            _build_ligand_pdb(lig.smiles, lig.code, ligand_pdb)
            _run(
                [
                    python_prep,
                    str(args.cgmap_script),
                    str(ligand_pdb),
                    "--format",
                    "json",
                    "--vdm-database",
                    str(args.vdm_database),
                    "-o",
                    str(cg_atommap),
                ],
                cwd=root,
            )

            backend = str(args.polar_backend)
            if backend == "auto":
                try:
                    _run([python_prep, str(extract_py), str(ligand_pdb), "-o", str(polar_sites), "--backend", "openbabel"], cwd=root)
                    backend = "openbabel"
                except Exception:
                    _run([python_run, str(extract_py), str(ligand_pdb), "-o", str(polar_sites), "--backend", "rdkit"], cwd=root)
                    backend = "rdkit"
            else:
                _run([python_run, str(extract_py), str(ligand_pdb), "-o", str(polar_sites), "--backend", backend], cwd=root)
            row["polar_backend"] = backend

            _run([python_run, str(cov_py), "--polar-sites", str(polar_sites), "--cg-atommap", str(cg_atommap), "-o", str(polar_cov)], cwd=root)
            _run(
                [
                    python_run,
                    str(frames_py),
                    str(ligand_pdb),
                    "--cg-atommap",
                    str(cg_atommap),
                    "--polar-sites",
                    str(polar_sites),
                    "--cg-frame-defs",
                    str(cg_frame_defs),
                    "--add-bb-cnh-for-uncovered-donors",
                    "--add-ccn-for-cations",
                    "-o",
                    str(site_frames),
                ],
                cwd=root,
            )

            frame_meta = _load_json(site_frames)
            row["n_frames"] = int(frame_meta.get("n_sites", 0))
            row["uncovered_polar_atoms"] = list(frame_meta.get("uncovered_polar_atoms") or [])
            if row["n_frames"] <= 0:
                raise RuntimeError("Prepared site_frames has no sites.")

            sweep_tag = f"{args.tag}_{lig.code.lower()}"
            sweep_cmd = [
                python_run,
                str(sweep_py),
                "--repo-root",
                str(root),
                "--tag",
                sweep_tag,
                "--python-bin",
                str(python_run),
                "--ligand-tag",
                lig.code,
                "--ligand-resname",
                lig.code,
                "--ligand-pdb",
                str(ligand_pdb),
                "--polar-sites",
                str(polar_sites),
                "--site-frames",
                str(site_frames),
                "--vdxform-dir",
                str(vdxform_dir),
                "--solver",
                str(args.solver),
                "--top-per-site",
                str(int(args.top_per_site)),
                "--top-per-site-per-atom",
                str(int(args.top_per_site_per_atom)),
                "--chunk-size",
                str(int(args.chunk_size)),
                "--time-limit-s",
                str(float(args.time_limit_s)),
                "--ca-prefilter",
                str(float(args.ca_prefilter)),
                "--acceptor-model",
                str(args.acceptor_model),
                "--score-w-coverage-list",
                str(args.score_w_coverage_list),
                "--score-w-contact-list",
                str(args.score_w_contact_list),
                "--score-w-shell-list",
                str(args.score_w_shell_list),
                "--target-res-list",
                str(args.target_res_list),
                "--min-cover-per-polar-list",
                str(args.min_cover_per_polar_list),
                "--max-runs",
                str(int(args.max_runs)),
            ]
            if bool(args.run_plip):
                sweep_cmd.extend(["--run-plip", "--plip-bin", str(args.plip_bin)])
            if bool(args.require_plip_success):
                sweep_cmd.append("--require-plip-success")

            sweep_returncode = 0
            try:
                _run(sweep_cmd, cwd=root)
            except subprocess.CalledProcessError as e:
                sweep_returncode = int(e.returncode)
                if sweep_returncode != 1:
                    raise

            sweep_summary = root / "processed/99_harness" / _safe_summary_name(lig.code, sweep_tag) / "summary.json"
            if not sweep_summary.exists():
                raise FileNotFoundError(f"Sweep summary not found: {sweep_summary}")
            summary = _load_json(sweep_summary)
            sweep_summaries_by_ligand[lig.code] = summary

            run_by_id = {str(r.get("run_id")): r for r in summary.get("runs", [])}
            rec_id = str(summary.get("recommended_run_id") or "")
            rec_row = run_by_id.get(rec_id, {})
            row.update(
                {
                    "status": "ok" if sweep_returncode == 0 else "no_feasible",
                    "sweep_returncode": int(sweep_returncode),
                    "prepared": {
                        "ligand_pdb": str(ligand_pdb),
                        "cg_atommap": str(cg_atommap),
                        "polar_sites": str(polar_sites),
                        "polar_coverage": str(polar_cov),
                        "site_frames": str(site_frames),
                    },
                    "sweep_summary": str(sweep_summary),
                    "n_runs": int(summary.get("sweep", {}).get("n_runs", 0)),
                    "recommended_run_id": summary.get("recommended_run_id"),
                    "recommended_feasible": bool(summary.get("recommended_feasible", False)),
                    "recommended_metrics": rec_row.get("metrics", {}),
                }
            )
        except Exception as e:
            row["error"] = repr(e)
            if bool(args.stop_on_error):
                raise
        per_ligand.append(row)

    ok_ligands = [r for r in per_ligand if r.get("sweep_summary")]
    n_ligands = len(ok_ligands)
    param_stats: dict[str, dict[str, Any]] = {}
    for row in ok_ligands:
        lig = str(row["ligand"])
        summary = sweep_summaries_by_ligand[lig]
        for run in summary.get("runs", []):
            params = {
                "score_w_coverage": float(run.get("params", {}).get("score_w_coverage", 0.0)),
                "score_w_contact": float(run.get("params", {}).get("score_w_contact", 0.0)),
                "score_w_shell": float(run.get("params", {}).get("score_w_shell", 0.0)),
                "target_res": int(run.get("params", {}).get("target_res", 0)),
                "min_cover_per_polar": int(run.get("params", {}).get("min_cover_per_polar", 0)),
            }
            key = json.dumps(params, sort_keys=True)
            p = param_stats.setdefault(
                key,
                {
                    "params": params,
                    "pass_count": 0,
                    "quality_scores": [],
                    "total_s": [],
                },
            )
            if bool(run.get("feasible")):
                p["pass_count"] += 1
            p["quality_scores"].append(float(run.get("quality_score", 0.0)))
            p["total_s"].append(float(run.get("timing_s", {}).get("total_s", 0.0)))

    param_summary: list[dict[str, Any]] = []
    for p in param_stats.values():
        q = p["quality_scores"]
        t = p["total_s"]
        row = {
            "params": p["params"],
            "pass_count": int(p["pass_count"]),
            "n_ligands": int(n_ligands),
            "pass_all": bool(n_ligands > 0 and int(p["pass_count"]) == n_ligands),
            "avg_quality_score": float(statistics.mean(q)) if q else 0.0,
            "avg_total_s": float(statistics.mean(t)) if t else 0.0,
        }
        param_summary.append(row)
    param_summary.sort(key=lambda x: (int(x["pass_count"]), float(x["avg_quality_score"])), reverse=True)
    best_param_set = param_summary[0]["params"] if param_summary else None

    summary = {
        "tag": str(args.tag),
        "inputs": {
            "repo_root": str(root),
            "python_run": str(python_run),
            "python_prep": str(python_prep),
            "vdxform_dir": str(vdxform_dir),
            "cgmap_script": str(args.cgmap_script),
            "vdm_database": str(args.vdm_database),
            "ligands": [f"{x.code}:{x.smiles}" for x in ligands],
            "polar_backend": str(args.polar_backend),
            "solver": str(args.solver),
            "top_per_site": int(args.top_per_site),
            "top_per_site_per_atom": int(args.top_per_site_per_atom),
            "chunk_size": int(args.chunk_size),
            "time_limit_s": float(args.time_limit_s),
            "ca_prefilter": float(args.ca_prefilter),
            "acceptor_model": str(args.acceptor_model),
            "score_w_coverage_list": str(args.score_w_coverage_list),
            "score_w_contact_list": str(args.score_w_contact_list),
            "score_w_shell_list": str(args.score_w_shell_list),
            "target_res_list": str(args.target_res_list),
            "min_cover_per_polar_list": str(args.min_cover_per_polar_list),
            "max_runs": int(args.max_runs),
            "run_plip": bool(args.run_plip),
            "require_plip_success": bool(args.require_plip_success),
            "plip_bin": str(args.plip_bin),
        },
        "per_ligand": per_ligand,
        "param_summary": param_summary,
        "best_param_set": best_param_set,
    }
    report_json = outdir / "summary.json"
    report_md = outdir / "summary.md"
    _write_json(report_json, summary)
    _write_md(report_md, summary)
    print(f"[done] summary: {report_json}")
    if best_param_set is not None:
        print(f"[done] best_param_set: {json.dumps(best_param_set, sort_keys=True)}")


if __name__ == "__main__":
    main()
