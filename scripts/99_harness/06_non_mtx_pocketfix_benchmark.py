#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass(frozen=True)
class LigandSpec:
    code: str
    smiles: str


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_md(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        f"# Non-MTX Pocketfix Benchmark ({summary['tag']})",
        "",
        "## Inputs",
        f"- repo_root: `{summary['inputs']['repo_root']}`",
        f"- python_run: `{summary['inputs']['python_run']}`",
        f"- python_prep: `{summary['inputs']['python_prep']}`",
        f"- vdxform_dir: `{summary['inputs']['vdxform_dir']}`",
        f"- ligands: `{','.join(summary['inputs']['ligands'])}`",
        f"- configs: `{summary['inputs']['configs']}`",
        "",
    ]
    best = summary.get("best_config")
    if best:
        lines.append(f"## Best Config: `{best}`")
        lines.append("")
    lines.extend(
        [
            "## Config Aggregate",
            "",
            "| config | pass_count | total | pass_all | avg_pocket_mean | avg_internal_worst | avg_solve_s |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.get("config_summary", []):
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.get("config")),
                    str(row.get("pass_count")),
                    str(row.get("n_ligands")),
                    str(row.get("pass_all")),
                    f"{float(row.get('avg_pocket_mean', 0.0)):.3f}",
                    f"{float(row.get('avg_internal_worst', 0.0)):.3f}",
                    f"{float(row.get('avg_solve_s', 0.0)):.3f}",
                ]
            )
            + " |"
        )
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_ligands(raw: str) -> list[LigandSpec]:
    out: list[LigandSpec] = []
    for token in [x.strip() for x in str(raw).split(",") if x.strip()]:
        if ":" not in token:
            raise ValueError(f"Invalid ligand token (expect CODE:SMILES): {token}")
        code, smiles = token.split(":", 1)
        code = code.strip().upper()
        if len(code) > 3:
            raise ValueError(f"Ligand code must be <=3 chars for PDB residue name, got: {code}")
        if not code:
            raise ValueError(f"Empty ligand code in token: {token}")
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
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        raise RuntimeError(f"RDKit embedding failed for SMILES: {smiles}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    mol = Chem.RemoveHs(mol)

    for i, atom in enumerate(mol.GetAtoms(), start=1):
        info = Chem.AtomPDBResidueInfo()
        info.SetName((atom.GetSymbol().upper() + str(i))[:4].rjust(4))
        info.SetResidueName(resname)
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(out_pdb))


def _get_valid_configs(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for row in report.get("results", []):
        name = str(row.get("name"))
        if not name:
            continue
        passed = bool(
            row.get("status") == "ok"
            and row.get("coverage_complete")
            and row.get("polar_all_satisfied")
            and row.get("ligand_clash_ok")
            and row.get("internal_clash_ok")
            and not row.get("empty_candidates")
        )
        if passed:
            out[name] = row
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Benchmark pocketfix transferability on non-MTX ligands.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--python-run", type=str, default=sys.executable, help="Python for candidate/solver/sweep scripts.")
    ap.add_argument("--python-prep", type=str, default="python3", help="Python for cgmap/openbabel prep scripts.")
    ap.add_argument(
        "--ligands",
        type=str,
        default="NIC:NC(=O)c1ccncc1,BZA:O=C(O)c1ccccc1",
        help="Comma-separated CODE:SMILES tokens for non-MTX ligands.",
    )
    ap.add_argument("--configs", type=str, default="baseline_legacy,pocketfix_default,pocketfix_weighted,pocketfix_strict")
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--solver", type=str, default="cp_sat", choices=["cp_sat", "greedy"])
    ap.add_argument("--time-limit-s", type=float, default=30.0)
    ap.add_argument("--max-internal-overlap-rank", type=float, default=0.25)
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
    ap.add_argument("--vdxform-dir", type=Path, default=Path("processed/03_vdxform"))
    args = ap.parse_args()

    root = args.repo_root.resolve()
    vdxform_dir = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()
    ligands = _parse_ligands(args.ligands)
    outdir = root / "processed/99_harness" / f"non_mtx_pocketfix_benchmark_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    sweep_py = root / "scripts/99_harness/05_mtx_pocketfix_sweep.py"
    extract_py = root / "scripts/02_polar_sites/01_extract_polar_sites.py"
    cov_py = root / "scripts/02_polar_sites/02_check_polar_coverage_vs_cgmap.py"
    frames_py = root / "scripts/02_polar_sites/03_build_ligand_site_frames.py"
    cg_frame_defs = root / "configs/cg_frame_defs.json"

    for need in [sweep_py, extract_py, cov_py, frames_py, args.cgmap_script, args.vdm_database, vdxform_dir]:
        if not Path(need).exists():
            raise FileNotFoundError(f"Missing required path: {need}")

    per_ligand_reports: list[dict[str, Any]] = []
    for lig in ligands:
        lig_dir = outdir / lig.code
        lig_dir.mkdir(parents=True, exist_ok=True)
        ligand_pdb = lig_dir / f"{lig.code}.pdb"
        cg_atommap = lig_dir / f"{lig.code}_cg_atommap.json"
        polar_sites = lig_dir / f"{lig.code}_polar_sites.json"
        polar_cov = lig_dir / f"{lig.code}_polar_coverage.json"
        site_frames = lig_dir / f"{lig.code}_site_frames.json"

        row: dict[str, Any] = {"ligand": lig.code, "smiles": lig.smiles}
        try:
            _build_ligand_pdb(lig.smiles, lig.code, ligand_pdb)

            _run(
                [
                    str(args.python_prep),
                    str(args.cgmap_script),
                    str(ligand_pdb),
                    "--format",
                    "json",
                    "--vdm-database",
                    str(args.vdm_database),
                    "-o",
                    str(cg_atommap),
                ]
            )

            try:
                _run([str(args.python_prep), str(extract_py), str(ligand_pdb), "-o", str(polar_sites), "--backend", "openbabel"])
                row["polar_backend"] = "openbabel"
            except Exception:
                _run([str(args.python_run), str(extract_py), str(ligand_pdb), "-o", str(polar_sites), "--backend", "rdkit"])
                row["polar_backend"] = "rdkit"

            _run([str(args.python_run), str(cov_py), "--polar-sites", str(polar_sites), "--cg-atommap", str(cg_atommap), "-o", str(polar_cov)])
            _run(
                [
                    str(args.python_run),
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
                ]
            )

            site_meta = _load_json(site_frames)
            if int(site_meta.get("n_sites", 0)) <= 0:
                raise RuntimeError("Prepared site_frames has no sites.")

            sweep_tag = f"{args.tag}_{lig.code.lower()}"
            _run(
                [
                    str(args.python_run),
                    str(sweep_py),
                    "--repo-root",
                    str(root),
                    "--tag",
                    sweep_tag,
                    "--python",
                    str(args.python_run),
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
                    "--configs",
                    str(args.configs),
                    "--top-per-site",
                    str(int(args.top_per_site)),
                    "--solver",
                    str(args.solver),
                    "--time-limit-s",
                    str(float(args.time_limit_s)),
                    "--max-internal-overlap-rank",
                    str(float(args.max_internal_overlap_rank)),
                ]
            )
            report_path = root / "processed/99_harness" / f"mtx_pocketfix_sweep_{sweep_tag}" / "report.json"
            report = _load_json(report_path)
            row["status"] = "ok"
            row["sweep_report"] = str(report_path)
            row["best"] = report.get("best", {}).get("name")
            row["valid_configs"] = sorted(_get_valid_configs(report).keys())
            row["prepared"] = {
                "ligand_pdb": str(ligand_pdb),
                "cg_atommap": str(cg_atommap),
                "polar_sites": str(polar_sites),
                "site_frames": str(site_frames),
            }
        except Exception as e:
            row["status"] = "error"
            row["error"] = repr(e)
        per_ligand_reports.append(row)

    config_names = [x.strip() for x in str(args.configs).split(",") if x.strip()]
    config_summary: list[dict[str, Any]] = []
    best_config = None
    best_key = None
    for cfg in config_names:
        pass_count = 0
        pocket_vals: list[float] = []
        internal_vals: list[float] = []
        solve_vals: list[float] = []
        for row in per_ligand_reports:
            if row.get("status") != "ok":
                continue
            report = _load_json(Path(str(row["sweep_report"])))
            valid = _get_valid_configs(report)
            item = valid.get(cfg)
            if item is None:
                continue
            pass_count += 1
            pocket_vals.append(float(item.get("selected_pocket_score_mean", 0.0)))
            internal_vals.append(float(item.get("internal_worst_overlap", 999.0)))
            solve_vals.append(float(item.get("solve_s", 9999.0)))

        n_lig = len([r for r in per_ligand_reports if r.get("status") == "ok"])
        avg_pocket = float(np.mean(pocket_vals)) if pocket_vals else 0.0
        avg_internal = float(np.mean(internal_vals)) if internal_vals else 999.0
        avg_solve = float(np.mean(solve_vals)) if solve_vals else 9999.0
        pass_all = bool(n_lig > 0 and pass_count == n_lig)
        row = {
            "config": cfg,
            "pass_count": int(pass_count),
            "n_ligands": int(n_lig),
            "pass_all": pass_all,
            "avg_pocket_mean": avg_pocket,
            "avg_internal_worst": avg_internal,
            "avg_solve_s": avg_solve,
        }
        config_summary.append(row)
        key = (int(pass_count), avg_pocket, -avg_internal, -avg_solve)
        if best_key is None or key > best_key:
            best_key = key
            best_config = cfg

    summary = {
        "tag": str(args.tag),
        "inputs": {
            "repo_root": str(root),
            "python_run": str(args.python_run),
            "python_prep": str(args.python_prep),
            "vdxform_dir": str(vdxform_dir),
            "cgmap_script": str(args.cgmap_script),
            "vdm_database": str(args.vdm_database),
            "ligands": [f"{x.code}:{x.smiles}" for x in ligands],
            "configs": str(args.configs),
            "top_per_site": int(args.top_per_site),
            "solver": str(args.solver),
            "time_limit_s": float(args.time_limit_s),
            "max_internal_overlap_rank": float(args.max_internal_overlap_rank),
        },
        "per_ligand": per_ligand_reports,
        "config_summary": config_summary,
        "best_config": best_config,
    }

    report_json = outdir / "report.json"
    report_md = outdir / "report.md"
    _write_json(report_json, summary)
    _write_md(report_md, summary)
    print(f"[done] report: {report_json}")
    if best_config is not None:
        print(f"[done] best_config: {best_config}")


if __name__ == "__main__":
    main()
