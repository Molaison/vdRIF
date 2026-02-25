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
        f"# Non-MTX Cation Typing Validation ({summary['tag']})",
        "",
        "## Inputs",
        f"- repo_root: `{summary['inputs']['repo_root']}`",
        f"- ligands: `{','.join(summary['inputs']['ligands'])}`",
        f"- openbabel_python: `{summary['inputs']['openbabel_python']}`",
        f"- plip_bin: `{summary['inputs']['plip_bin']}`",
        "",
        "## Summary",
        f"- tested: {summary['metrics']['n_tested']}",
        f"- solved: {summary['metrics']['n_solved']}",
        f"- plip_compared: {summary['metrics']['n_plip_compared']}",
        f"- ligands_with_dropped_cation_atoms: {summary['metrics']['n_with_dropped_cation_atoms']}",
        f"- old_fail_new_pass: {summary['metrics']['n_old_fail_new_pass']}",
        "",
        "## Per Ligand",
        "",
        "| ligand | status | cand_n | solve_ok | dropped_cation_atoms | plip_new | plip_old | clash_ok | internal_ok |",
        "|---|---|---:|---:|---:|---|---|---:|---:|",
    ]
    for row in summary.get("per_ligand", []):
        plip_new = row.get("plip_new")
        plip_old = row.get("plip_old")
        def _fmt_plip(x: Any) -> str:
            if not isinstance(x, dict):
                return "-"
            return f"{int(x.get('n_satisfied', 0))}/{int(x.get('n_polar_atoms', 0))}"

        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.get("ligand")),
                    str(row.get("status")),
                    str(int(row.get("cand_n", 0))),
                    str(bool(row.get("solve_ok"))).lower(),
                    str(int(len(row.get("dropped_cation_atoms") or []))),
                    _fmt_plip(plip_new),
                    _fmt_plip(plip_old),
                    str(bool(row.get("clash_ok"))).lower(),
                    str(bool(row.get("internal_ok"))).lower(),
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
        smiles = smiles.strip()
        if not code:
            raise ValueError(f"Empty ligand code in token: {token}")
        if len(code) > 3:
            raise ValueError(f"Ligand code must be <=3 chars for PDB residue name, got: {code}")
        if not smiles:
            raise ValueError(f"Empty SMILES in token: {token}")
        out.append(LigandSpec(code=code, smiles=smiles))
    if not out:
        raise ValueError("No ligands provided.")
    return out


def _build_ligand_pdb_and_metadata(smiles: str, resname: str, out_pdb: Path) -> tuple[Chem.Mol, list[dict[str, Any]]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        raise RuntimeError(f"RDKit embedding failed for SMILES: {smiles}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    mol = Chem.RemoveHs(mol)

    conf = mol.GetConformer()
    excluded_nplus_lt3: list[dict[str, Any]] = []
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        info = Chem.AtomPDBResidueInfo()
        atom_name = (atom.GetSymbol().upper() + str(i))[:4].rjust(4)
        info.SetName(atom_name)
        info.SetResidueName(resname)
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

        if int(atom.GetFormalCharge()) > 0 and int(atom.GetAtomicNum()) == 7 and int(atom.GetDegree()) < 3:
            p = conf.GetAtomPosition(atom.GetIdx())
            excluded_nplus_lt3.append(
                {
                    "atom_name": atom_name.strip(),
                    "atom_idx": int(atom.GetIdx()) + 1,  # OpenBabel index is 1-based.
                    "element": atom.GetSymbol(),
                    "formal_charge": int(atom.GetFormalCharge()),
                    "xyz": [float(p.x), float(p.y), float(p.z)],
                    "degree": int(atom.GetDegree()),
                }
            )

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(out_pdb))
    return mol, excluded_nplus_lt3


def _collect_formal_positive_atoms(mol: Chem.Mol) -> list[dict[str, Any]]:
    conf = mol.GetConformer()
    out: list[dict[str, Any]] = []
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        chg = int(atom.GetFormalCharge())
        if chg <= 0:
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        out.append(
            {
                "atom_name": (atom.GetSymbol().upper() + str(i))[:4].rjust(4).strip(),
                "atom_idx": int(atom.GetIdx()) + 1,  # OpenBabel index is 1-based.
                "element": atom.GetSymbol(),
                "formal_charge": chg,
                "xyz": [float(p.x), float(p.y), float(p.z)],
            }
        )
    return out


def _build_legacy_polar_from_current(
    new_polar_path: Path,
    formal_positive_atoms: list[dict[str, Any]],
    out_path: Path,
) -> tuple[dict[str, Any], list[str]]:
    payload = _load_json(new_polar_path)
    sites = list(payload.get("sites", []))
    by_name: dict[str, dict[str, Any]] = {str(s["atom_name"]): s for s in sites}
    dropped_names: list[str] = []
    for atom in formal_positive_atoms:
        nm = str(atom["atom_name"])
        if nm in by_name:
            roles = set(str(x) for x in by_name[nm].get("roles", []))
            if "cation" not in roles:
                dropped_names.append(nm)
            roles.add("cation")
            by_name[nm]["roles"] = sorted(roles)
            continue
        dropped_names.append(nm)
        by_name[nm] = {
            "atom_name": nm,
            "atom_idx": int(atom["atom_idx"]),
            "element": str(atom["element"]),
            "formal_charge": int(atom["formal_charge"]),
            "xyz": [float(x) for x in atom["xyz"]],
            "roles": ["cation"],
        }

    out_sites = sorted(by_name.values(), key=lambda s: (str(s.get("atom_name")), int(s.get("atom_idx", 0))))
    out_payload = {
        "input": dict(payload.get("input", {})),
        "n_sites": len(out_sites),
        "sites": out_sites,
        "legacy_cation_patch_atoms": dropped_names,
        "legacy_cation_rule": "formal_charge>0 implies cation on all atoms",
    }
    _write_json(out_path, out_payload)
    return out_payload, sorted(set(dropped_names))


def _count_candidates(npz_path: Path) -> int:
    z = np.load(npz_path, allow_pickle=True)
    return int(z["cand_id_u64"].shape[0])


def _safe_run(cmd: list[str]) -> tuple[bool, str | None]:
    try:
        _run(cmd)
        return True, None
    except subprocess.CalledProcessError as e:
        return False, f"CalledProcessError(returncode={e.returncode})"


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate conservative cation typing on non-MTX ionic ligands.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument(
        "--ligands",
        type=str,
        default="PYM:C[n+]1ccccc1,PYA:c1cc[n+](C)cc1,TRM:C[N+](C)(C)C,BIM:C[n+]1cnccn1,GUA:[NH2+]C(=N)N",
        help="Comma-separated CODE:SMILES tokens.",
    )
    ap.add_argument("--python-run", type=str, default=sys.executable, help="Python for uv-run scripts.")
    ap.add_argument(
        "--openbabel-python",
        type=str,
        default=str(Path.home() / "miniconda/bin/python"),
        help="Python interpreter with working openbabel+rdkit for --backend openbabel extraction.",
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
    ap.add_argument("--vdxform-dir", type=Path, default=Path("processed/03_vdxform"))
    ap.add_argument("--plip-bin", type=str, default=str(Path.home() / "miniconda/bin/plip"))
    ap.add_argument("--top-per-site", type=int, default=120)
    ap.add_argument("--chunk-size", type=int, default=2000)
    ap.add_argument("--min-res", type=int, default=4)
    ap.add_argument("--max-res", type=int, default=12)
    ap.add_argument("--time-limit-s", type=float, default=20.0)
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"non_mtx_cation_typing_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)
    ligands = _parse_ligands(args.ligands)

    extract_py = root / "scripts/02_polar_sites/01_extract_polar_sites.py"
    frames_py = root / "scripts/02_polar_sites/03_build_ligand_site_frames.py"
    cands_py = root / "scripts/04_candidates/01_generate_candidates.py"
    solve_py = root / "scripts/05_solver/01_solve_motif.py"
    v_clash_py = root / "scripts/05_solver/04_validate_motif_clashes.py"
    v_internal_py = root / "scripts/05_solver/05_validate_motif_internal_clashes.py"
    v_plip_py = root / "scripts/05_solver/06_validate_motif_plip.py"

    vdx = (root / args.vdxform_dir).resolve() if not args.vdxform_dir.is_absolute() else args.vdxform_dir.resolve()
    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        p = vdx / cg / f"vdxform_{cg}.npz"
        if not p.exists():
            raise FileNotFoundError(f"Missing required vdxform file: {p}")

    for need in [extract_py, frames_py, cands_py, solve_py, v_clash_py, v_internal_py, v_plip_py, args.cgmap_script]:
        if not Path(need).exists():
            raise FileNotFoundError(f"Missing required script/path: {need}")

    per_ligand: list[dict[str, Any]] = []
    for lig in ligands:
        lig_dir = outdir / lig.code
        lig_dir.mkdir(parents=True, exist_ok=True)
        row: dict[str, Any] = {"ligand": lig.code, "smiles": lig.smiles, "status": "init"}

        ligand_pdb = lig_dir / f"{lig.code}.pdb"
        cgmap = lig_dir / f"{lig.code}_cgmap.json"
        polar_new = lig_dir / f"{lig.code}_polar_new.json"
        polar_old = lig_dir / f"{lig.code}_polar_old.json"
        frames = lig_dir / f"{lig.code}_frames.json"
        cands_prefix = lig_dir / f"{lig.code}_candidates"
        motif_json = lig_dir / f"{lig.code}_motif.json"
        motif_pdb = lig_dir / f"{lig.code}_motif.pdb"
        clash_json = lig_dir / f"{lig.code}_motif_clash_validation.json"
        internal_json = lig_dir / f"{lig.code}_motif_internal_clash_validation.json"
        plip_new_json = lig_dir / f"{lig.code}_motif_plip_new.json"
        plip_old_json = lig_dir / f"{lig.code}_motif_plip_old.json"

        try:
            mol, excluded_atoms = _build_ligand_pdb_and_metadata(lig.smiles, lig.code, ligand_pdb)
            row["excluded_nplus_degree_lt3_atoms"] = excluded_atoms
            formal_positive_atoms = _collect_formal_positive_atoms(mol)

            _run(
                [
                    "python3",
                    str(args.cgmap_script),
                    str(ligand_pdb),
                    "--format",
                    "json",
                    "--vdm-database",
                    str(args.vdm_database),
                    "-o",
                    str(cgmap),
                ]
            )

            _run(
                [
                    str(args.openbabel_python),
                    str(extract_py),
                    str(ligand_pdb),
                    "--backend",
                    "openbabel",
                    "-o",
                    str(polar_new),
                ]
            )

            _, dropped_names = _build_legacy_polar_from_current(polar_new, formal_positive_atoms, polar_old)
            row["dropped_cation_atoms"] = dropped_names

            _run(
                [
                    str(args.python_run),
                    str(frames_py),
                    str(ligand_pdb),
                    "--cg-atommap",
                    str(cgmap),
                    "--polar-sites",
                    str(polar_new),
                    "--add-bb-cnh-for-uncovered-donors",
                    "--add-ccn-for-cations",
                    "-o",
                    str(frames),
                ]
            )

            _run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(cands_py),
                    "--ligand-pdb",
                    str(ligand_pdb),
                    "--polar-sites",
                    str(polar_new),
                    "--site-frames",
                    str(frames),
                    "--vdxform-dir",
                    str(vdx),
                    "--out-prefix",
                    str(cands_prefix),
                    "--top-per-site",
                    str(int(args.top_per_site)),
                    "--chunk-size",
                    str(int(args.chunk_size)),
                    "--acceptor-model",
                    "plip",
                    "--min-sidechain-ligand-contacts",
                    "0",
                    "--max-sidechain-min-dist",
                    "6.0",
                    "--score-weight-contact-count",
                    "0.0",
                ]
            )

            cand_n = _count_candidates(cands_prefix.with_suffix(".npz"))
            row["cand_n"] = cand_n
            if cand_n <= 0:
                row["status"] = "empty_candidates"
                per_ligand.append(row)
                continue

            ok_solve, solve_err = _safe_run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(solve_py),
                    "--solver",
                    "greedy",
                    "--candidates-npz",
                    str(cands_prefix.with_suffix(".npz")),
                    "--candidates-meta",
                    str(cands_prefix.with_suffix(".json")),
                    "--ligand-pdb",
                    str(ligand_pdb),
                    "--out-json",
                    str(motif_json),
                    "--out-pdb",
                    str(motif_pdb),
                    "--min-res",
                    str(int(args.min_res)),
                    "--max-res",
                    str(int(args.max_res)),
                    "--num-workers",
                    "1",
                    "--time-limit-s",
                    str(float(args.time_limit_s)),
                ]
            )
            row["solve_ok"] = bool(ok_solve)
            if not ok_solve:
                row["status"] = "solve_failed"
                row["solve_error"] = solve_err
                per_ligand.append(row)
                continue

            _run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(v_clash_py),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    lig.code,
                    "-o",
                    str(clash_json),
                ]
            )
            _run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(v_internal_py),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--ligand-resname",
                    lig.code,
                    "-o",
                    str(internal_json),
                ]
            )
            row["clash_ok"] = bool(_load_json(clash_json).get("ok"))
            row["internal_ok"] = bool(_load_json(internal_json).get("ok"))

            _run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(v_plip_py),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--polar-sites",
                    str(polar_new),
                    "--plip-bin",
                    str(args.plip_bin),
                    "-o",
                    str(plip_new_json),
                ]
            )
            _run(
                [
                    "uv",
                    "run",
                    "-p",
                    "3.11",
                    "python",
                    str(v_plip_py),
                    "--motif-pdb",
                    str(motif_pdb),
                    "--polar-sites",
                    str(polar_old),
                    "--plip-bin",
                    str(args.plip_bin),
                    "-o",
                    str(plip_old_json),
                ]
            )

            newj = _load_json(plip_new_json)
            oldj = _load_json(plip_old_json)
            row["plip_new"] = {
                "all_satisfied": bool(newj.get("all_satisfied")),
                "n_polar_atoms": int(newj.get("n_polar_atoms", 0)),
                "n_satisfied": int(newj.get("n_satisfied", 0)),
                "unsatisfied_polar_atoms": list(newj.get("unsatisfied_polar_atoms") or []),
            }
            row["plip_old"] = {
                "all_satisfied": bool(oldj.get("all_satisfied")),
                "n_polar_atoms": int(oldj.get("n_polar_atoms", 0)),
                "n_satisfied": int(oldj.get("n_satisfied", 0)),
                "unsatisfied_polar_atoms": list(oldj.get("unsatisfied_polar_atoms") or []),
            }
            row["status"] = "ok"
        except Exception as e:
            row["status"] = "error"
            row["error"] = repr(e)

        per_ligand.append(row)

    solved_rows = [r for r in per_ligand if bool(r.get("solve_ok"))]
    plip_rows = [r for r in solved_rows if isinstance(r.get("plip_new"), dict) and isinstance(r.get("plip_old"), dict)]
    summary = {
        "tag": str(args.tag),
        "inputs": {
            "repo_root": str(root),
            "ligands": [f"{x.code}:{x.smiles}" for x in ligands],
            "python_run": str(args.python_run),
            "openbabel_python": str(args.openbabel_python),
            "cgmap_script": str(args.cgmap_script),
            "vdm_database": str(args.vdm_database),
            "vdxform_dir": str(vdx),
            "plip_bin": str(args.plip_bin),
        },
        "per_ligand": per_ligand,
        "metrics": {
            "n_tested": len(per_ligand),
            "n_solved": len(solved_rows),
            "n_plip_compared": len(plip_rows),
            "n_with_dropped_cation_atoms": int(sum(1 for r in per_ligand if len(r.get("dropped_cation_atoms") or []) > 0)),
            "n_old_fail_new_pass": int(
                sum(
                    1
                    for r in plip_rows
                    if bool(r["plip_new"].get("all_satisfied")) and not bool(r["plip_old"].get("all_satisfied"))
                )
            ),
        },
    }

    report_json = outdir / "report.json"
    report_md = outdir / "report.md"
    _write_json(report_json, summary)
    _write_md(report_md, summary)
    print(f"[done] report: {report_json}")


if __name__ == "__main__":
    main()
