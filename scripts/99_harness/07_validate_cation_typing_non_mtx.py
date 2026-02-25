#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem


CORE_LIGANDS: list[dict[str, str]] = [
    {"code": "GDM", "smiles": "NC(=[NH2+])N"},
    {"code": "BZA", "smiles": "NC(=[NH2+])c1ccccc1"},
    {"code": "BTM", "smiles": "C[N+](C)(C)Cc1ccccc1"},
]

DROP_PROBE: dict[str, str] = {"code": "DZN", "smiles": "c1ccccc1[N+]#N"}


def _run(cmd: list[str], *, cwd: Path) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd))


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_pdb_from_smiles(smiles: str, code: str, out_pdb: Path) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES for {code}: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0) != 0:
        raise RuntimeError(f"RDKit embed failed for {code}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=300)
    block = Chem.MolToPDBBlock(mol)

    lines: list[str] = []
    serial = 1
    resname = code[:3].upper()
    for line in block.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            line = "HETATM" + f"{serial:5d}" + line[11:17] + f"{resname:>3}" + line[20:]
            serial += 1
        lines.append(line)
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _extract_legacy_openbabel(pdb_path: Path) -> dict[str, Any]:
    conv = ob.OBConversion()
    if not conv.SetInFormat("pdb"):
        raise RuntimeError("OpenBabel: failed to set PDB input format")
    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(pdb_path)):
        raise ValueError(f"OpenBabel: failed to read ligand PDB: {pdb_path}")
    if mol.NumAtoms() <= 0:
        raise ValueError(f"OpenBabel: ligand PDB has no atoms: {pdb_path}")

    sites: list[dict[str, Any]] = []
    for atom in ob.OBMolAtomIter(mol):
        if int(atom.GetAtomicNum()) == 1:
            continue
        res = atom.GetResidue()
        atom_name = res.GetAtomID(atom).strip() if res is not None else f"ATOM{atom.GetIdx()}"
        roles: set[str] = set()
        if atom.IsHbondDonor():
            roles.add("donor")
        if atom.IsHbondAcceptor():
            roles.add("acceptor")
        chg = int(atom.GetFormalCharge())
        if chg > 0:
            roles.add("cation")
        elif chg < 0:
            roles.add("anion")
        if not roles:
            continue
        sites.append(
            {
                "atom_name": atom_name,
                "atom_idx": int(atom.GetIdx()),
                "element": ob.GetSymbol(int(atom.GetAtomicNum())),
                "formal_charge": chg,
                "xyz": [float(atom.GetX()), float(atom.GetY()), float(atom.GetZ())],
                "roles": sorted(roles),
            }
        )
    sites = sorted(sites, key=lambda s: (str(s["atom_name"]), int(s["atom_idx"])))
    return {"input": {"ligand_pdb": str(pdb_path), "backend": "openbabel_legacy"}, "n_sites": len(sites), "sites": sites}


def _role_map(payload: dict[str, Any]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for s in payload.get("sites", []):
        out[str(s["atom_name"])] = set(str(r) for r in s.get("roles", []))
    return out


def _compare_role_payloads(old_payload: dict[str, Any], new_payload: dict[str, Any]) -> dict[str, Any]:
    old_map = _role_map(old_payload)
    new_map = _role_map(new_payload)
    atoms = sorted(set(old_map) | set(new_map))

    dropped_cation_atoms: list[str] = []
    added_cation_atoms: list[str] = []
    diffs: list[dict[str, Any]] = []

    for atom in atoms:
        old_roles = sorted(old_map.get(atom, set()))
        new_roles = sorted(new_map.get(atom, set()))
        if old_roles != new_roles:
            diffs.append({"atom_name": atom, "old_roles": old_roles, "new_roles": new_roles})
        if "cation" in old_map.get(atom, set()) and "cation" not in new_map.get(atom, set()):
            dropped_cation_atoms.append(atom)
        if "cation" not in old_map.get(atom, set()) and "cation" in new_map.get(atom, set()):
            added_cation_atoms.append(atom)

    return {
        "n_old_sites": int(old_payload.get("n_sites", 0)),
        "n_new_sites": int(new_payload.get("n_sites", 0)),
        "n_role_diffs": len(diffs),
        "role_diffs": diffs,
        "dropped_cation_atoms": dropped_cation_atoms,
        "added_cation_atoms": added_cation_atoms,
        "old_roles_by_atom": {k: sorted(v) for k, v in sorted(old_map.items())},
        "new_roles_by_atom": {k: sorted(v) for k, v in sorted(new_map.items())},
    }


def _extract_new_polar_openbabel(repo_root: Path, ligand_pdb: Path, out_json: Path) -> None:
    _run(
        [
            sys.executable,
            str(repo_root / "scripts/02_polar_sites/01_extract_polar_sites.py"),
            str(ligand_pdb),
            "-o",
            str(out_json),
            "--backend",
            "openbabel",
        ],
        cwd=repo_root,
    )


def _generate_cgmap(
    ligand_pdb: Path,
    out_json: Path,
    *,
    cgmap_script: Path,
    vdm_db: Path,
    repo_root: Path,
) -> None:
    _run(
        [
            sys.executable,
            str(cgmap_script),
            str(ligand_pdb),
            "--format",
            "json",
            "--vdm-database",
            str(vdm_db),
            "-o",
            str(out_json),
        ],
        cwd=repo_root,
    )


def _run_pipeline_with_polar(
    *,
    repo_root: Path,
    py311: Path,
    vdxform_dir: Path,
    ligand_pdb: Path,
    cgmap_json: Path,
    polar_json: Path,
    run_dir: Path,
    ligand_code: str,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    frames = run_dir / f"{ligand_code}_site_frames.json"
    cand_prefix = run_dir / f"{ligand_code}_candidates"
    motif_json = run_dir / f"{ligand_code}_motif.json"
    motif_pdb = run_dir / f"{ligand_code}_motif.pdb"
    plip_json = run_dir / f"{ligand_code}_plip.json"

    _run(
        [
            str(py311),
            str(repo_root / "scripts/02_polar_sites/03_build_ligand_site_frames.py"),
            str(ligand_pdb),
            "--cg-atommap",
            str(cgmap_json),
            "--polar-sites",
            str(polar_json),
            "--cg-frame-defs",
            str(repo_root / "configs/cg_frame_defs.json"),
            "--add-bb-cnh-for-uncovered-donors",
            "--add-ccn-for-cations",
            "-o",
            str(frames),
        ],
        cwd=repo_root,
    )

    _run(
        [
            str(py311),
            str(repo_root / "scripts/04_candidates/01_generate_candidates.py"),
            "--ligand-pdb",
            str(ligand_pdb),
            "--polar-sites",
            str(polar_json),
            "--site-frames",
            str(frames),
            "--vdxform-dir",
            str(vdxform_dir),
            "--out-prefix",
            str(cand_prefix),
            "--chunk-size",
            "5000",
            "--top-per-site",
            "200",
            "--top-per-site-per-atom",
            "50",
            "--acceptor-model",
            "plip",
            "--min-lig-donor-angle-deg",
            "60",
            "--pocket-contact-cutoff",
            "4.5",
            "--min-pocket-sidechain-contacts",
            "0",
            "--pocket-contact-weight",
            "0.2",
            "--clash-tol",
            "0.5",
            "--exclude-aa3",
            "PRO,CYS",
            "--require-sidechain-facing",
            "--min-sidechain-facing-dot",
            "0.2",
            "--min-sidechain-centroid-dot",
            "0.0",
            "--allow-backbone-hbonds",
            "--require-full-coverage",
        ],
        cwd=repo_root,
    )

    _run(
        [
            str(py311),
            str(repo_root / "scripts/05_solver/01_solve_motif.py"),
            "--candidates-npz",
            str(cand_prefix.with_suffix(".npz")),
            "--candidates-meta",
            str(cand_prefix.with_suffix(".json")),
            "--ligand-pdb",
            str(ligand_pdb),
            "--out-json",
            str(motif_json),
            "--out-pdb",
            str(motif_pdb),
            "--solver",
            "greedy",
            "--min-res",
            "3",
            "--max-res",
            "12",
            "--time-limit-s",
            "120",
            "--num-workers",
            "1",
            "--grid-size",
            "4.0",
            "--ca-prefilter",
            "12.0",
            "--clash-tol",
            "0.5",
        ],
        cwd=repo_root,
    )

    _run(
        [
            str(py311),
            str(repo_root / "scripts/05_solver/06_validate_motif_plip.py"),
            "--motif-pdb",
            str(motif_pdb),
            "--polar-sites",
            str(polar_json),
            "-o",
            str(plip_json),
            "--maxthreads",
            "1",
            "--timeout-s",
            "180",
        ],
        cwd=repo_root,
    )

    cand = _read_json(cand_prefix.with_suffix(".json"))
    motif = _read_json(motif_json)
    plip = _read_json(plip_json)
    return {
        "status": "ok",
        "paths": {
            "site_frames": str(frames),
            "candidates_json": str(cand_prefix.with_suffix(".json")),
            "candidates_npz": str(cand_prefix.with_suffix(".npz")),
            "motif_json": str(motif_json),
            "motif_pdb": str(motif_pdb),
            "plip_json": str(plip_json),
        },
        "metrics": {
            "n_candidates": int(cand.get("n_candidates", 0)),
            "n_selected": int(motif.get("n_selected", 0)),
            "coverage_complete": bool(motif.get("coverage_complete", False)),
            "plip_all_satisfied": bool(plip.get("all_satisfied", False)),
            "plip_n_satisfied": int(plip.get("n_satisfied", 0)),
            "plip_n_polar": int(plip.get("n_polar_atoms", 0)),
            "plip_unsatisfied_polar_atoms": list(plip.get("unsatisfied_polar_atoms", [])),
        },
    }


def _make_single_atom_polar(in_json: Path, atom_name: str, out_json: Path) -> None:
    payload = _read_json(in_json)
    sites = [s for s in payload.get("sites", []) if str(s.get("atom_name")) == atom_name]
    if len(sites) != 1:
        raise ValueError(f"Expected exactly 1 site for atom {atom_name} in {in_json}, got {len(sites)}")
    out = {"input": payload.get("input", {}), "n_sites": 1, "sites": sites}
    _write_json(out_json, out)


def _parse_int_list_text(text: str) -> list[int]:
    vals: list[int] = []
    for x in text.split(","):
        x = x.strip()
        if not x:
            continue
        vals.append(int(x))
    return vals


def _read_lig_idx_list(parent: ET.Element, tag: str) -> list[int]:
    el = parent.find(tag)
    if el is None:
        return []
    kids = el.findall("idx")
    if kids:
        out: list[int] = []
        for k in kids:
            if k.text is None:
                continue
            txt = k.text.strip()
            if txt:
                out.append(int(txt))
        return out
    if el.text is None:
        return []
    txt = el.text.strip()
    if not txt:
        return []
    return _parse_int_list_text(txt)


def _saltbridge_atom_names_from_plip(plip_json_path: Path) -> list[str]:
    pl = _read_json(plip_json_path)
    report_xml = Path(pl["plip"]["report_xml"])
    plipfixed = Path(pl["plip"]["plipfixed_pdb"])
    lig_info = pl.get("ligand", {})
    lig_resname = str(lig_info.get("plip_resname", "")).strip()
    lig_chain = str(lig_info.get("plip_chain", "")).strip()
    lig_resnum = int(lig_info.get("plip_resnum", 0))

    serial_to_name: dict[int, str] = {}
    for line in plipfixed.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resnum = int(line[22:26].strip())
        if resname == lig_resname and chain == lig_chain and resnum == lig_resnum:
            serial_to_name[serial] = atom_name

    tree = ET.parse(str(report_xml))
    root = tree.getroot()
    salt_serials: list[int] = []
    for sb in root.findall(".//salt_bridge"):
        salt_serials.extend(_read_lig_idx_list(sb, "lig_idx_list"))

    names = sorted({serial_to_name[i] for i in salt_serials if i in serial_to_name})
    return names


def _write_summary_md(summary: dict[str, Any], out_md: Path) -> None:
    core = summary["core_ligands"]
    lines: list[str] = [
        "# Non-MTX Conservative Cation Typing Validation",
        "",
        f"- tag: `{summary['tag']}`",
        f"- outdir: `{summary['outdir']}`",
        f"- py311: `{summary['config']['python311']}`",
        f"- vdxform_dir: `{summary['config']['vdxform_dir']}`",
        "",
        "## Core Non-MTX PLIP Validation (new rule)",
        "",
        "| ligand | n_old_sites | n_new_sites | dropped_cation_atoms | n_candidates | n_selected | coverage_complete | plip_all_satisfied | plip_satisfied / n_polar |",
        "|---|---:|---:|---|---:|---:|---|---|---|",
    ]
    for rec in core:
        cmp = rec["role_compare"]
        run = rec["new_pipeline"]
        met = run.get("metrics", {})
        lines.append(
            "| "
            + f"{rec['code']} | {cmp['n_old_sites']} | {cmp['n_new_sites']} | "
            + f"{','.join(cmp['dropped_cation_atoms']) if cmp['dropped_cation_atoms'] else '-'} | "
            + f"{met.get('n_candidates', 0)} | {met.get('n_selected', 0)} | "
            + f"{met.get('coverage_complete', False)} | {met.get('plip_all_satisfied', False)} | "
            + f"{met.get('plip_n_satisfied', 0)} / {met.get('plip_n_polar', 0)} |"
        )

    probe = summary["dropped_cation_probe"]
    lines.extend(
        [
            "",
            "## Dropped-Cation Probe",
            "",
            f"- probe ligand: `{probe['code']}`",
            f"- dropped_cation_atoms: `{probe['role_compare']['dropped_cation_atoms']}`",
        ]
    )
    if probe.get("focus_atom") is not None:
        lines.extend(
            [
                f"- focus_atom: `{probe['focus_atom']}`",
                f"- legacy focus roles: `{probe['focus_legacy_roles']}`",
                f"- new focus roles: `{probe['focus_new_roles']}`",
                f"- legacy salt-bridge ligand atoms: `{probe['legacy_focus_saltbridge_atoms']}`",
                f"- new salt-bridge ligand atoms: `{probe['new_focus_saltbridge_atoms']}`",
                f"- dropped atoms with legacy salt-bridge evidence: `{probe['dropped_atoms_with_legacy_saltbridge']}`",
            ]
        )
    else:
        lines.append("- No dropped cation atom found in probe ligand.")

    acc = summary["acceptance_check"]
    lines.extend(
        [
            "",
            "## Acceptance Check",
            "",
            f"- plip_runs_on_core_non_mtx: {acc['plip_runs_on_core_non_mtx']}",
            f"- plip_all_satisfied_on_core_non_mtx: {acc['plip_all_satisfied_on_core_non_mtx']}",
            f"- has_minimum_three_plip_runs: {acc['has_minimum_three_plip_runs']}",
            f"- total_dropped_cation_atoms_across_core: {acc['total_dropped_cation_atoms_across_core']}",
            f"- dropped_cation_has_true_saltbridge_evidence: {acc['dropped_cation_has_true_saltbridge_evidence']}",
        ]
    )

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate conservative cation typing on non-MTX ligands.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=time.strftime("%Y%m%d-%H%M%S"))
    ap.add_argument("--python311", type=Path, required=True, help="Python 3.11 env for candidate/solver scripts.")
    ap.add_argument("--vdxform-dir", type=Path, required=True)
    ap.add_argument("--cgmap-script", type=Path, required=True)
    ap.add_argument("--vdm-db", type=Path, required=True)
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    py311 = args.python311 if args.python311.is_absolute() else (repo_root / args.python311)
    vdxform_dir = args.vdxform_dir.resolve()
    cgmap_script = args.cgmap_script.resolve()
    vdm_db = args.vdm_db.resolve()

    outdir = repo_root / "processed/99_harness" / f"non_mtx_cation_validate_{args.tag}"
    lig_dir = outdir / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)

    core_records: list[dict[str, Any]] = []
    for lig in CORE_LIGANDS:
        code = lig["code"]
        smiles = lig["smiles"]
        ldir = lig_dir / code
        ldir.mkdir(parents=True, exist_ok=True)
        ligand_pdb = ldir / f"{code}.pdb"
        cgmap_json = ldir / f"{code}_cgmap.json"
        polar_new = ldir / f"{code}_polar_new.json"
        polar_legacy = ldir / f"{code}_polar_legacy.json"

        _write_pdb_from_smiles(smiles, code, ligand_pdb)
        _generate_cgmap(ligand_pdb, cgmap_json, cgmap_script=cgmap_script, vdm_db=vdm_db, repo_root=repo_root)
        _extract_new_polar_openbabel(repo_root, ligand_pdb, polar_new)
        _write_json(polar_legacy, _extract_legacy_openbabel(ligand_pdb))

        role_compare = _compare_role_payloads(_read_json(polar_legacy), _read_json(polar_new))
        _write_json(ldir / f"{code}_role_compare.json", role_compare)

        try:
            new_pipeline = _run_pipeline_with_polar(
                repo_root=repo_root,
                py311=py311,
                vdxform_dir=vdxform_dir,
                ligand_pdb=ligand_pdb,
                cgmap_json=cgmap_json,
                polar_json=polar_new,
                run_dir=ldir / "new_pipeline",
                ligand_code=code,
            )
        except Exception as e:
            new_pipeline = {"status": "error", "error": repr(e)}

        core_records.append(
            {
                "code": code,
                "smiles": smiles,
                "paths": {
                    "ligand_pdb": str(ligand_pdb),
                    "cgmap_json": str(cgmap_json),
                    "polar_new_json": str(polar_new),
                    "polar_legacy_json": str(polar_legacy),
                },
                "role_compare": role_compare,
                "new_pipeline": new_pipeline,
            }
        )

    # Dropped-cation probe: use one synthetic ligand with an N+ center expected to be dropped.
    probe_code = DROP_PROBE["code"]
    probe_smiles = DROP_PROBE["smiles"]
    pdir = lig_dir / probe_code
    pdir.mkdir(parents=True, exist_ok=True)
    probe_pdb = pdir / f"{probe_code}.pdb"
    probe_cgmap = pdir / f"{probe_code}_cgmap.json"
    probe_new = pdir / f"{probe_code}_polar_new.json"
    probe_legacy = pdir / f"{probe_code}_polar_legacy.json"

    _write_pdb_from_smiles(probe_smiles, probe_code, probe_pdb)
    _generate_cgmap(probe_pdb, probe_cgmap, cgmap_script=cgmap_script, vdm_db=vdm_db, repo_root=repo_root)
    _extract_new_polar_openbabel(repo_root, probe_pdb, probe_new)
    _write_json(probe_legacy, _extract_legacy_openbabel(probe_pdb))
    probe_cmp = _compare_role_payloads(_read_json(probe_legacy), _read_json(probe_new))

    probe_focus_atom: str | None = None
    focus_new_roles: list[str] = []
    focus_legacy_roles: list[str] = []
    legacy_focus_salt_atoms: list[str] = []
    new_focus_salt_atoms: list[str] = []
    dropped_atoms_with_legacy_saltbridge: list[str] = []
    probe_new_focus_run: dict[str, Any] | None = None
    probe_legacy_focus_run: dict[str, Any] | None = None

    if probe_cmp["dropped_cation_atoms"]:
        probe_focus_atom = str(probe_cmp["dropped_cation_atoms"][0])
        probe_new_single = pdir / f"{probe_code}_polar_new_{probe_focus_atom}.json"
        probe_legacy_single = pdir / f"{probe_code}_polar_legacy_{probe_focus_atom}.json"
        _make_single_atom_polar(probe_new, probe_focus_atom, probe_new_single)
        _make_single_atom_polar(probe_legacy, probe_focus_atom, probe_legacy_single)
        focus_new_roles = list(_read_json(probe_new_single)["sites"][0]["roles"])
        focus_legacy_roles = list(_read_json(probe_legacy_single)["sites"][0]["roles"])

        probe_new_focus_run = _run_pipeline_with_polar(
            repo_root=repo_root,
            py311=py311,
            vdxform_dir=vdxform_dir,
            ligand_pdb=probe_pdb,
            cgmap_json=probe_cgmap,
            polar_json=probe_new_single,
            run_dir=pdir / "focus_new_pipeline",
            ligand_code=probe_code,
        )
        probe_legacy_focus_run = _run_pipeline_with_polar(
            repo_root=repo_root,
            py311=py311,
            vdxform_dir=vdxform_dir,
            ligand_pdb=probe_pdb,
            cgmap_json=probe_cgmap,
            polar_json=probe_legacy_single,
            run_dir=pdir / "focus_legacy_pipeline",
            ligand_code=probe_code,
        )

        legacy_focus_salt_atoms = _saltbridge_atom_names_from_plip(
            Path(probe_legacy_focus_run["paths"]["plip_json"])
        )
        new_focus_salt_atoms = _saltbridge_atom_names_from_plip(Path(probe_new_focus_run["paths"]["plip_json"]))
        dropped_atoms_with_legacy_saltbridge = [a for a in probe_cmp["dropped_cation_atoms"] if a in legacy_focus_salt_atoms]

    n_plip_runs = sum(1 for r in core_records if r["new_pipeline"].get("status") == "ok")
    n_plip_all_sat = sum(
        1
        for r in core_records
        if r["new_pipeline"].get("status") == "ok"
        and bool((r["new_pipeline"].get("metrics") or {}).get("plip_all_satisfied", False))
    )
    total_dropped_core = sum(len(r["role_compare"]["dropped_cation_atoms"]) for r in core_records)

    summary: dict[str, Any] = {
        "tag": args.tag,
        "outdir": str(outdir),
        "config": {
            "repo_root": str(repo_root),
            "python_current": sys.executable,
            "python311": str(py311),
            "vdxform_dir": str(vdxform_dir),
            "cgmap_script": str(cgmap_script),
            "vdm_db": str(vdm_db),
        },
        "core_ligands": core_records,
        "dropped_cation_probe": {
            "code": probe_code,
            "smiles": probe_smiles,
            "paths": {
                "ligand_pdb": str(probe_pdb),
                "cgmap_json": str(probe_cgmap),
                "polar_new_json": str(probe_new),
                "polar_legacy_json": str(probe_legacy),
            },
            "role_compare": probe_cmp,
            "focus_atom": probe_focus_atom,
            "focus_new_roles": focus_new_roles,
            "focus_legacy_roles": focus_legacy_roles,
            "focus_new_pipeline": probe_new_focus_run,
            "focus_legacy_pipeline": probe_legacy_focus_run,
            "new_focus_saltbridge_atoms": new_focus_salt_atoms,
            "legacy_focus_saltbridge_atoms": legacy_focus_salt_atoms,
            "dropped_atoms_with_legacy_saltbridge": dropped_atoms_with_legacy_saltbridge,
        },
        "acceptance_check": {
            "plip_runs_on_core_non_mtx": n_plip_runs,
            "plip_all_satisfied_on_core_non_mtx": n_plip_all_sat,
            "has_minimum_three_plip_runs": bool(n_plip_runs >= 3),
            "total_dropped_cation_atoms_across_core": total_dropped_core,
            "dropped_cation_has_true_saltbridge_evidence": bool(len(dropped_atoms_with_legacy_saltbridge) > 0),
        },
    }

    out_json = outdir / "summary.json"
    out_md = outdir / "summary.md"
    _write_json(out_json, summary)
    _write_summary_md(summary, out_md)
    print(str(out_json))
    print(str(out_md))


if __name__ == "__main__":
    main()
