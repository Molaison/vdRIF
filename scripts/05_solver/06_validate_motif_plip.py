#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    atom_name: str
    resname: str
    chain: str
    resnum: int


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_pdb_atoms(pdb_path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        # PDB fixed-width columns
        serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resnum = int(line[22:26].strip())
        atoms.append(PdbAtom(serial=serial, atom_name=atom_name, resname=resname, chain=chain, resnum=resnum))
    return atoms


def _choose_ligand_residue(motif_pdb: Path, polar_atom_names: set[str]) -> tuple[str, str, int]:
    atoms = _parse_pdb_atoms(motif_pdb)
    # Prefer the residue that overlaps most with the polar atom names.
    by_res: dict[tuple[str, str, int], set[str]] = {}
    for a in atoms:
        if motif_pdb.name and a.resname == "HOH":
            continue
        by_res.setdefault((a.resname, a.chain, a.resnum), set()).add(a.atom_name)

    best_key: tuple[str, str, int] | None = None
    best_overlap = -1
    for key, anames in by_res.items():
        overlap = len(anames & polar_atom_names)
        if overlap > best_overlap:
            best_overlap = overlap
            best_key = key
    if best_key is None or best_overlap <= 0:
        # Fallback: first HETATM residue.
        for line in motif_pdb.read_text(encoding="utf-8").splitlines():
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resnum = int(line[22:26].strip())
                return resname, chain, resnum
        raise ValueError(f"Failed to identify ligand residue in motif PDB: {motif_pdb}")

    return best_key


def _find_latest(glob_parent: Path, pattern: str) -> Path:
    matches = sorted(glob_parent.glob(pattern), key=lambda p: p.stat().st_mtime)
    if not matches:
        raise FileNotFoundError(f"Expected at least 1 match for {pattern} in {glob_parent}, got 0")
    return matches[-1]


def _parse_int_list(s: str) -> list[int]:
    # Examples seen in PLIP XML: "27,28" or "27, 28"
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in re.split(r"\s*,\s*", s) if x]


def _read_plip_idx_list(parent: ET.Element, tag: str) -> list[int]:
    """
    PLIP XML uses two formats depending on interaction type/version:
    - <lig_idx_list>27,28</lig_idx_list>
    - <lig_idx_list><idx id="1">27</idx><idx id="2">28</idx></lig_idx_list>
    """
    el = parent.find(tag)
    if el is None:
        return []
    # Nested <idx> children
    idx_children = el.findall("idx")
    if idx_children:
        out: list[int] = []
        for c in idx_children:
            if c.text is None:
                continue
            try:
                out.append(int(c.text.strip()))
            except ValueError:
                continue
        return out
    # Flat text list
    if el.text:
        return _parse_int_list(el.text)
    return []


def _collect_ligand_atom_serials_from_plip(plipfixed_pdb: Path, ligand_key: tuple[str, str, int]) -> dict[int, str]:
    resname, chain, resnum = ligand_key
    out: dict[int, str] = {}
    for a in _parse_pdb_atoms(plipfixed_pdb):
        if a.resname == resname and a.chain == chain and a.resnum == resnum:
            out[a.serial] = a.atom_name
    if not out:
        raise ValueError(f"PLIP fixed PDB does not contain ligand residue {ligand_key}: {plipfixed_pdb}")
    return out


def _collect_interacting_ligand_serials(report_xml: Path, ligand_serials: set[int]) -> set[int]:
    tree = ET.parse(str(report_xml))
    root = tree.getroot()
    hit: set[int] = set()

    # Hydrogen bonds: donoridx / acceptoridx
    for hb in root.findall(".//hydrogen_bond"):
        d = hb.findtext("donoridx")
        a = hb.findtext("acceptoridx")
        for x in (d, a):
            if x is None:
                continue
            try:
                idx = int(x.strip())
            except ValueError:
                continue
            if idx in ligand_serials:
                hit.add(idx)

    # Salt bridges: lig_idx_list
    for sb in root.findall(".//salt_bridge"):
        for idx in _read_plip_idx_list(sb, "lig_idx_list"):
            if idx in ligand_serials:
                hit.add(idx)

    # Water bridges (optional): lig_idx_list (PLIP 3.x uses this naming)
    for wb in root.findall(".//water_bridge"):
        for idx in _read_plip_idx_list(wb, "lig_idx_list"):
            if idx in ligand_serials:
                hit.add(idx)

    return hit


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate ligand polar-atom satisfaction using PLIP.")
    ap.add_argument("--motif-pdb", type=Path, required=True, help="Complex PDB containing ligand + motif residues.")
    ap.add_argument("--polar-sites", type=Path, required=True, help="Polar atom list to be satisfied (AtomSite JSON).")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--maxthreads", type=int, default=1)
    ap.add_argument("--timeout-s", type=float, default=120.0)
    ap.add_argument("--keep-plip-outdir", action="store_true")
    ap.add_argument("--plip-outdir", type=Path, default=None, help="Optional explicit PLIP output dir (for debugging).")
    args = ap.parse_args()

    polar = _read_json(args.polar_sites)
    polar_atom_names = {str(s["atom_name"]) for s in polar.get("sites", [])}
    if not polar_atom_names:
        raise ValueError(f"No polar sites in {args.polar_sites}")

    ligand_key = _choose_ligand_residue(args.motif_pdb, polar_atom_names)

    outdir = args.plip_outdir
    if outdir is None:
        outdir = args.out.parent / f"plip_{args.motif_pdb.stem}"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "plip",
        "-f",
        str(args.motif_pdb),
        "-x",
        "-t",
        "--maxthreads",
        str(int(args.maxthreads)),
        "-o",
        str(outdir),
        "--silent",
    ]
    subprocess.run(cmd, check=True, timeout=float(args.timeout_s))

    report_xml = outdir / f"{args.motif_pdb.stem}_report.xml"
    if not report_xml.exists():
        # Some PLIP versions rename output based on input basename, keep a robust glob fallback.
        report_xml = _find_latest(outdir, "*_report.xml")

    plipfixed = _find_latest(outdir, "plipfixed.*.pdb")

    ligand_serial_to_name = _collect_ligand_atom_serials_from_plip(plipfixed, ligand_key)
    interacting_serials = _collect_interacting_ligand_serials(report_xml, set(ligand_serial_to_name.keys()))

    satisfied_by_name = {ligand_serial_to_name[i] for i in interacting_serials}
    # Restrict to the polar atom list (PLIP may also report interactions for non-polar atoms).
    satisfied_polar = sorted(satisfied_by_name & polar_atom_names)
    unsatisfied_polar = sorted(polar_atom_names - set(satisfied_polar))

    payload = {
        "all_satisfied": len(unsatisfied_polar) == 0,
        "inputs": {"motif_pdb": str(args.motif_pdb), "polar_sites": str(args.polar_sites)},
        "ligand": {"resname": ligand_key[0], "chain": ligand_key[1], "resnum": ligand_key[2]},
        "n_polar_atoms": len(polar_atom_names),
        "n_satisfied": len(satisfied_polar),
        "satisfied_polar_atoms": satisfied_polar,
        "unsatisfied_polar_atoms": unsatisfied_polar,
        "plip": {"report_xml": str(report_xml), "plipfixed_pdb": str(plipfixed), "outdir": str(outdir)},
    }
    _write_json(args.out, payload)

    if not args.keep_plip_outdir:
        # Keep only the JSON by default; PLIP output can be large / noisy.
        # Users can re-run with --keep-plip-outdir if they need details.
        for p in outdir.iterdir():
            if p == report_xml or p == plipfixed:
                continue
            # leave report_xml + plipfixed for minimal debugging
            try:
                p.unlink()
            except IsADirectoryError:
                pass


if __name__ == "__main__":
    main()
