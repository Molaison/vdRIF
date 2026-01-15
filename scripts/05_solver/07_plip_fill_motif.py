#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
import importlib.util
import sys


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_text(path: Path, s: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s, encoding="utf-8")


def _atom_name(atom: Chem.Atom) -> str:
    info = atom.GetPDBResidueInfo()
    if info is None:
        return f"ATOM{atom.GetIdx()}"
    return info.GetName().strip()


def _load_pdb_mol(pdb_path: Path) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"Failed to read PDB: {pdb_path}")
    if mol.GetNumConformers() != 1:
        raise ValueError(f"Expected 1 conformer in {pdb_path}, got {mol.GetNumConformers()}")
    return mol


def _ligand_donor_h_xyz_by_atom_name(mol: Chem.Mol) -> dict[str, np.ndarray]:
    conf = mol.GetConformer()
    out: dict[str, list[list[float]]] = {}
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        hs = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "H"]
        if not hs:
            continue
        an = _atom_name(a)
        coords: list[list[float]] = []
        for h in hs:
            p = conf.GetAtomPosition(h.GetIdx())
            coords.append([p.x, p.y, p.z])
        out[an] = coords
    return {k: np.asarray(v, dtype=np.float64) for k, v in out.items()}


@dataclass(frozen=True)
class PdbAtom:
    record: str  # ATOM/HETATM
    serial: int
    atom_name: str
    resname: str
    chain: str
    resnum: int
    x: float
    y: float
    z: float
    element: str

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=np.float64)


def _parse_pdb_atoms(pdb_path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        record = line[0:6].strip()
        serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resnum = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        element = (line[76:78].strip() or atom_name[0]).strip()
        atoms.append(
            PdbAtom(
                record=record,
                serial=serial,
                atom_name=atom_name,
                resname=resname,
                chain=chain,
                resnum=resnum,
                x=x,
                y=y,
                z=z,
                element=element,
            )
        )
    return atoms


def _format_pdb_atom(a: PdbAtom) -> str:
    # Minimal but column-correct PDB ATOM/HETATM writer.
    # Columns (1-based):
    # 1-6 record, 7-11 serial, 13-16 atom, 17 altLoc, 18-20 resName, 22 chain,
    # 23-26 resSeq, 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 temp, 77-78 element
    altloc = " "
    icode = " "
    occ = 1.00
    temp = 0.00
    return (
        f"{a.record:<6}{a.serial:>5} {a.atom_name:>4}{altloc}{a.resname:>3} {a.chain:1}"
        f"{a.resnum:>4}{icode}   {a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}{occ:>6.2f}{temp:>6.2f}          {a.element:>2}"
    )


def _vdw_radius(element: str) -> float:
    vdw_by_elem = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
    return float(vdw_by_elem.get(element, 1.70))


def _unpack_xform12(x12: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Packed row-major 3x4 (R|t) with translation interleaved per-row:
    # [R00,R01,R02,t0, R10,R11,R12,t1, R20,R21,R22,t2]
    x12 = np.asarray(x12, dtype=np.float64).reshape(12)
    R = np.array(
        [
            [x12[0], x12[1], x12[2]],
            [x12[4], x12[5], x12[6]],
            [x12[8], x12[9], x12[10]],
        ],
        dtype=np.float64,
    )
    t = np.array([x12[3], x12[7], x12[11]], dtype=np.float64)
    return R, t


def _candidate_atoms_world(
    atom_order: list[str],
    aa3: str,
    chain: str,
    resnum: int,
    xform_world_stub_12: np.ndarray,
    center_xyz_stub: np.ndarray,
) -> list[PdbAtom]:
    R, t = _unpack_xform12(xform_world_stub_12)
    xyz_stub = center_xyz_stub.astype(np.float64)
    xyz_world = (xyz_stub @ R.T) + t[None, :]

    out: list[PdbAtom] = []
    for atom_name, (x, y, z) in zip(atom_order, xyz_world, strict=True):
        if not np.isfinite(x) or not np.isfinite(y) or not np.isfinite(z):
            continue
        element = "S" if atom_name.startswith("S") else atom_name[0]
        out.append(
            PdbAtom(
                record="ATOM",
                serial=0,
                atom_name=str(atom_name),
                resname=str(aa3),
                chain=str(chain),
                resnum=int(resnum),
                x=float(x),
                y=float(y),
                z=float(z),
                element=element,
            )
        )
    return out


def _clashes(
    atoms_a: list[PdbAtom],
    atoms_b: list[PdbAtom],
    tolerance: float,
    heavy_only: bool = True,
) -> bool:
    if heavy_only:
        atoms_a = [a for a in atoms_a if a.element != "H"]
        atoms_b = [b for b in atoms_b if b.element != "H"]
    if not atoms_a or not atoms_b:
        return False
    xa = np.stack([a.xyz for a in atoms_a], axis=0)
    xb = np.stack([b.xyz for b in atoms_b], axis=0)
    ra = np.array([_vdw_radius(a.element) for a in atoms_a], dtype=np.float64)
    rb = np.array([_vdw_radius(b.element) for b in atoms_b], dtype=np.float64)
    thr = (ra[:, None] + rb[None, :] - float(tolerance)) ** 2
    d = xa[:, None, :] - xb[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", d, d)
    return bool((d2 < thr).any())


def _score_candidate_for_target(
    atom_order: list[str],
    center_xyz_stub: np.ndarray,
    xform_world_stub_12: np.ndarray,
    lig_atom_name: str,
    lig_xyz: np.ndarray,
    lig_roles: set[str],
    lig_donor_h_xyz_by_name: dict[str, np.ndarray],
) -> float:
    """
    Heuristic ranking score for a candidate vs a single ligand atom:
    - For ligand acceptors: favor close protein donors with good B-D-A angle.
    - For ligand donors: favor close protein acceptors aligned with ligand H.
    Higher is better.
    """
    # donor base mapping (same as in candidate generation)
    donor_atoms = {"N", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "ND1", "ND2", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD", "SG"}
    donor_base_name = {
        "N": "CA",
        "NZ": "CE",
        "NE": "CD",
        "NH1": "CZ",
        "NH2": "CZ",
        "ND2": "CG",
        "NE2": "CD",
        "ND1": "CG",
        "NE1": "CD1",
        "OG": "CB",
        "OG1": "CB",
        "OH": "CZ",
    }

    # Build world coords for all atoms once.
    R, t = _unpack_xform12(xform_world_stub_12)
    xyz_world = (center_xyz_stub.astype(np.float64) @ R.T) + t[None, :]
    name_to_xyz = {nm: xyz_world[i] for i, nm in enumerate(atom_order) if np.all(np.isfinite(xyz_world[i]))}

    best = -1e9

    if "acceptor" in lig_roles:
        # Need a protein donor.
        for dnm in donor_atoms:
            dxyz = name_to_xyz.get(dnm)
            bnm = donor_base_name.get(dnm)
            bxyz = name_to_xyz.get(bnm) if bnm else None
            if dxyz is None or bxyz is None:
                continue
            da = lig_xyz - dxyz
            dist = float(np.linalg.norm(da))
            if dist > 3.6:
                continue
            v_db = bxyz - dxyz
            v_da = lig_xyz - dxyz
            n1 = float(np.linalg.norm(v_db))
            n2 = float(np.linalg.norm(v_da))
            if n1 <= 0 or n2 <= 0:
                continue
            cos = float(np.clip(np.dot(v_db, v_da) / (n1 * n2), -1.0, 1.0))
            ang = math.degrees(math.acos(cos))
            if ang < 110.0:
                continue
            # Score: closer + straighter
            best = max(best, 5.0 - dist + (ang - 110.0) / 100.0)

    if "donor" in lig_roles:
        # Need a protein acceptor aligned with ligand H.
        hs = lig_donor_h_xyz_by_name.get(lig_atom_name)
        if hs is None or hs.size == 0:
            return best
        for anm in acceptor_atoms:
            axyz = name_to_xyz.get(anm)
            if axyz is None:
                continue
            da = axyz - lig_xyz
            dist = float(np.linalg.norm(da))
            if dist > 3.6:
                continue
            for hxyz in hs:
                dh = hxyz - lig_xyz  # D->H
                n1 = float(np.linalg.norm(dh))
                n2 = float(np.linalg.norm(da))
                if n1 <= 0 or n2 <= 0:
                    continue
                cos = float(np.clip(np.dot(dh, da) / (n1 * n2), -1.0, 1.0))
                ang = math.degrees(math.acos(cos))
                ha = float(np.linalg.norm(axyz - hxyz))
                if ang < 90.0 or ha > 3.2:
                    continue
                best = max(best, 5.0 - dist + (ang - 90.0) / 120.0)

    return best


def main() -> None:
    ap = argparse.ArgumentParser(description="Iteratively add candidates until PLIP satisfies all ligand polar atoms.")
    ap.add_argument("--motif-pdb", type=Path, required=True)
    ap.add_argument("--candidates-npz", type=Path, required=True)
    ap.add_argument("--candidates-meta", type=Path, required=True)
    ap.add_argument("--polar-sites", type=Path, required=True)
    ap.add_argument("--ligand-pdb", type=Path, required=True)
    ap.add_argument("--out-pdb", type=Path, required=True)
    ap.add_argument("--out-report-json", type=Path, required=True)
    ap.add_argument("--max-res", type=int, default=15)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument("--top-try-per-atom", type=int, default=80)
    args = ap.parse_args()

    polar = _read_json(args.polar_sites)
    polar_by_name = {str(s["atom_name"]): set(s["roles"]) for s in polar.get("sites", [])}
    polar_xyz_by_name = {str(s["atom_name"]): np.asarray(s["xyz"], dtype=np.float64) for s in polar.get("sites", [])}
    polar_atom_names = sorted(polar_by_name.keys())

    cand_meta = _read_json(args.candidates_meta)
    atom_order = list(cand_meta.get("atom_order") or [])
    polar_atoms_order = list(cand_meta.get("polar_atoms") or [])
    if not atom_order:
        raise ValueError(f"Missing atom_order in candidates meta: {args.candidates_meta}")
    if not polar_atoms_order:
        raise ValueError(f"Missing polar_atoms in candidates meta: {args.candidates_meta}")
    polar_to_bit = {nm: i for i, nm in enumerate(polar_atoms_order)}

    z = np.load(args.candidates_npz)
    aa3 = z["aa3"].astype("U3")
    vdm_id = z["vdm_id_u64"].astype(np.uint64)
    score = z["score_f32"].astype(np.float64)
    cover = z["cover_mask_u16"].astype(np.uint16)
    xw12 = z["xform_world_stub_12_f32"].astype(np.float64)
    center_stub = z["center_atom_xyz_stub_f32"].astype(np.float64)

    ligand_mol = _load_pdb_mol(args.ligand_pdb)
    lig_donor_h_xyz = _ligand_donor_h_xyz_by_atom_name(ligand_mol)

    # Start from the existing motif PDB (ligand + current residues).
    base_atoms = _parse_pdb_atoms(args.motif_pdb)
    # Identify current maximum resnum on chain M.
    max_resnum_m = max((a.resnum for a in base_atoms if a.chain == "M"), default=0)
    current_res_count = len({(a.chain, a.resnum) for a in base_atoms if a.chain == "M"})

    # Separate ligand vs protein atoms for clash checks.
    # Ligand: best overlap with polar names.
    by_res: dict[tuple[str, str, int], set[str]] = {}
    for a in base_atoms:
        by_res.setdefault((a.resname, a.chain, a.resnum), set()).add(a.atom_name)
    ligand_key = max(by_res.items(), key=lambda kv: len(kv[1] & set(polar_atom_names)))[0]
    ligand_atoms = [a for a in base_atoms if (a.resname, a.chain, a.resnum) == ligand_key]
    protein_atoms = [a for a in base_atoms if a.record == "ATOM" and a.chain == "M"]

    report: dict[str, Any] = {
        "inputs": {
            "motif_pdb": str(args.motif_pdb),
            "candidates_npz": str(args.candidates_npz),
            "candidates_meta": str(args.candidates_meta),
            "polar_sites": str(args.polar_sites),
        },
        "steps": [],
    }

    # Iterate until PLIP says all satisfied.
    validator_path = Path(__file__).with_name("06_validate_motif_plip.py")
    spec = importlib.util.spec_from_file_location("_plip_validate", validator_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load PLIP validator module: {validator_path}")
    _plip_mod = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = _plip_mod
    spec.loader.exec_module(_plip_mod)  # type: ignore[attr-defined]
    _plip_validate_main = getattr(_plip_mod, "main")

    def plip_validate(tmp_pdb: Path) -> dict[str, Any]:
        out_json = tmp_pdb.with_suffix(".plip.json")
        # Call the validator as a module entrypoint by emulating argv.
        import sys

        argv0 = sys.argv[:]
        try:
            sys.argv = [
                "06_validate_motif_plip.py",
                "--motif-pdb",
                str(tmp_pdb),
                "--polar-sites",
                str(args.polar_sites),
                "-o",
                str(out_json),
                "--maxthreads",
                "1",
                "--timeout-s",
                "120",
            ]
            _plip_validate_main()
        finally:
            sys.argv = argv0
        return _read_json(out_json)

    with tempfile.TemporaryDirectory(prefix="plip_fill_") as td:
        tmp_dir = Path(td)
        tmp_pdb = tmp_dir / "current.pdb"

        # Start with the current motif.
        cur_atoms = list(base_atoms)
        for i, a in enumerate(cur_atoms, start=1):
            cur_atoms[i - 1] = PdbAtom(**{**a.__dict__, "serial": i})
        _write_text(tmp_pdb, "\n".join(_format_pdb_atom(a) for a in cur_atoms) + "\nEND\n")

        val = plip_validate(tmp_pdb)
        missing = list(val.get("unsatisfied_polar_atoms") or [])

        while missing:
            if current_res_count >= int(args.max_res):
                break

            target = str(sorted(missing)[0])
            bit = polar_to_bit.get(target)
            if bit is None:
                raise ValueError(f"Target polar atom {target} not found in candidates meta polar_atoms")

            # Candidate indices that claim to cover this target.
            idxs = np.nonzero((cover & (np.uint16(1 << bit))) != 0)[0]
            if idxs.size == 0:
                report["steps"].append({"target": target, "status": "no_candidates"})
                break

            # Rank candidates deterministically with a geometry heuristic (higher is better),
            # tie-break by score desc then vdm_id asc then cand index.
            scored: list[tuple[float, float, int, int]] = []
            lig_xyz = polar_xyz_by_name[target]
            lig_roles = polar_by_name[target]
            for i in idxs:
                geom = _score_candidate_for_target(
                    atom_order=atom_order,
                    center_xyz_stub=center_stub[i],
                    xform_world_stub_12=xw12[i],
                    lig_atom_name=target,
                    lig_xyz=lig_xyz,
                    lig_roles=lig_roles,
                    lig_donor_h_xyz_by_name=lig_donor_h_xyz,
                )
                scored.append((float(geom), float(score[i]), int(vdm_id[i]), int(i)))
            scored.sort(key=lambda x: (-x[0], -x[1], x[2], x[3]))

            chosen_idx: int | None = None
            for geom, sc, v_id, i in scored[: int(args.top_try_per_atom)]:
                new_resnum = max_resnum_m + 1
                cand_atoms = _candidate_atoms_world(
                    atom_order=atom_order,
                    aa3=str(aa3[i]),
                    chain="M",
                    resnum=new_resnum,
                    xform_world_stub_12=xw12[i],
                    center_xyz_stub=center_stub[i],
                )
                # Quick clash checks vs ligand and existing motif protein.
                if _clashes(cand_atoms, ligand_atoms, tolerance=float(args.clash_tol), heavy_only=True):
                    continue
                if _clashes(cand_atoms, protein_atoms, tolerance=float(args.clash_tol), heavy_only=True):
                    continue

                # Write a temporary PDB and ask PLIP whether this improves coverage.
                test_atoms = list(cur_atoms) + cand_atoms
                # Renumber serials deterministically
                test_atoms2: list[PdbAtom] = []
                for srl, a in enumerate(test_atoms, start=1):
                    test_atoms2.append(PdbAtom(**{**a.__dict__, "serial": srl}))
                _write_text(tmp_pdb, "\n".join(_format_pdb_atom(a) for a in test_atoms2) + "\nEND\n")

                test_val = plip_validate(tmp_pdb)
                test_missing = list(test_val.get("unsatisfied_polar_atoms") or [])
                if len(test_missing) < len(missing) and target not in test_missing:
                    chosen_idx = int(i)
                    # Accept
                    cur_atoms = test_atoms2
                    protein_atoms.extend(cand_atoms)
                    max_resnum_m = new_resnum
                    current_res_count += 1
                    missing = test_missing
                    report["steps"].append(
                        {
                            "target": target,
                            "chosen_candidate_index": chosen_idx,
                            "aa3": str(aa3[i]),
                            "vdm_id_u64": int(v_id),
                            "score": float(sc),
                            "geom_score": float(geom),
                            "remaining_unsatisfied": list(missing),
                        }
                    )
                    break

            if chosen_idx is None:
                report["steps"].append({"target": target, "status": "no_improving_candidate", "remaining": list(missing)})
                break

        # Write final
        _write_text(args.out_pdb, "\n".join(_format_pdb_atom(a) for a in cur_atoms) + "\nEND\n")
        report["final_plip"] = plip_validate(args.out_pdb)

    args.out_report_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_report_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
