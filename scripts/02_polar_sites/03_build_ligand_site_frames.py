#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _load_pdb_mol(pdb_path: Path) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"Failed to read PDB: {pdb_path}")
    if mol.GetNumConformers() != 1:
        raise ValueError(f"Expected 1 conformer in {pdb_path}, got {mol.GetNumConformers()}")
    return mol


def _atom_name(atom: Chem.Atom) -> str:
    info = atom.GetPDBResidueInfo()
    if info is None:
        return f"ATOM{atom.GetIdx()}"
    return info.GetName().strip()


def _name_to_idx(mol: Chem.Mol) -> dict[str, int]:
    d: dict[str, int] = {}
    for a in mol.GetAtoms():
        name = _atom_name(a)
        # Ligand atom names are assumed unique within the ligand residue.
        d[name] = int(a.GetIdx())
    return d


def _xyz(mol: Chem.Mol, atom_idx: int) -> np.ndarray:
    p = mol.GetConformer().GetAtomPosition(atom_idx)
    return np.array([p.x, p.y, p.z], dtype=np.float64)


def _pick_frame_indices(correspond_names: list[str]) -> tuple[int, int, int]:
    non_h = [i for i, n in enumerate(correspond_names) if not n.startswith("H")]
    if len(non_h) >= 3:
        return (non_h[0], non_h[1], non_h[2])
    if len(correspond_names) < 3:
        raise ValueError(f"Need >=3 correspond_names to define a frame, got: {correspond_names}")
    return (0, 1, 2)


def _pick_frame_indices_from_defs(
    correspond_names: list[str], cg: str, cg_frame_defs: dict[str, list[list[str]]] | None
) -> tuple[int, int, int] | None:
    if not cg_frame_defs:
        return None
    alts = cg_frame_defs.get(cg)
    if not alts:
        return None
    for want in alts:
        if len(want) != 3:
            continue
        idxs: list[int] = []
        ok = True
        for atom_name in want:
            try:
                idxs.append(correspond_names.index(atom_name))
            except ValueError:
                ok = False
                break
        if ok:
            return (int(idxs[0]), int(idxs[1]), int(idxs[2]))
    return None


def _frame_from_points(a0: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Deterministic right-handed frame.

    origin: a1
    x-axis: a1->a2
    y-axis: component of a1->a0 orthogonal to x
    z-axis: x cross y

    Returns (R, t) where X_world(p_local) = R @ p_local + t.
    """
    x = a2 - a1
    nx = np.linalg.norm(x)
    if nx < 1e-8:
        raise ValueError("Degenerate frame: a1 and a2 coincide")
    x = x / nx

    v = a0 - a1
    v = v - float(np.dot(v, x)) * x
    nv = np.linalg.norm(v)
    if nv < 1e-8:
        raise ValueError("Degenerate frame: a0 is colinear with a1->a2")
    y = v / nv

    z = np.cross(x, y)
    nz = np.linalg.norm(z)
    if nz < 1e-8:
        raise ValueError("Degenerate frame: x and y not linearly independent")
    z = z / nz

    R = np.stack([x, y, z], axis=1)  # columns are basis vectors
    t = a1
    return R, t


def _neighbors_by_name(mol: Chem.Mol, name_to_idx: dict[str, int], atom_name: str) -> list[str]:
    idx = name_to_idx[atom_name]
    a = mol.GetAtomWithIdx(idx)
    names = []
    for nb in a.GetNeighbors():
        names.append(_atom_name(nb))
    return names


def _attached_h_names(mol: Chem.Mol, name_to_idx: dict[str, int], atom_name: str) -> list[str]:
    idx = name_to_idx[atom_name]
    a = mol.GetAtomWithIdx(idx)
    hs = []
    for nb in a.GetNeighbors():
        if nb.GetSymbol() == "H":
            hs.append(_atom_name(nb))
    return sorted(hs)


def _heavy_neighbor_name(mol: Chem.Mol, name_to_idx: dict[str, int], atom_name: str) -> str:
    idx = name_to_idx[atom_name]
    a = mol.GetAtomWithIdx(idx)
    heavy = [nb for nb in a.GetNeighbors() if nb.GetSymbol() != "H"]
    if len(heavy) != 1:
        raise ValueError(f"Expected exactly 1 heavy neighbor for {atom_name}, got {len(heavy)}")
    return _atom_name(heavy[0])


def _pick_two_heavy_neighbors_for_ccn_frame(mol: Chem.Mol, name_to_idx: dict[str, int], atom_name: str) -> tuple[str, str]:
    """
    Build a lysine-like (CD, CE, NZ) local frame for a cationic atom in an arbitrary ligand.

    Returns (a0, a1) heavy neighbor names to be used as (CD, CE) with a2 being the cation atom (NZ).
    Deterministic tie-break by atom name.
    """
    idx = name_to_idx[atom_name]
    a = mol.GetAtomWithIdx(idx)
    heavy = sorted([nb for nb in a.GetNeighbors() if nb.GetSymbol() != "H"], key=lambda x: _atom_name(x))
    if len(heavy) >= 2:
        # Use the two nearest heavy neighbors.
        a1 = _atom_name(heavy[0])
        a0 = _atom_name(heavy[1])
        return a0, a1
    if len(heavy) == 1:
        a1_atom = heavy[0]
        # Find a second heavy atom by stepping one bond further (exclude the cation itself).
        heavy2 = sorted(
            [nb for nb in a1_atom.GetNeighbors() if nb.GetSymbol() != "H" and int(nb.GetIdx()) != idx],
            key=lambda x: _atom_name(x),
        )
        if not heavy2:
            raise ValueError(f"Cannot build ccn frame for {atom_name}: only one heavy neighbor and no second-step heavy atom.")
        a1 = _atom_name(a1_atom)
        a0 = _atom_name(heavy2[0])
        return a0, a1
    raise ValueError(f"Cannot build ccn frame for {atom_name}: no heavy neighbors.")


@dataclass(frozen=True)
class SiteFrame:
    site_id: str
    cg: str
    lgd_sel: tuple[str, ...]
    correspond_resname: str
    correspond_names: tuple[str, ...]
    frame_atom_names: tuple[str, str, str]  # (a0,a1,a2) from lgd_sel
    R: np.ndarray  # 3x3
    t: np.ndarray  # 3
    covers_polar_atoms: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "site_id": self.site_id,
            "cg": self.cg,
            "lgd_sel": list(self.lgd_sel),
            "correspond_resname": self.correspond_resname,
            "correspond_names": list(self.correspond_names),
            "frame_atom_names": list(self.frame_atom_names),
            "R": self.R.tolist(),
            "t": self.t.tolist(),
            "covers_polar_atoms": list(self.covers_polar_atoms),
        }


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Build deterministic ligand site frames from cg_atommap, and optionally "
            "add fallback bb_cnh frames for uncovered donor atoms."
        )
    )
    ap.add_argument("ligand_pdb", type=Path)
    ap.add_argument("--cg-atommap", type=Path, required=True)
    ap.add_argument("--polar-sites", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument(
        "--cg-frame-defs",
        type=Path,
        default=Path("configs/cg_frame_defs.json"),
        help="JSON mapping cg_type -> [a0,a1,a2] residue-atom names used to define the iFG frame.",
    )
    ap.add_argument("--add-bb-cnh-for-uncovered-donors", action="store_true")
    ap.add_argument("--add-ccn-for-cations", action="store_true")
    args = ap.parse_args()

    mol = _load_pdb_mol(args.ligand_pdb)
    name_to_idx = _name_to_idx(mol)

    cgmap = _load_json(args.cg_atommap)
    polar = _load_json(args.polar_sites)
    polar_sites = {s["atom_name"]: set(s["roles"]) for s in polar["sites"]}
    polar_atom_names = set(polar_sites.keys())

    if not isinstance(cgmap, dict):
        raise TypeError("Expected cg_atommap JSON to be a dict keyed by 'i_j'.")

    cg_frame_defs: dict[str, list[list[str]]] | None = None
    if args.cg_frame_defs is not None and args.cg_frame_defs.exists():
        cg_frame_defs = _load_json(args.cg_frame_defs)

    sites: list[SiteFrame] = []
    covered_by_sites: set[str] = set()

    # 1) Add frames from cg_atommap entries
    for key in sorted(cgmap.keys()):
        entry = cgmap[key]
        lgd_sel = list(entry["lgd_sel"])
        corr_names = list(entry["correspond_names"])
        if len(lgd_sel) != len(corr_names):
            raise ValueError(f"Length mismatch in {key}: lgd_sel vs correspond_names")

        cg = str(entry["cg"])
        idxs = _pick_frame_indices_from_defs(corr_names, cg=cg, cg_frame_defs=cg_frame_defs)
        if idxs is None:
            idxs = _pick_frame_indices(corr_names)
        i0, i1, i2 = idxs
        a0n, a1n, a2n = lgd_sel[i0], lgd_sel[i1], lgd_sel[i2]
        c2n = corr_names[i2]
        for nm in (a0n, a1n, a2n):
            if nm not in name_to_idx:
                raise KeyError(f"Ligand atom name not found in PDB: {nm}")

        R, t = _frame_from_points(_xyz(mol, name_to_idx[a0n]), _xyz(mol, name_to_idx[a1n]), _xyz(mol, name_to_idx[a2n]))
        covers = sorted(set(lgd_sel) & polar_atom_names)
        covered_by_sites.update(covers)

        sites.append(
            SiteFrame(
                site_id=f"cgmap:{key}",
                cg=cg,
                lgd_sel=tuple(lgd_sel),
                correspond_resname=str(entry["correspond_resname"]),
                correspond_names=tuple(corr_names),
                frame_atom_names=(a0n, a1n, a2n),
                R=R,
                t=t,
                covers_polar_atoms=tuple(covers),
            )
        )

        # Symmetry handling (minimal MVP):
        # Some cg types have symmetric atoms (e.g., coo: OD1/OD2 and OE1/OE2).
        # In proteins these labels are well-defined, but in ligands O1/O2 style names are arbitrary.
        # To avoid mirroring artifacts and missed/misplaced interactions, emit a swapped frame too.
        if cg == "coo":
            alt_c2 = None
            if c2n == "OD1" and "OD2" in corr_names:
                alt_c2 = "OD2"
            elif c2n == "OD2" and "OD1" in corr_names:
                alt_c2 = "OD1"
            elif c2n == "OE1" and "OE2" in corr_names:
                alt_c2 = "OE2"
            elif c2n == "OE2" and "OE1" in corr_names:
                alt_c2 = "OE1"

            if alt_c2 is not None:
                j2 = corr_names.index(alt_c2)
                alt_a2n = lgd_sel[j2]
                if alt_a2n != a2n:
                    if alt_a2n not in name_to_idx:
                        raise KeyError(f"Ligand atom name not found in PDB (symmetry swap): {alt_a2n}")
                    R2, t2 = _frame_from_points(
                        _xyz(mol, name_to_idx[a0n]),
                        _xyz(mol, name_to_idx[a1n]),
                        _xyz(mol, name_to_idx[alt_a2n]),
                    )
                    sites.append(
                        SiteFrame(
                            site_id=f"cgmap:{key}:swap_{c2n}_to_{alt_c2}",
                            cg=cg,
                            lgd_sel=tuple(lgd_sel),
                            correspond_resname=str(entry["correspond_resname"]),
                            correspond_names=tuple(corr_names),
                            frame_atom_names=(a0n, a1n, alt_a2n),
                            R=R2,
                            t=t2,
                            covers_polar_atoms=tuple(covers),
                        )
                    )

    # 2) Optional fallback for uncovered donor atoms: map to bb_cnh (CA,N,H)
    if args.add_bb_cnh_for_uncovered_donors:
        uncovered = sorted(polar_atom_names - covered_by_sites)
        for atom_name in uncovered:
            roles = polar_sites.get(atom_name, set())
            if "donor" not in roles:
                continue
            if atom_name not in name_to_idx:
                continue

            h_names = _attached_h_names(mol, name_to_idx, atom_name)
            if not h_names:
                raise ValueError(f"Uncovered donor has no attached H atoms (cannot build bb_cnh frame): {atom_name}")
            h_name = h_names[0]  # deterministic
            ca_name = _heavy_neighbor_name(mol, name_to_idx, atom_name)

            # bb_cnh correspond atom order is (CA, N, H)
            lgd_sel = (ca_name, atom_name, h_name)
            for nm in lgd_sel:
                if nm not in name_to_idx:
                    raise KeyError(f"Fallback bb_cnh atom not found in ligand PDB: {nm}")

            R, t = _frame_from_points(_xyz(mol, name_to_idx[ca_name]), _xyz(mol, name_to_idx[atom_name]), _xyz(mol, name_to_idx[h_name]))
            sites.append(
                SiteFrame(
                    site_id=f"fallback:bb_cnh:{atom_name}",
                    cg="bb_cnh",
                    lgd_sel=lgd_sel,
                    correspond_resname="*",
                    correspond_names=("CA", "N", "H"),
                    frame_atom_names=lgd_sel,
                    R=R,
                    t=t,
                    covers_polar_atoms=(atom_name,),
                )
            )
            covered_by_sites.add(atom_name)

    # 3) Optional role-based frames for ligand cations: map to ccn-like frame (CD,CE,NZ)
    if args.add_ccn_for_cations:
        for atom_name in sorted(polar_atom_names):
            roles = polar_sites.get(atom_name, set())
            if "cation" not in roles:
                continue
            if atom_name not in name_to_idx:
                continue
            try:
                a0_name, a1_name = _pick_two_heavy_neighbors_for_ccn_frame(mol, name_to_idx, atom_name)
            except ValueError:
                # Skip if we can't construct a stable frame for this atom.
                continue

            lgd_sel = (a0_name, a1_name, atom_name)  # (CD, CE, NZ)
            for nm in lgd_sel:
                if nm not in name_to_idx:
                    raise KeyError(f"Fallback ccn atom not found in ligand PDB: {nm}")

            R, t = _frame_from_points(_xyz(mol, name_to_idx[a0_name]), _xyz(mol, name_to_idx[a1_name]), _xyz(mol, name_to_idx[atom_name]))
            sites.append(
                SiteFrame(
                    site_id=f"fallback:ccn:{atom_name}",
                    cg="ccn",
                    lgd_sel=lgd_sel,
                    correspond_resname="LYS",
                    correspond_names=("CD", "CE", "NZ"),
                    frame_atom_names=lgd_sel,
                    R=R,
                    t=t,
                    covers_polar_atoms=(atom_name,),
                )
            )
            covered_by_sites.add(atom_name)

    payload = {
        "inputs": {
            "ligand_pdb": str(args.ligand_pdb),
            "cg_atommap": str(args.cg_atommap),
            "polar_sites": str(args.polar_sites),
        },
        "n_sites": len(sites),
        "covered_polar_atoms": sorted(covered_by_sites),
        "uncovered_polar_atoms": sorted(polar_atom_names - covered_by_sites),
        "sites": [s.to_json() for s in sites],
    }

    _write_json(args.out, payload)


if __name__ == "__main__":
    main()
