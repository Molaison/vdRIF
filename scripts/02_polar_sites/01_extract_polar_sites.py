#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdDetermineBonds


@dataclass(frozen=True)
class AtomSite:
    atom_name: str
    atom_idx: int
    element: str
    formal_charge: int
    x: float
    y: float
    z: float
    roles: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "atom_idx": self.atom_idx,
            "element": self.element,
            "formal_charge": self.formal_charge,
            "xyz": [self.x, self.y, self.z],
            "roles": list(self.roles),
        }


def _feature_factory() -> ChemicalFeatures.FreeChemicalFeatureFactory:
    fdef = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
    return ChemicalFeatures.BuildFeatureFactory(str(fdef))


def _load_pdb_ligand(pdb_path: Path) -> Chem.Mol:
    # PDB has no bond orders; relying on default RDKit perception can create spurious implicit Hs,
    # which breaks donor/acceptor typing. Instead:
    # 1) read without sanitization
    # 2) infer bonds/bond orders from 3D geometry (rdDetermineBonds)
    # 3) sanitize
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(f"Failed to read ligand PDB: {pdb_path}")
    rdDetermineBonds.DetermineBondOrders(mol)
    Chem.SanitizeMol(mol)
    if mol.GetNumConformers() != 1:
        raise ValueError(f"Expected exactly 1 conformer in {pdb_path}, got {mol.GetNumConformers()}")
    return mol


def _atom_name(atom: Chem.Atom) -> str:
    info = atom.GetPDBResidueInfo()
    if info is None:
        return f"ATOM{atom.GetIdx()}"
    return info.GetName().strip()


def _coords(mol: Chem.Mol, atom_idx: int) -> np.ndarray:
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(atom_idx)
    return np.array([p.x, p.y, p.z], dtype=np.float64)


def _collect_feature_roles(mol: Chem.Mol) -> dict[int, set[str]]:
    ff = _feature_factory()
    roles_by_atom: dict[int, set[str]] = {}

    for feat in ff.GetFeaturesForMol(mol):
        fam = feat.GetFamily()
        if fam == "Donor":
            role = "donor"
        elif fam == "Acceptor":
            role = "acceptor"
        else:
            continue

        for aidx in feat.GetAtomIds():
            roles_by_atom.setdefault(aidx, set()).add(role)

    # Add formal-charge-derived ion roles (deterministic)
    for atom in mol.GetAtoms():
        chg = int(atom.GetFormalCharge())
        if chg > 0:
            roles_by_atom.setdefault(atom.GetIdx(), set()).add("cation")
        elif chg < 0:
            roles_by_atom.setdefault(atom.GetIdx(), set()).add("anion")

    return roles_by_atom


def _iter_sites(mol: Chem.Mol, include_h: bool) -> Iterable[AtomSite]:
    roles_by_atom = _collect_feature_roles(mol)
    for atom in mol.GetAtoms():
        if not include_h and atom.GetSymbol() == "H":
            continue
        roles_set = set(roles_by_atom.get(atom.GetIdx(), set()))
        # Guardrail: treat "donor" as requiring at least one attached H (explicit or implicit).
        if "donor" in roles_set and int(atom.GetTotalNumHs(includeNeighbors=True)) <= 0:
            roles_set.remove("donor")
        roles = tuple(sorted(roles_set))
        if not roles:
            continue
        xyz = _coords(mol, atom.GetIdx())
        yield AtomSite(
            atom_name=_atom_name(atom),
            atom_idx=int(atom.GetIdx()),
            element=atom.GetSymbol(),
            formal_charge=int(atom.GetFormalCharge()),
            x=float(xyz[0]),
            y=float(xyz[1]),
            z=float(xyz[2]),
            roles=roles,
        )


def main() -> None:
    ap = argparse.ArgumentParser(description="Extract ligand polar sites (HBD/HBA/ion) using RDKit features.")
    ap.add_argument("ligand_pdb", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--include-h", action="store_true", help="Include hydrogen atoms if they carry roles (rare).")
    args = ap.parse_args()

    mol = _load_pdb_ligand(args.ligand_pdb)
    sites = sorted(_iter_sites(mol, include_h=bool(args.include_h)), key=lambda s: (s.atom_name, s.atom_idx))

    payload = {
        "input": {"ligand_pdb": str(args.ligand_pdb)},
        "n_sites": len(sites),
        "sites": [s.to_json() for s in sites],
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
