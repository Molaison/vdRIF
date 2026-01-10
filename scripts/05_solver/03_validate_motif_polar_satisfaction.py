#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_pdb_atoms(pdb_path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21:22].strip()
        resi = int(line[22:26].strip())
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        elem = line[76:78].strip() if len(line) >= 78 else name[0]
        atoms.append({"name": name, "resn": resn, "chain": chain, "resi": resi, "elem": elem, "xyz": np.array([x, y, z])})
    return atoms


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate that each ligand polar atom has a nearby complementary motif atom.")
    ap.add_argument("--polar-sites", type=Path, required=True)
    ap.add_argument("--motif-pdb", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--hbond-dist", type=float, default=3.5, help="Heavy-atom distance cutoff for Hbond-like satisfaction.")
    ap.add_argument("--ion-dist", type=float, default=4.0, help="Distance cutoff for ionic satisfaction.")
    args = ap.parse_args()

    polar = _load_json(args.polar_sites)
    atoms = _parse_pdb_atoms(args.motif_pdb)

    # Separate ligand atoms (chain/resseq from MTX input; we identify by residue name MTX)
    ligand_atoms = {a["name"]: a for a in atoms if a["resn"] == "MTX"}
    motif_atoms = [a for a in atoms if a["resn"] != "MTX"]

    # Very simple atom typing for motif residues (heavy atoms only).
    donor_atoms = {"N", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "ND1", "ND2", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD", "SG"}
    cation_atoms = {"NZ", "NH1", "NH2"}  # minimal
    anion_atoms = {"OD1", "OD2", "OE1", "OE2"}

    hbd2 = float(args.hbond_dist) ** 2
    ion2 = float(args.ion_dist) ** 2

    results = []
    ok_all = True
    for s in polar["sites"]:
        name = s["atom_name"]
        roles = set(s["roles"])
        if name not in ligand_atoms:
            results.append({"atom_name": name, "roles": sorted(roles), "found_in_pdb": False, "satisfied": False})
            ok_all = False
            continue
        p = ligand_atoms[name]["xyz"]

        # Determine what we need in the motif
        needs: list[tuple[str, set[str], float]] = []
        if "acceptor" in roles:
            needs.append(("donor", donor_atoms, hbd2))
        if "donor" in roles:
            needs.append(("acceptor", acceptor_atoms, hbd2))
        if "cation" in roles:
            needs.append(("anion", anion_atoms, ion2))
        if "anion" in roles:
            needs.append(("cation", cation_atoms, ion2))

        sat = True
        sat_details = []
        for need_name, allowed, cutoff2 in needs:
            best = None
            best_d2 = None
            for a in motif_atoms:
                if a["name"] not in allowed:
                    continue
                d = a["xyz"] - p
                d2 = float(d @ d)
                if d2 <= cutoff2 and (best_d2 is None or d2 < best_d2):
                    best_d2 = d2
                    best = a
            if best is None:
                sat = False
                sat_details.append({"need": need_name, "found": False})
            else:
                sat_details.append(
                    {
                        "need": need_name,
                        "found": True,
                        "motif_atom": f"{best['resn']} {best['chain']}{best['resi']}:{best['name']}",
                        "dist": float(np.sqrt(best_d2)),
                    }
                )

        ok_all = ok_all and sat
        results.append({"atom_name": name, "roles": sorted(roles), "found_in_pdb": True, "satisfied": sat, "details": sat_details})

    _write_json(
        args.out,
        {
            "inputs": {"polar_sites": str(args.polar_sites), "motif_pdb": str(args.motif_pdb)},
            "params": {"hbond_dist": float(args.hbond_dist), "ion_dist": float(args.ion_dist)},
            "all_satisfied": bool(ok_all),
            "results": results,
        },
    )

    if not ok_all:
        raise SystemExit(f"Polar satisfaction FAILED; see {args.out}")


if __name__ == "__main__":
    main()

