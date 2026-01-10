#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _load_pdb_atoms(pdb_path: Path) -> list[dict[str, Any]]:
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
        elem = line[76:78].strip() if len(line) >= 78 else ""
        if not elem:
            elem = name[0]
        atoms.append({"name": name, "resn": resn, "chain": chain, "resi": resi, "elem": elem, "xyz": np.array([x, y, z])})
    return atoms


def _vdw_radius(elem: str) -> float:
    e = elem.strip().upper()
    # common 2-letter elements in PDB element column
    if len(e) == 2 and e[1].islower():
        e = e[0] + e[1].upper()
    vdw = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
    return float(vdw.get(e, 1.70))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate that the motif has no severe ligand–residue clashes (vdW overlap).")
    ap.add_argument("--motif-pdb", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--ligand-resname", type=str, default="MTX")
    ap.add_argument("--tol", type=float, default=0.5, help="Tolerance subtracted from vdw sum; matches candidate generator.")
    ap.add_argument("--max-report", type=int, default=50)
    ap.add_argument("--fail-overlap", type=float, default=0.5, help="Fail if any pair overlap exceeds this (Å).")
    args = ap.parse_args()

    atoms = _load_pdb_atoms(args.motif_pdb)
    lig = [a for a in atoms if a["resn"] == args.ligand_resname and a["elem"].upper() != "H"]
    mot = [a for a in atoms if a["resn"] != args.ligand_resname and a["elem"].upper() != "H"]

    pairs = []
    worst = 0.0
    for la in lig:
        rl = _vdw_radius(la["elem"])
        for ma in mot:
            rm = _vdw_radius(ma["elem"])
            d = ma["xyz"] - la["xyz"]
            dist = float(np.sqrt(d @ d))
            thr = rl + rm - float(args.tol)
            overlap = float(thr - dist)
            if overlap > 0.0:
                worst = max(worst, overlap)
                pairs.append(
                    {
                        "overlap": overlap,
                        "dist": dist,
                        "thr": thr,
                        "ligand_atom": la["name"],
                        "motif_atom": f"{ma['resn']} {ma['chain']}{ma['resi']}:{ma['name']}",
                    }
                )

    pairs_sorted = sorted(pairs, key=lambda x: (-float(x["overlap"]), float(x["dist"])))[: int(args.max_report)]
    ok = bool(worst <= float(args.fail_overlap))
    _write_json(
        args.out,
        {
            "inputs": {"motif_pdb": str(args.motif_pdb), "ligand_resname": str(args.ligand_resname)},
            "params": {"tol": float(args.tol), "fail_overlap": float(args.fail_overlap)},
            "ok": ok,
            "worst_overlap": float(worst),
            "n_clashing_pairs": int(len(pairs)),
            "top_clashes": pairs_sorted,
        },
    )
    if not ok:
        raise SystemExit(f"Clash validation FAILED (worst_overlap={worst:.3f}); see {args.out}")


if __name__ == "__main__":
    main()

