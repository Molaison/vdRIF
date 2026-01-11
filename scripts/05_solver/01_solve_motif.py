#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from ortools.sat.python import cp_model


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _stable_u64(s: str) -> np.uint64:
    h = hashlib.blake2b(s.encode("utf-8"), digest_size=8).digest()
    return np.frombuffer(h, dtype="<u8")[0]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _unpack_xform12(x12: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x12 = np.asarray(x12)
    n = x12.shape[0]
    R = np.zeros((n, 3, 3), dtype=np.float64)
    t = np.zeros((n, 3), dtype=np.float64)
    R[:, 0, 0] = x12[:, 0]
    R[:, 0, 1] = x12[:, 1]
    R[:, 0, 2] = x12[:, 2]
    t[:, 0] = x12[:, 3]
    R[:, 1, 0] = x12[:, 4]
    R[:, 1, 1] = x12[:, 5]
    R[:, 1, 2] = x12[:, 6]
    t[:, 1] = x12[:, 7]
    R[:, 2, 0] = x12[:, 8]
    R[:, 2, 1] = x12[:, 9]
    R[:, 2, 2] = x12[:, 10]
    t[:, 2] = x12[:, 11]
    return R, t


def _xform_apply(R: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
    return (R @ p) + t


def _candidate_stub_atoms(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Return (4,3) coordinates of [N,CA,C,CB] in world.
    These are the same local coordinates as used in rifdock BackboneActor.hh.
    """
    local = np.array(
        [
            [2.80144, -0.992889, -1.52486],  # N
            [1.95280, 0.220007, -1.52486],  # CA
            [2.87767, 1.43290, -1.52486],  # C
            [1.0264273, 0.25245885, -0.308907],  # CB
        ],
        dtype=np.float64,
    )
    return np.stack([_xform_apply(R, t, p) for p in local], axis=0)


def _pair_clash(a_xyz: np.ndarray, b_xyz: np.ndarray, tol: float) -> bool:
    """
    Very fast approximate residue-residue clash using N/CA/C/CB only.
    """
    # vdw radii (Angstrom) for [N,CA,C,CB] ~ [N,C,C,C]
    a_vdw = np.array([1.55, 1.70, 1.70, 1.70], dtype=np.float64)
    b_vdw = a_vdw
    thr = (a_vdw[:, None] + b_vdw[None, :] - tol) ** 2  # (4,4)
    d = a_xyz[:, None, :] - b_xyz[None, :, :]
    d2 = np.einsum("ijx,ijx->ij", d, d)
    return bool((d2 < thr).any())


def _atom_order_vdw(atom_order: list[str]) -> np.ndarray:
    def elem_from_atom_name(nm: str) -> str:
        # ATOM_ORDER uses standard PDB atom names; for amino acids we only expect C/N/O/S.
        return "S" if nm.startswith("S") else nm[0]

    vdw_by_elem = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
    elems = [elem_from_atom_name(nm) for nm in atom_order]
    return np.asarray([vdw_by_elem.get(e, 1.70) for e in elems], dtype=np.float64)


def _pair_clash_full_atom(
    a_xyz: np.ndarray,  # (m,3) with NaNs for missing atoms
    b_xyz: np.ndarray,  # (m,3) with NaNs for missing atoms
    vdw: np.ndarray,  # (m,)
    tol: float,
) -> bool:
    a_ok = np.isfinite(a_xyz).all(axis=1)
    b_ok = np.isfinite(b_xyz).all(axis=1)
    if not np.any(a_ok) or not np.any(b_ok):
        return False
    aa = a_xyz[a_ok]
    bb = b_xyz[b_ok]
    av = vdw[a_ok]
    bv = vdw[b_ok]
    d = aa[:, None, :] - bb[None, :, :]
    d2 = np.einsum("ijx,ijx->ij", d, d)
    thr = (av[:, None] + bv[None, :] - float(tol)) ** 2
    return bool((d2 < thr).any())


def _build_conflicts_grid(
    R: np.ndarray,
    t: np.ndarray,
    ca_only: np.ndarray,
    center_atom_xyz_stub: np.ndarray | None,
    atom_order: list[str] | None,
    grid_size: float,
    ca_prefilter: float,
    tol: float,
) -> list[tuple[int, int]]:
    """
    Build (i,j) conflicts with a uniform grid over CA positions for near-neighbor pruning.
    """
    assert ca_only.shape[0] == R.shape[0] == t.shape[0]
    n = ca_only.shape[0]
    inv = 1.0 / float(grid_size)
    use_full = center_atom_xyz_stub is not None and atom_order is not None
    vdw = _atom_order_vdw(atom_order) if use_full else None
    world_all = None
    if use_full:
        loc = center_atom_xyz_stub.astype(np.float64)
        world_all = np.einsum("nij,nkj->nki", R, loc) + t[:, None, :]

    def cell_key(p: np.ndarray) -> tuple[int, int, int]:
        return (int(math.floor(p[0] * inv)), int(math.floor(p[1] * inv)), int(math.floor(p[2] * inv)))

    cells: dict[tuple[int, int, int], list[int]] = {}
    for i in range(n):
        k = cell_key(ca_only[i])
        cells.setdefault(k, []).append(i)

    # Need to search enough neighboring cells so that any pair within ca_prefilter is considered.
    # Range in cells is ceil(ca_prefilter / grid_size) per dimension.
    r = int(math.ceil(float(ca_prefilter) / float(grid_size)))
    neighbors = [(dx, dy, dz) for dx in range(-r, r + 1) for dy in range(-r, r + 1) for dz in range(-r, r + 1)]

    conflicts: list[tuple[int, int]] = []
    ca2 = float(ca_prefilter) ** 2

    for key, idxs in cells.items():
        # Compare within this cell and neighbor cells, but avoid duplicates by enforcing i<j and ordering by cell key.
        for dx, dy, dz in neighbors:
            nk = (key[0] + dx, key[1] + dy, key[2] + dz)
            if nk not in cells:
                continue
            jdxs = cells[nk]
            for i in idxs:
                for j in jdxs:
                    if j <= i:
                        continue
                    # CA prefilter
                    d = ca_only[i] - ca_only[j]
                    if float(d @ d) > ca2:
                        continue
                    if use_full:
                        assert world_all is not None and vdw is not None
                        if _pair_clash_full_atom(world_all[i], world_all[j], vdw=vdw, tol=tol):
                            conflicts.append((i, j))
                    else:
                        a_atoms = _candidate_stub_atoms(R[i], t[i])
                        b_atoms = _candidate_stub_atoms(R[j], t[j])
                        if _pair_clash(a_atoms, b_atoms, tol=tol):
                            conflicts.append((i, j))

    return conflicts


def _bit_covered(mask: np.ndarray, bit: int) -> np.ndarray:
    return (mask & (1 << bit)) != 0


@dataclass(frozen=True)
class SolveParams:
    min_res: int
    max_res: int
    time_limit_s: float
    num_workers: int
    grid_size: float
    ca_prefilter: float
    clash_tol: float


def _dump_motif_pdb(
    ligand_pdb: Path,
    out_pdb: Path,
    selected: list[int],
    aa3: np.ndarray,
    xform12: np.ndarray,
    center_atom_xyz_stub: np.ndarray | None,
    atom_order: list[str] | None,
) -> None:
    """
    PDB: ligand coordinates + selected vdM interaction residues (chain X resnum==10),
    reconstructed from `center_atom_xyz_stub_f32` in the rifdock BackboneActor stub frame.
    """
    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    # Keep ligand as-is, then append motif residues
    ligand_lines = ligand_pdb.read_text(encoding="utf-8").splitlines()
    # strip trailing END if present to append cleanly
    while ligand_lines and ligand_lines[-1].strip() == "END":
        ligand_lines.pop()

    R, t = _unpack_xform12(xform12)

    atom_serial = 1
    # renumber after ligand atoms if possible
    for line in ligand_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                atom_serial = max(atom_serial, int(line[6:11].strip()) + 1)
            except Exception:
                pass

    lines = list(ligand_lines)
    chain = "M"
    if center_atom_xyz_stub is None or atom_order is None:
        raise ValueError(
            "Missing full-atom residue coordinates (center_atom_xyz_stub_f32) or atom_order; "
            "regenerate `processed/03_vdxform/<cg>/atom_order.json` and candidates NPZ."
        )
    for i_out, i in enumerate(selected, start=1):
        resn = str(aa3[i])
        Ri = R[i]
        ti = t[i]
        # Full-atom center residue coordinates (heavy atoms only), skipping NaNs.
        for an, p in zip(atom_order, center_atom_xyz_stub[i]):
            if not np.isfinite(p).all():
                continue
            xyz = (Ri @ p.astype(np.float64)) + ti
            elem = an[0]
            if elem == "C":
                ae = "C"
            elif elem == "N":
                ae = "N"
            elif elem == "O":
                ae = "O"
            elif elem == "S":
                ae = "S"
            else:
                ae = "C"
            lines.append(
                f"ATOM  {atom_serial:5d} {an:>4s} {resn:>3s} {chain}{i_out:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {ae:>2s}"
            )
            atom_serial += 1
    lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Deterministic 8–15 residue motif solver (CP-SAT).")
    ap.add_argument("--candidates-npz", type=Path, required=True)
    ap.add_argument("--candidates-meta", type=Path, required=True)
    ap.add_argument("--ligand-pdb", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-pdb", type=Path, required=True)
    ap.add_argument("--min-res", type=int, default=8)
    ap.add_argument("--max-res", type=int, default=15)
    ap.add_argument("--time-limit-s", type=float, default=60.0)
    ap.add_argument("--num-workers", type=int, default=1, help="Set to 1 for determinism.")
    ap.add_argument("--grid-size", type=float, default=4.0)
    ap.add_argument("--ca-prefilter", type=float, default=12.0)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    args = ap.parse_args()

    params = SolveParams(
        min_res=int(args.min_res),
        max_res=int(args.max_res),
        time_limit_s=float(args.time_limit_s),
        num_workers=int(args.num_workers),
        grid_size=float(args.grid_size),
        ca_prefilter=float(args.ca_prefilter),
        clash_tol=float(args.clash_tol),
    )
    if params.num_workers != 1:
        raise ValueError("For determinism, require --num-workers=1.")

    meta = _load_json(args.candidates_meta)
    polar_atoms: list[str] = meta["polar_atoms"]
    atom_order = meta.get("atom_order")
    if len(polar_atoms) > 16:
        raise ValueError("MVP expects <=16 polar atoms.")

    z = np.load(args.candidates_npz, allow_pickle=True)
    cand_id = z["cand_id_u64"].astype(np.uint64)
    cover = z["cover_mask_u16"].astype(np.uint16)
    score_f = z["score_f32"].astype(np.float64)
    aa3 = z["aa3"]
    xform12 = z["xform_world_stub_12_f32"]
    center_atom_xyz_stub = z["center_atom_xyz_stub_f32"] if "center_atom_xyz_stub_f32" in z.files else None

    n = int(cand_id.shape[0])
    if n == 0:
        raise ValueError("No candidates provided.")

    # Deterministic ranking for tie-break.
    order = sorted(range(n), key=lambda i: (-float(score_f[i]), int(cand_id[i])))
    rank = np.empty((n,), dtype=np.int32)
    for r, i in enumerate(order):
        rank[i] = r

    # Build conflicts (motif internal clashes) using CA grid for pruning.
    R, t = _unpack_xform12(xform12)
    # IMPORTANT: Do not assume a specific stub-frame CA local coordinate unless we are sure our
    # stub convention matches rifdock exactly. Prefer the CA atom coordinate from the same
    # `center_atom_xyz_stub_f32` used to write the motif PDB.
    if center_atom_xyz_stub is not None and atom_order is not None and "CA" in atom_order:
        ca_i = int(atom_order.index("CA"))
        ca_loc = center_atom_xyz_stub[:, ca_i, :].astype(np.float64)
        ca_xyz = np.einsum("nij,nj->ni", R, ca_loc) + t
    else:
        ca_local = np.array([1.95280, 0.220007, -1.52486], dtype=np.float64)
        ca_xyz = np.einsum("nij,j->ni", R, ca_local) + t
    conflicts = _build_conflicts_grid(
        R=R,
        t=t,
        ca_only=ca_xyz,
        center_atom_xyz_stub=center_atom_xyz_stub,
        atom_order=atom_order,
        grid_size=params.grid_size,
        ca_prefilter=params.ca_prefilter,
        tol=params.clash_tol,
    )

    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x_{i}") for i in range(n)]

    # Cardinality
    model.Add(sum(x) >= params.min_res)
    model.Add(sum(x) <= params.max_res)

    # Coverage constraints
    for b, atom_name in enumerate(polar_atoms):
        idxs = [i for i in range(n) if int(cover[i]) & (1 << b)]
        if not idxs:
            raise ValueError(f"No candidates cover polar atom {atom_name} (bit {b}).")
        model.Add(sum(x[i] for i in idxs) >= 1)

    # Clash constraints
    for i, j in conflicts:
        model.Add(x[i] + x[j] <= 1)

    # Objective: lexicographic (encoded as 1 linear objective):
    # 1) minimize number of residues (prefer smaller motifs; satisfies user's \"8-15\" as a bound)
    # 2) maximize total score (within the minimal count)
    # 3) deterministic tie-break by candidate rank
    #
    # Encode as: maximize sum( gain_i - B ) * x_i where B is large enough that any +1 residue loses
    # against any possible gain difference, thus enforcing \"minimize count\" first.
    score_int = np.round(np.clip(score_f, -1e6, 1e6) * 1000.0).astype(np.int64)
    tie_int = (n - 1 - rank).astype(np.int64)  # earlier rank => larger tie
    gain = score_int + tie_int
    span = int(params.max_res) * int(gain.max() - gain.min())
    B = span + 1
    coeff = gain - B
    model.Maximize(sum(int(coeff[i]) * x[i] for i in range(n)))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = params.time_limit_s
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.log_search_progress = False

    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(f"CP-SAT failed with status {solver.StatusName(status)}")

    selected = [i for i in range(n) if solver.Value(x[i]) == 1]
    selected_sorted = sorted(selected, key=lambda i: (-float(score_f[i]), int(cand_id[i])))

    # Coverage report
    cov = 0
    for i in selected_sorted:
        cov |= int(cover[i])
    need = (1 << len(polar_atoms)) - 1

    report = {
        "inputs": {
            "ligand_pdb": str(args.ligand_pdb),
            "candidates_npz_sha256": _sha256(args.candidates_npz),
            "candidates_meta_sha256": _sha256(args.candidates_meta),
        },
        "params": {
            "min_res": params.min_res,
            "max_res": params.max_res,
            "time_limit_s": params.time_limit_s,
            "num_workers": params.num_workers,
            "grid_size": params.grid_size,
            "ca_prefilter": params.ca_prefilter,
            "clash_tol": params.clash_tol,
        },
        "solver": {
            "status": solver.StatusName(status),
            "objective_value": float(solver.ObjectiveValue()),
            "n_conflicts": len(conflicts),
        },
        "n_candidates": n,
        "n_selected": len(selected_sorted),
        "coverage_bitmask": int(cov),
        "coverage_complete": bool(cov == need),
        "polar_atoms": polar_atoms,
        "selected": [
            {
                "cand_id_u64": int(cand_id[i]),
                "score": float(score_f[i]),
                "aa3": str(aa3[i]),
                "cover_mask_u16": int(cover[i]),
            }
            for i in selected_sorted
        ],
    }

    _write_json(args.out_json, report)
    _dump_motif_pdb(
        args.ligand_pdb,
        args.out_pdb,
        selected_sorted,
        aa3=aa3,
        xform12=xform12,
        center_atom_xyz_stub=center_atom_xyz_stub,
        atom_order=atom_order,
    )


if __name__ == "__main__":
    main()
