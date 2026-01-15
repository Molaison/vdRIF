#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import heapq
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem


def _stable_u64(s: str) -> np.uint64:
    h = hashlib.blake2b(s.encode("utf-8"), digest_size=8).digest()
    return np.frombuffer(h, dtype="<u8")[0]


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


def _ligand_heavy_atoms(mol: Chem.Mol) -> tuple[list[str], np.ndarray, np.ndarray]:
    conf = mol.GetConformer()
    names: list[str] = []
    xyz: list[list[float]] = []
    vdw: list[float] = []

    # Reasonable defaults for quick clash filtering (Angstrom)
    vdw_by_elem = {
        "H": 1.20,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Br": 1.85,
        "I": 1.98,
    }

    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        p = conf.GetAtomPosition(a.GetIdx())
        names.append(_atom_name(a))
        xyz.append([p.x, p.y, p.z])
        vdw.append(vdw_by_elem.get(a.GetSymbol(), 1.70))

    return names, np.asarray(xyz, dtype=np.float64), np.asarray(vdw, dtype=np.float64)


def _unpack_xform12(x12: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    x12: shape (n,12) float32/float64, packed row-major 3x4 (R|t)
    Returns (R,t) with shapes (n,3,3) and (n,3)
    """
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


def _pack_xform12(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    n = R.shape[0]
    x12 = np.zeros((n, 12), dtype=np.float32)
    x12[:, 0] = R[:, 0, 0]
    x12[:, 1] = R[:, 0, 1]
    x12[:, 2] = R[:, 0, 2]
    x12[:, 3] = t[:, 0]
    x12[:, 4] = R[:, 1, 0]
    x12[:, 5] = R[:, 1, 1]
    x12[:, 6] = R[:, 1, 2]
    x12[:, 7] = t[:, 1]
    x12[:, 8] = R[:, 2, 0]
    x12[:, 9] = R[:, 2, 1]
    x12[:, 10] = R[:, 2, 2]
    x12[:, 11] = t[:, 2]
    return x12


def _xform_apply(R: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
    # R: (n,3,3) t: (n,3) p: (3,)
    return np.einsum("nij,j->ni", R, p) + t


def _compose_xforms(
    R_a: np.ndarray,
    t_a: np.ndarray,
    R_b: np.ndarray,
    t_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compose X = A * B.
    A: (R_a, t_a), B: (R_b, t_b), each with batch dimension (n,...).
    Returns (R, t).
    """
    R = np.einsum("nij,njk->nik", R_a, R_b)
    t = np.einsum("nij,nj->ni", R_a, t_b) + t_a
    return R, t


def _is_finite(x: np.ndarray) -> np.ndarray:
    return np.isfinite(x)


def _prior_score(cluster_rank: np.ndarray, c_score: np.ndarray) -> np.ndarray:
    """
    Deterministic, simple prior:
    - prefer lower cluster_rank (rank 1 best) -> map to positive score
    - else fall back to c_score if present
    """
    score = np.zeros_like(c_score, dtype=np.float64)
    have_rank = _is_finite(cluster_rank) & (cluster_rank > 0)
    # Use -log(rank) so rank=1 gives 0, rank=10 gives ~-2.3; then negate to make higher better.
    score[have_rank] = -np.log(cluster_rank[have_rank])

    need = ~have_rank
    have_c = need & _is_finite(c_score)
    score[have_c] = c_score[have_c]
    return score


def _clash_mask(
    residue_atoms_xyz: np.ndarray,  # (n,k,3)
    residue_atoms_vdw: np.ndarray,  # (k,)
    ligand_xyz: np.ndarray,  # (m,3)
    ligand_vdw: np.ndarray,  # (m,)
    tolerance: float,
) -> np.ndarray:
    """
    Returns boolean mask (n,) where True means "no clash".
    """
    # dist^2: (n,k,m)
    d = residue_atoms_xyz[:, :, None, :] - ligand_xyz[None, None, :, :]
    d2 = np.einsum("nkmj,nkmj->nkm", d, d)
    # threshold: (k,m)
    thr = residue_atoms_vdw[:, None] + ligand_vdw[None, :] - float(tolerance)
    thr2 = thr * thr
    clash = (d2 < thr2[None, :, :]).any(axis=(1, 2))
    return ~clash


def _clash_mask_full_fast(
    residue_atoms_xyz: np.ndarray,  # (n,k,3)
    residue_atoms_vdw: np.ndarray,  # (k,)
    ligand_xyz: np.ndarray,  # (m,3)
    ligand_vdw: np.ndarray,  # (m,)
    tolerance: float,
) -> np.ndarray:
    """
    Returns boolean mask (n,) where True means "no clash".

    Full-atom clash filter is a major hotspot; this implementation avoids allocating
    a 4D (n,k,m,3) tensor by using the identity:
      ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a·b
    which can be computed via a single matmul on (n*k,3)×(3,m).

    NaNs (missing atoms) propagate through and do not count as clashes.
    """
    residue_atoms_xyz = np.asarray(residue_atoms_xyz, dtype=np.float32)
    ligand_xyz = np.asarray(ligand_xyz, dtype=np.float32)
    residue_atoms_vdw = np.asarray(residue_atoms_vdw, dtype=np.float32)
    ligand_vdw = np.asarray(ligand_vdw, dtype=np.float32)

    n = int(residue_atoms_xyz.shape[0])
    if n == 0:
        return np.zeros((0,), dtype=bool)
    k = int(residue_atoms_xyz.shape[1])
    m = int(ligand_xyz.shape[0])

    res = residue_atoms_xyz.reshape(n * k, 3)
    lig = ligand_xyz

    res2 = np.sum(res * res, axis=1, dtype=np.float32)[:, None]  # (n*k,1)
    lig2 = np.sum(lig * lig, axis=1, dtype=np.float32)[None, :]  # (1,m)
    dot = res @ lig.T  # (n*k,m); NaNs propagate
    d2 = (res2 + lig2) - (2.0 * dot)
    d2 = d2.reshape(n, k, m)

    thr = residue_atoms_vdw[:, None] + ligand_vdw[None, :] - float(tolerance)  # (k,m)
    thr2 = thr * thr
    clash = (d2 < thr2[None, :, :]).any(axis=(1, 2))
    return ~clash


@dataclass(frozen=True)
class Site:
    site_index: int
    site_id: str
    cg: str
    R_ifg: np.ndarray  # 3x3
    t_ifg: np.ndarray  # 3
    cover_atom_names: tuple[str, ...]


def _needed_roles(roles: set[str]) -> list[str]:
    need: list[str] = []
    if "acceptor" in roles:
        need.append("donor")
    if "donor" in roles:
        need.append("acceptor")
    if "cation" in roles:
        need.append("anion")
    if "anion" in roles:
        need.append("cation")
    return need


def _to_sites(site_frames_json: dict[str, Any], polar_atom_to_bit: dict[str, int]) -> list[Site]:
    sites: list[Site] = []
    for i, s in enumerate(site_frames_json["sites"]):
        R = np.asarray(s["R"], dtype=np.float64)
        t = np.asarray(s["t"], dtype=np.float64)
        cover_atoms = tuple(str(a) for a in (s.get("covers_polar_atoms", []) or []) if a in polar_atom_to_bit)
        sites.append(
            Site(
                site_index=i,
                site_id=str(s["site_id"]),
                cg=str(s["cg"]),
                R_ifg=R,
                t_ifg=t,
                cover_atom_names=cover_atoms,
            )
        )
    return sites


def _load_vdxform_npz(path: Path) -> dict[str, np.ndarray]:
    z = np.load(path, allow_pickle=True)
    return {k: z[k] for k in z.files}


def _as_aa3_str(aa_arr: np.ndarray) -> np.ndarray:
    # Accept object array, bytes array, or unicode
    if aa_arr.dtype.kind == "S":
        return aa_arr.astype("U3")
    if aa_arr.dtype.kind == "U":
        return aa_arr.astype("U3")
    # object
    return aa_arr.astype("U3")


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate vdM-based residue placement candidates around a fixed ligand.")
    ap.add_argument("--ligand-pdb", type=Path, required=True)
    ap.add_argument("--polar-sites", type=Path, required=True)
    ap.add_argument("--site-frames", type=Path, required=True)
    ap.add_argument("--vdxform-dir", type=Path, required=True, help="Directory containing per-cg vdxform_<cg>.npz files.")
    ap.add_argument("--out-prefix", type=Path, required=True, help="Prefix for output files (npz + json).")
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--top-per-site", type=int, default=2000)
    ap.add_argument("--clash-tol", type=float, default=0.5)
    ap.add_argument(
        "--exclude-aa3",
        type=str,
        default="PRO,CYS",
        help="Comma-separated 3-letter amino acids to exclude (default: PRO,CYS).",
    )
    ap.add_argument(
        "--sidechain-only-satisfaction",
        action="store_true",
        help="If set, exclude backbone N/O atoms from donor/acceptor satisfaction tests (sidechain-only).",
    )
    ap.add_argument(
        "--require-sidechain-facing",
        action="store_true",
        default=True,
        help=(
            "Require the residue sidechain direction to point toward the satisfied ligand polar site(s). "
            "This filters out candidates where only the backbone interacts and the sidechain points away."
        ),
    )
    ap.add_argument(
        "--min-sidechain-facing-dot",
        type=float,
        default=0.2,
        help="Minimum dot( sidechain_dir , ligand_dir ) to accept a candidate (default 0.2).",
    )
    ap.add_argument(
        "--min-sidechain-centroid-dot",
        type=float,
        default=0.0,
        help="Minimum dot( sidechain_dir , direction-to-ligand-centroid ) to accept a candidate (default 0.0).",
    )
    ap.add_argument("--require-full-coverage", action="store_true")
    args = ap.parse_args()

    polar = _load_json(args.polar_sites)
    polar_atoms = [s["atom_name"] for s in polar["sites"]]
    if len(polar_atoms) > 16:
        raise ValueError(f"MVP uses uint16 coverage bitmask; got {len(polar_atoms)} polar atoms (>16).")
    polar_atom_to_bit = {name: i for i, name in enumerate(polar_atoms)}

    frames = _load_json(args.site_frames)
    sites = _to_sites(frames, polar_atom_to_bit)

    polar_by_name = {s["atom_name"]: set(s["roles"]) for s in polar["sites"]}
    polar_xyz_by_name = {s["atom_name"]: np.asarray(s["xyz"], dtype=np.float64) for s in polar["sites"]}

    ligand = _load_pdb_mol(args.ligand_pdb)
    lig_names, lig_xyz, lig_vdw = _ligand_heavy_atoms(ligand)
    lig_centroid = lig_xyz.mean(axis=0)

    # Residue atoms used for clash filtering, in rifdock BackboneActor local coords (Angstrom)
    # See external/rifdock/schemelib/scheme/actor/BackboneActor.hh get_n_ca_c/get_cb.
    res_atom_names = ["N", "CA", "C", "CB"]
    res_atom_elem = ["N", "C", "C", "C"]
    res_atom_local = np.array(
        [
            [2.80144, -0.992889, -1.52486],  # N
            [1.95280, 0.220007, -1.52486],  # CA
            [2.87767, 1.43290, -1.52486],  # C
            [1.0264273, 0.25245885, -0.308907],  # CB
        ],
        dtype=np.float64,
    )
    vdw_by_elem = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20}
    res_atom_vdw = np.array([vdw_by_elem[e] for e in res_atom_elem], dtype=np.float64)

    out_records: list[dict[str, Any]] = []

    # Candidate arrays to write
    cand_id: list[np.uint64] = []
    site_index: list[np.uint16] = []
    vdm_id: list[np.uint64] = []
    aa3: list[str] = []
    score: list[np.float32] = []
    cover_mask: list[np.uint16] = []
    xform_world_stub_12: list[np.ndarray] = []
    center_atom_xyz_stub: list[np.ndarray] = []
    cluster_number: list[np.int32] = []

    # Union of polar atoms that are coverable by the site-frames (not geometric satisfaction).
    coverage_union = 0
    # Union of polar atoms actually satisfied by emitted candidates (geometric, via coarse distance tests).
    satisfied_union = 0
    missing_cg: set[str] = set()

    atom_order: list[str] | None = None

    # Motif atom typing based on atom names (MVP; heavy atoms only)
    donor_atoms = {"N", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "ND1", "ND2", "OG", "OG1", "OH"}
    acceptor_atoms = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD", "SG"}
    cation_atoms = {"NZ", "NH1", "NH2"}
    anion_atoms = {"OD1", "OD2", "OE1", "OE2"}
    backbone_atoms = {"N", "C", "O", "CA"}

    exclude_aa3 = {x.strip().upper() for x in str(args.exclude_aa3).split(",") if x.strip()}

    for site in sites:
        # Only polar atoms explicitly associated with this ligand site-frame are considered coverable from it.
        site_polar = []
        site_cov = 0
        for an in site.cover_atom_names:
            roles = polar_by_name.get(an)
            pxyz = polar_xyz_by_name.get(an)
            if roles is None or pxyz is None:
                continue
            site_polar.append((an, polar_atom_to_bit[an], roles, pxyz))
            site_cov |= 1 << polar_atom_to_bit[an]
        coverage_union |= site_cov

        vdx_path = args.vdxform_dir / site.cg / f"vdxform_{site.cg}.npz"
        if not vdx_path.exists():
            missing_cg.add(site.cg)
            continue

        vdx = _load_vdxform_npz(vdx_path)
        vdm_ids = vdx["vdm_id_u64"].astype(np.uint64)
        cluster_number_i32 = vdx.get("cluster_number_i32")
        if cluster_number_i32 is None:
            raise KeyError(f"vdXform missing cluster_number_i32: {vdx_path}")
        cluster_number_i32 = cluster_number_i32.astype(np.int32)
        aa_arr = _as_aa3_str(vdx["aa"])
        x_ifg_to_stub_12 = vdx["xform_ifg_to_stub_12_f32"]
        R_ifg_to_stub, t_ifg_to_stub = _unpack_xform12(x_ifg_to_stub_12)
        center_xyz_stub = vdx.get("center_atom_xyz_stub_f32")
        if center_xyz_stub is None:
            raise KeyError(f"vdXform missing center_atom_xyz_stub_f32: {vdx_path}")

        # Load atom order once (shared across all cg in this repo).
        if atom_order is None:
            atom_order_path = args.vdxform_dir / site.cg / "atom_order.json"
            if atom_order_path.exists():
                atom_order = _load_json(atom_order_path)["atom_order"]
        if atom_order is None:
            raise RuntimeError("atom_order missing; expected processed/03_vdxform/<cg>/atom_order.json")

        def _elem_from_atom_name(nm: str) -> str:
            # ATOM_ORDER uses standard PDB atom names; element is the first letter for our set (C/N/O/S).
            return "S" if nm.startswith("S") else nm[0]

        atom_order_elems = [_elem_from_atom_name(nm) for nm in atom_order]
        atom_order_vdw = np.array([vdw_by_elem.get(e, 1.70) for e in atom_order_elems], dtype=np.float64)

        def atom_indices(names: set[str]) -> list[int]:
            return [i for i, nm in enumerate(atom_order) if nm in names]

        if args.sidechain_only_satisfaction:
            donor_idx = atom_indices({nm for nm in donor_atoms if nm not in backbone_atoms})
            acceptor_idx = atom_indices({nm for nm in acceptor_atoms if nm not in backbone_atoms})
        else:
            donor_idx = atom_indices(donor_atoms)
            acceptor_idx = atom_indices(acceptor_atoms)
        cation_idx = atom_indices(cation_atoms)
        anion_idx = atom_indices(anion_atoms)

        # Sidechain direction indices (for orientation filtering)
        ca_i = int(atom_order.index("CA")) if "CA" in atom_order else -1
        sidechain_idx = [i for i, nm in enumerate(atom_order) if nm not in backbone_atoms]

        prior = _prior_score(
            vdx.get("cluster_rank_ABPLE_A_f32", np.full((vdm_ids.shape[0],), np.nan, dtype=np.float32)).astype(np.float64),
            vdx.get("C_score_ABPLE_A_f32", np.full((vdm_ids.shape[0],), np.nan, dtype=np.float32)).astype(np.float64),
        )
        prior[~np.isfinite(prior)] = -1e9

        # Keep top candidates for this site (min-heap style arrays)
        top_k = int(args.top_per_site)
        # Deterministic top-k:
        # - primary: higher score is better
        # - tie: smaller vdm_id is better
        # We implement this with a min-heap whose root is the *worst* element:
        # key = (score, -vdm_id)  => for equal score, larger vdm_id yields smaller key (worse).
        heap: list[tuple[float, int, int, int]] = []  # (score, -vdm_id, idx, cov)

        def maybe_add(idx: int, sc: float, cov: int) -> None:
            v = int(vdm_ids[idx])
            item = (float(sc), -v, int(idx), int(cov))
            if len(heap) < top_k:
                heapq.heappush(heap, item)
                return
            if item[:2] > heap[0][:2]:
                heapq.heapreplace(heap, item)

        n = vdm_ids.shape[0]
        chunk = int(args.chunk_size)
        for start in range(0, n, chunk):
            end = min(n, start + chunk)
            Rb = R_ifg_to_stub[start:end]
            tb = t_ifg_to_stub[start:end]

            # Compose X_world_stub = X_world_ifg * X_ifg_to_stub
            Ra = np.broadcast_to(site.R_ifg, (end - start, 3, 3))
            ta = np.broadcast_to(site.t_ifg, (end - start, 3))
            Rw, tw = _compose_xforms(Ra, ta, Rb, tb)

            # Compute residue backbone atom coords for clash filtering
            # xyz = Rw @ p_local + tw
            xyz_atoms = np.stack([_xform_apply(Rw, tw, p) for p in res_atom_local], axis=1)  # (n,k,3)

            ok = _clash_mask(xyz_atoms, res_atom_vdw, lig_xyz, lig_vdw, tolerance=float(args.clash_tol))
            if not np.any(ok):
                continue

            ok_local = np.nonzero(ok)[0]
            if ok_local.size == 0:
                continue

            # We intentionally defer full-atom ligand clash filtering until after
            # *geometric satisfaction* pruning. Only a small fraction of vdMs will
            # satisfy a given polar site, so doing full-atom ligand clash checks
            # on all `ok_local` candidates is a huge performance cost.
            #
            # Step 1: transform only the relevant polar atom sets to evaluate satisfaction.
            Rw_ok = Rw[ok_local]
            tw_ok = tw[ok_local]

            def world_coords_subset(idxs: list[int]) -> np.ndarray:
                if not idxs:
                    return np.zeros((ok_local.size, 0, 3), dtype=np.float64)
                loc = center_xyz_stub[start:end][ok_local][:, idxs, :].astype(np.float64)
                return np.einsum("nij,nkj->nki", Rw_ok, loc) + tw_ok[:, None, :]

            donors_w = world_coords_subset(donor_idx)
            acceptors_w = world_coords_subset(acceptor_idx)
            cations_w = world_coords_subset(cation_idx)
            anions_w = world_coords_subset(anion_idx)

            hbd2 = 3.5 * 3.5
            ion2 = 4.0 * 4.0

            def _min_d2_all(points_xyz: np.ndarray, target_xyz: np.ndarray) -> np.ndarray:
                # points_xyz: (n,k,3); may contain NaNs for missing atoms.
                if points_xyz.shape[1] == 0:
                    return np.full((points_xyz.shape[0],), float("inf"), dtype=np.float64)
                d = points_xyz - target_xyz[None, None, :]
                d2 = np.einsum("nkj,nkj->nk", d, d)  # NaNs propagate
                d2 = np.where(np.isfinite(d2), d2, float("inf"))
                return d2.min(axis=1)

            # Vectorized coverage mask for ok_local candidates.
            cov_arr = np.zeros((ok_local.size,), dtype=np.uint16)
            for _an, bit, roles, pxyz in site_polar:
                need = _needed_roles(roles)
                sat_mask = np.ones((ok_local.size,), dtype=bool)
                for need_role in need:
                    if need_role == "donor":
                        sat_mask &= _min_d2_all(donors_w, pxyz) <= hbd2
                    elif need_role == "acceptor":
                        sat_mask &= _min_d2_all(acceptors_w, pxyz) <= hbd2
                    elif need_role == "cation":
                        sat_mask &= _min_d2_all(cations_w, pxyz) <= ion2
                    elif need_role == "anion":
                        sat_mask &= _min_d2_all(anions_w, pxyz) <= ion2
                cov_arr |= (sat_mask.astype(np.uint16) << np.uint16(bit))

            # Early prune: only candidates that satisfy at least one polar atom can proceed.
            pass_mask = cov_arr != 0
            if not np.any(pass_mask):
                continue

            pass_local = ok_local[pass_mask]
            pass_cov = cov_arr[pass_mask].astype(np.uint16)
            pass_idx = (start + pass_local).astype(np.int64)
            pass_aa = aa_arr[pass_idx]

            # Exclude amino acids early (cheap).
            ok_aa = np.array([str(x).upper() not in exclude_aa3 for x in pass_aa], dtype=bool)
            if not np.any(ok_aa):
                continue
            pass_local = pass_local[ok_aa]
            pass_cov = pass_cov[ok_aa]
            pass_idx = pass_idx[ok_aa]

            # Step 2: for the (much smaller) satisfaction-passing set, compute full-atom world coords,
            # then apply facing + full-atom ligand clash filter.
            loc_all = center_xyz_stub[start:end][pass_local].astype(np.float64)  # (n_pass,atoms,3)
            world_all = np.einsum("nij,nkj->nki", Rw[pass_local], loc_all) + tw[pass_local][:, None, :]

            if args.require_sidechain_facing:
                if ca_i < 0 or not sidechain_idx:
                    continue

                ca_xyz = world_all[:, ca_i, :]
                ca_ok = np.isfinite(ca_xyz).all(axis=1)
                if not np.any(ca_ok):
                    continue
                world_all = world_all[ca_ok]
                pass_local = pass_local[ca_ok]
                pass_cov = pass_cov[ca_ok]
                pass_idx = pass_idx[ca_ok]

                # Prefer CA->CB; fall back to sidechain COM if CB is missing/NaN.
                v_sc = None
                cb_i = int(atom_order.index("CB")) if "CB" in atom_order else -1
                if cb_i >= 0:
                    cb_xyz = world_all[:, cb_i, :]
                    cb_ok = np.isfinite(cb_xyz).all(axis=1)
                    v_sc = np.where(cb_ok[:, None], cb_xyz - ca_xyz[ca_ok], np.nan)
                if v_sc is None:
                    v_sc = np.full_like(ca_xyz[ca_ok], np.nan)

                # For candidates without a valid CB vector, compute sidechain COM vector.
                missing = ~np.isfinite(v_sc).all(axis=1)
                if np.any(missing):
                    sc_xyz = world_all[missing][:, sidechain_idx, :]
                    sc_ok = np.isfinite(sc_xyz).all(axis=2)
                    # If no sidechain atoms are present (e.g., GLY), reject.
                    has_sc = sc_ok.any(axis=1)
                    v_tmp = np.full((int(missing.sum()), 3), np.nan, dtype=np.float64)
                    if np.any(has_sc):
                        sc_com = np.array(
                            [
                                sc_xyz[i][sc_ok[i]].mean(axis=0) if has_sc[i] else np.array([np.nan, np.nan, np.nan])
                                for i in range(sc_xyz.shape[0])
                            ],
                            dtype=np.float64,
                        )
                        v_tmp[has_sc] = sc_com[has_sc] - ca_xyz[ca_ok][missing][has_sc]
                    v_sc[missing] = v_tmp

                n_sc = np.linalg.norm(v_sc, axis=1)
                ok_sc = np.isfinite(n_sc) & (n_sc > 1e-8)
                if not np.any(ok_sc):
                    continue

                # Use the mean position of the site’s covered ligand polar atoms as the target direction.
                tgt = np.mean(np.stack([p for _an, _bit, _roles, p in site_polar], axis=0), axis=0)
                v_t = tgt[None, :] - ca_xyz[ca_ok]
                n_t = np.linalg.norm(v_t, axis=1)
                ok_t = np.isfinite(n_t) & (n_t > 1e-8)

                v_c = lig_centroid[None, :] - ca_xyz[ca_ok]
                n_c = np.linalg.norm(v_c, axis=1)
                ok_c = np.isfinite(n_c) & (n_c > 1e-8)

                ok_dir = ok_sc & ok_t & ok_c
                if not np.any(ok_dir):
                    continue

                v_sc_u = v_sc[ok_dir] / n_sc[ok_dir, None]
                v_t_u = v_t[ok_dir] / n_t[ok_dir, None]
                v_c_u = v_c[ok_dir] / n_c[ok_dir, None]

                dot_t = np.einsum("ij,ij->i", v_sc_u, v_t_u)
                dot_c = np.einsum("ij,ij->i", v_sc_u, v_c_u)

                ok_face = (dot_t >= float(args.min_sidechain_facing_dot)) & (dot_c >= float(args.min_sidechain_centroid_dot))

                mask_face = np.zeros((ok_dir.shape[0],), dtype=bool)
                mask_face[np.nonzero(ok_dir)[0]] = ok_face

                if not np.any(mask_face):
                    continue

                world_all = world_all[mask_face]
                pass_local = pass_local[mask_face]
                pass_cov = pass_cov[mask_face]
                pass_idx = pass_idx[mask_face]

            ok_full = _clash_mask_full_fast(
                world_all,
                atom_order_vdw,
                lig_xyz,
                lig_vdw,
                tolerance=float(args.clash_tol),
            )
            if not np.any(ok_full):
                continue

            pass_local = pass_local[ok_full]
            pass_cov = pass_cov[ok_full]
            pass_idx = pass_idx[ok_full]

            # Finally, score and keep top-k for this site.
            for idx, cov in zip(pass_idx.tolist(), pass_cov.tolist(), strict=True):
                sc = float(prior[idx]) + 1e-3 * int(cov).bit_count()
                maybe_add(int(idx), sc, int(cov))

        # Emit this site's kept candidates
        kept = sorted(heap, key=lambda x: (-x[0], -x[1]))
        for sc, neg_vdm, idx, cov in kept:
            cid = _stable_u64(f"{site.site_id}|{int(vdm_ids[idx])}")
            cand_id.append(np.uint64(cid))
            site_index.append(np.uint16(site.site_index))
            vdm_id.append(np.uint64(vdm_ids[idx]))
            aa3.append(str(aa_arr[idx]))
            score.append(np.float32(sc))
            cover_mask.append(np.uint16(cov))
            cluster_number.append(np.int32(cluster_number_i32[idx]))
            satisfied_union |= int(cov)

            # Compute world stub xform12 for output
            Rb = R_ifg_to_stub[idx : idx + 1]
            tb = t_ifg_to_stub[idx : idx + 1]
            Ra = np.asarray(site.R_ifg, dtype=np.float64)[None, :, :]
            ta = np.asarray(site.t_ifg, dtype=np.float64)[None, :]
            Rw, tw = _compose_xforms(Ra, ta, Rb, tb)
            xform_world_stub_12.append(_pack_xform12(Rw, tw)[0])
            center_atom_xyz_stub.append(center_xyz_stub[idx].astype(np.float32))

    if args.require_full_coverage:
        need = (1 << len(polar_atoms)) - 1
        if coverage_union != need:
            raise ValueError(
                f"Site frames do not cover all polar atoms (coverage_union={bin(coverage_union)} need={bin(need)})."
            )
        if missing_cg:
            raise FileNotFoundError(f"Missing vdXform libraries for CG(s): {sorted(missing_cg)}")

    # Deduplicate by cand_id (keep best score deterministically)
    best: dict[int, int] = {}
    for i, cid in enumerate(cand_id):
        k = int(cid)
        if k not in best:
            best[k] = i
            continue
        j = best[k]
        if float(score[i]) > float(score[j]) or (float(score[i]) == float(score[j]) and int(vdm_id[i]) < int(vdm_id[j])):
            best[k] = i

    keep_idx = sorted(best.values(), key=lambda i: (-(float(score[i])), int(cand_id[i])))

    cand_id_arr = np.array([cand_id[i] for i in keep_idx], dtype=np.uint64)
    site_idx_arr = np.array([site_index[i] for i in keep_idx], dtype=np.uint16)
    vdm_id_arr = np.array([vdm_id[i] for i in keep_idx], dtype=np.uint64)
    aa3_arr = np.array([aa3[i] for i in keep_idx], dtype="U3")
    score_arr = np.array([score[i] for i in keep_idx], dtype=np.float32)
    cover_arr = np.array([cover_mask[i] for i in keep_idx], dtype=np.uint16)
    cluster_number_arr = np.array([cluster_number[i] for i in keep_idx], dtype=np.int32)
    xw12_arr = np.stack([xform_world_stub_12[i] for i in keep_idx], axis=0).astype(np.float32)
    center_stub_arr = np.stack([center_atom_xyz_stub[i] for i in keep_idx], axis=0).astype(np.float32)

    out_npz = args.out_prefix.with_suffix(".npz")
    out_json = args.out_prefix.with_suffix(".json")
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_npz,
        cand_id_u64=cand_id_arr,
        site_index_u16=site_idx_arr,
        vdm_id_u64=vdm_id_arr,
        cluster_number_i32=cluster_number_arr,
        aa3=aa3_arr,
        score_f32=score_arr,
        cover_mask_u16=cover_arr,
        xform_world_stub_12_f32=xw12_arr,
        center_atom_xyz_stub_f32=center_stub_arr,
    )

    meta = {
        "inputs": {
            "ligand_pdb": str(args.ligand_pdb),
            "polar_sites": str(args.polar_sites),
            "site_frames": str(args.site_frames),
            "vdxform_dir": str(args.vdxform_dir),
        },
        "params": {
            "chunk_size": int(args.chunk_size),
            "top_per_site": int(args.top_per_site),
            "clash_tol": float(args.clash_tol),
        },
        "ligand": {"n_heavy_atoms": int(lig_xyz.shape[0]), "heavy_atom_names": lig_names},
        "polar_atoms": polar_atoms,
        "atom_order": atom_order,
        "n_sites": len(sites),
        "n_candidates": int(cand_id_arr.shape[0]),
        "missing_cg": sorted(missing_cg),
        "frame_coverage_union_bitmask": int(coverage_union),
        "satisfaction_union_bitmask": int(satisfied_union),
    }
    _write_json(out_json, meta)


if __name__ == "__main__":
    main()
