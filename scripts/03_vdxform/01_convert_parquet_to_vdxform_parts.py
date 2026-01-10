#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyarrow.parquet as pq


def _stable_u64(s: str) -> np.uint64:
    h = hashlib.blake2b(s.encode("utf-8"), digest_size=8).digest()
    return np.frombuffer(h, dtype="<u8")[0]


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _inv_xform(X: np.ndarray) -> np.ndarray:
    R = X[:3, :3]
    t = X[:3, 3]
    Rt = R.T
    Y = np.eye(4, dtype=np.float64)
    Y[:3, :3] = Rt
    Y[:3, 3] = -(Rt @ t)
    return Y


def _xform(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    X = np.eye(4, dtype=np.float64)
    X[:3, :3] = R
    X[:3, 3] = t
    return X


def _frame_from_points(a0: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    # Same convention as scripts/02_polar_sites/03_build_ligand_site_frames.py
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

    R = np.stack([x, y, z], axis=1)
    X = np.eye(4, dtype=np.float64)
    X[:3, :3] = R
    X[:3, 3] = a1
    return X


def _rifdock_stub_from_n_ca_c(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Match rifdock's BackboneActor.from_n_ca_c convention.

    See: external/rifdock/schemelib/scheme/actor/BackboneActor.hh
    """
    e1 = (c + n) / 2.0 - ca
    ne1 = np.linalg.norm(e1)
    if ne1 < 1e-8:
        raise ValueError("Degenerate N/CA/C: e1 too small")
    e1 = e1 / ne1

    e3 = np.cross(e1, c - ca)
    ne3 = np.linalg.norm(e3)
    if ne3 < 1e-8:
        raise ValueError("Degenerate N/CA/C: e3 too small")
    e3 = e3 / ne3

    e2 = np.cross(e3, e1)
    ne2 = np.linalg.norm(e2)
    if ne2 < 1e-8:
        raise ValueError("Degenerate N/CA/C: e2 too small")
    e2 = e2 / ne2

    R = np.stack([e1, e2, e3], axis=1)
    # rifdock uses a fixed offset in this basis to define translation t
    offset = np.array([-1.952799123558066, -0.2200069625712990, 1.524857], dtype=np.float64)
    t = (R @ offset) + ca
    return _xform(R, t)


def _np12_from_xform(X: np.ndarray) -> np.ndarray:
    """Pack 4x4 into 12 floats: row-major 3x4 (R|t)."""
    out = np.zeros((12,), dtype=np.float32)
    out[0:4] = X[0, 0:4]
    out[4:8] = X[1, 0:4]
    out[8:12] = X[2, 0:4]
    return out


ATOM_ORDER = [
    # backbone
    "N",
    "CA",
    "C",
    "O",
    # aliphatic / common
    "CB",
    "CG",
    "CG1",
    "CG2",
    "CD",
    "CD1",
    "CD2",
    "CE",
    "CE1",
    "CE2",
    "CE3",
    "CZ",
    "CZ2",
    "CZ3",
    "CH2",
    # nitrogens
    "ND1",
    "ND2",
    "NE",
    "NE1",
    "NE2",
    "NH1",
    "NH2",
    "NZ",
    # oxygens
    "OD1",
    "OD2",
    "OE1",
    "OE2",
    "OG",
    "OG1",
    "OH",
    # sulfurs
    "SD",
    "SG",
]


def _apply_inv_xform(X: np.ndarray, p: np.ndarray) -> np.ndarray:
    R = X[:3, :3]
    t = X[:3, 3]
    return R.T @ (p - t)


@dataclass(frozen=True)
class ConvertOpts:
    cg: str
    parquet_path: Path
    out_npz: Path
    frame_defs: dict[str, list[list[str]]]
    max_instances: int | None


def _iter_instances(df: Any) -> Iterable[tuple[tuple[int, int, str], Any]]:
    # group key is (rota, CG, probe_name)
    for key, g in df.groupby(["rota", "CG", "probe_name"], sort=False, dropna=False):
        yield (int(key[0]), int(key[1]), str(key[2])), g


def _convert_one_file(opts: ConvertOpts) -> dict[str, Any]:
    # Minimal columns to compute X_ifg_to_stub + a basic prior
    needed = [
        "chain",
        "CG",
        "rota",
        "probe_name",
        "resnum",
        "pdb_segment",
        "pdb_chain",
        "pdb_resnum",
        "resname",
        "resname_rota",
        "name",
        "c_x",
        "c_y",
        "c_z",
        "cluster_number",
        "cluster_rank_ABPLE_A",
        "C_score_ABPLE_A",
    ]

    pf = pq.ParquetFile(opts.parquet_path)
    cols = set(pf.schema.names)
    cols_read = [c for c in needed if c in cols]

    # NOTE: for MVP we read the whole file (per-AA) into memory, then group.
    # This is resumable per file and parallelizable across files.
    df = pf.read(columns=cols_read).to_pandas()

    # Determine which 3 atoms define the iFG frame for this cg (may have multiple variants).
    if opts.cg not in opts.frame_defs:
        raise KeyError(f"cg '{opts.cg}' missing from frame_defs")
    frame_atom_alts = opts.frame_defs[opts.cg]
    if not frame_atom_alts:
        raise ValueError(f"frame_defs[{opts.cg}] has no alternatives")

    vdm_ids: list[np.uint64] = []
    rota_v: list[np.int32] = []
    CG_v: list[np.int32] = []
    aa_v: list[str] = []
    xform12_v: list[np.ndarray] = []
    center_atoms_stub_v: list[np.ndarray] = []
    cluster_number_v: list[np.int32] = []
    cluster_rank_v: list[np.float32] = []
    c_score_v: list[np.float32] = []

    n_done = 0
    for (rota, CG, probe_name), g in _iter_instances(df):
        if opts.max_instances is not None and n_done >= opts.max_instances:
            break

        # chain Y = iFG atoms
        yg = g[g["chain"] == "Y"]
        if len(yg) == 0:
            continue
        y_by_name = {str(r["name"]): np.array([r["c_x"], r["c_y"], r["c_z"]], dtype=np.float64) for _, r in yg.iterrows()}
        X_world_ifg = None
        for alt in frame_atom_alts:
            if len(alt) != 3:
                continue
            a0n, a1n, a2n = alt
            if a0n in y_by_name and a1n in y_by_name and a2n in y_by_name:
                X_world_ifg = _frame_from_points(y_by_name[a0n], y_by_name[a1n], y_by_name[a2n])
                break
        if X_world_ifg is None:
            continue

        # chain X encodes an extended 3-residue segment (i-1,i,i+1) with `resnum` 9/10/11.
        # The interaction residue is the center residue: chain X with resnum==10.
        xg = g[g["chain"] == "X"]
        if len(xg) == 0:
            continue
        if "resnum" not in xg.columns:
            continue
        xg = xg[xg["resnum"] == 10]
        if len(xg) == 0:
            continue

        # Build residue identifier
        xg = xg.copy()
        if {"pdb_segment", "pdb_chain", "pdb_resnum"}.issubset(xg.columns):
            xg["scrn"] = xg["pdb_segment"].astype(str) + "/" + xg["pdb_chain"].astype(str) + "/" + xg["pdb_resnum"].astype(str)
        else:
            # fallback: single residue (best effort)
            xg["scrn"] = "UNK/UNK/UNK"

        # Choose residues to place: exactly the center residue (chain X, resnum==10).
        scrn_list = list(xg["scrn"].drop_duplicates())
        if not scrn_list:
            continue
        if len(scrn_list) != 1:
            # Should not happen for normal vdMs schema; keep deterministic by skipping ambiguous instances.
            continue
        scrn = scrn_list[0]

        cluster_number_val = np.int32(int(g["cluster_number"].iloc[0])) if "cluster_number" in g.columns else np.int32(-1)
        cluster_rank_val = (
            np.float32(float(g["cluster_rank_ABPLE_A"].iloc[0])) if "cluster_rank_ABPLE_A" in g.columns else np.float32(np.nan)
        )
        c_score_val = np.float32(float(g["C_score_ABPLE_A"].iloc[0])) if "C_score_ABPLE_A" in g.columns else np.float32(np.nan)

        r = xg[xg["scrn"] == scrn]
        resname = str(r["resname"].iloc[0])
        atoms = {str(rr["name"]): np.array([rr["c_x"], rr["c_y"], rr["c_z"]], dtype=np.float64) for _, rr in r.iterrows()}
        if "N" not in atoms or "CA" not in atoms or "C" not in atoms:
            continue

        X_world_stub = _rifdock_stub_from_n_ca_c(atoms["N"], atoms["CA"], atoms["C"])
        X_ifg_to_stub = _inv_xform(X_world_ifg) @ X_world_stub

        vdm_id = _stable_u64(f"{opts.cg}|{rota}|{CG}|{probe_name}|{scrn}|{resname}")
        vdm_ids.append(vdm_id)
        rota_v.append(np.int32(rota))
        CG_v.append(np.int32(CG))
        aa_v.append(resname)
        xform12_v.append(_np12_from_xform(X_ifg_to_stub))

        inv_stub = _inv_xform(X_world_stub)
        atom_local = np.full((len(ATOM_ORDER), 3), np.nan, dtype=np.float32)
        for ai, aname in enumerate(ATOM_ORDER):
            if aname in atoms:
                p = atoms[aname]
                atom_local[ai, :] = _apply_inv_xform(inv_stub, p).astype(np.float32)
        center_atoms_stub_v.append(atom_local)

        cluster_number_v.append(cluster_number_val)
        cluster_rank_v.append(cluster_rank_val)
        c_score_v.append(c_score_val)

        n_done += 1

    out = {
        "vdm_id_u64": np.array(vdm_ids, dtype=np.uint64),
        "rota_i32": np.array(rota_v, dtype=np.int32),
        "CG_i32": np.array(CG_v, dtype=np.int32),
        "aa": np.array(aa_v, dtype=object),
        "xform_ifg_to_stub_12_f32": np.stack(xform12_v, axis=0) if xform12_v else np.zeros((0, 12), dtype=np.float32),
        "center_atom_xyz_stub_f32": np.stack(center_atoms_stub_v, axis=0)
        if center_atoms_stub_v
        else np.zeros((0, len(ATOM_ORDER), 3), dtype=np.float32),
        "cluster_number_i32": np.array(cluster_number_v, dtype=np.int32),
        "cluster_rank_ABPLE_A_f32": np.array(cluster_rank_v, dtype=np.float32),
        "C_score_ABPLE_A_f32": np.array(c_score_v, dtype=np.float32),
    }

    opts.out_npz.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=str(opts.out_npz.parent),
        prefix=opts.out_npz.name + ".",
        suffix=".tmp.npz",
        delete=False,
    ) as tf:
        tmp_path = Path(tf.name)
    try:
        np.savez_compressed(tmp_path, **out)
        os.replace(tmp_path, opts.out_npz)
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    return {
        "parquet": str(opts.parquet_path),
        "out_npz": str(opts.out_npz),
        "cg": opts.cg,
        "n_entries": int(out["vdm_id_u64"].shape[0]),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Combs2024 vdM parquet into vdXform parts (per AA file).")
    ap.add_argument("--vdm-db", type=Path, required=True, help="Path to Combs vdM db root (the directory containing CG folders).")
    ap.add_argument("--cg", type=str, required=True, help="CG folder name (e.g., coo, ph, conh2, bb_cnh).")
    ap.add_argument("--outdir", type=Path, required=True, help="Output dir for vdXform parts.")
    ap.add_argument("--frame-defs", type=Path, default=Path("configs/cg_frame_defs.json"))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--max-instances-per-file", type=int, default=0, help="For debugging: cap instances per AA file (0=all).")
    args = ap.parse_args()

    frame_defs = _load_json(args.frame_defs)

    cgdir = args.vdm_db / args.cg
    if not cgdir.exists():
        raise FileNotFoundError(f"CG folder not found: {cgdir}")
    in_files = sorted(cgdir.glob("*.parquet.gzip"))
    if not in_files:
        raise FileNotFoundError(f"No parquet files found in: {cgdir}")

    out_parts = args.outdir / args.cg / "parts"
    out_parts.mkdir(parents=True, exist_ok=True)

    # Record atom order once per cg for downstream tools.
    atom_order_path = args.outdir / args.cg / "atom_order.json"
    if not atom_order_path.exists():
        _write_json(atom_order_path, {"atom_order": ATOM_ORDER})

    max_instances = None if args.max_instances_per_file <= 0 else int(args.max_instances_per_file)

    jobs: list[ConvertOpts] = []
    for pf in in_files:
        out_npz = out_parts / (pf.name.replace(".parquet.gzip", "") + ".npz")
        if args.resume and out_npz.exists():
            continue
        jobs.append(
            ConvertOpts(
                cg=args.cg,
                parquet_path=pf,
                out_npz=out_npz,
                frame_defs=frame_defs,
                max_instances=max_instances,
            )
        )

    if not jobs:
        _write_json(args.outdir / args.cg / "manifest.json", {"cg": args.cg, "n_jobs": 0, "status": "nothing_to_do"})
        return

    if args.workers == 1:
        results = [_convert_one_file(j) for j in jobs]
    else:
        with Pool(processes=args.workers) as pool:
            results = list(pool.imap_unordered(_convert_one_file, jobs))

    manifest = {
        "cg": args.cg,
        "vdm_db": str(args.vdm_db),
        "frame_defs": str(args.frame_defs),
        "n_files_total": len(in_files),
        "n_files_converted": len(results),
        "results": sorted(results, key=lambda r: r["parquet"]),
    }
    _write_json(args.outdir / args.cg / "manifest.json", manifest)


if __name__ == "__main__":
    main()
