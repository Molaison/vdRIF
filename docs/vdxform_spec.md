# `vdXform` v0 (MVP) — compact vdM transform library

Purpose: make Combs2024 `vdMs/<cg>/<AA>.parquet.gzip` usable in a rifgen-style tight loop without reading parquet at runtime.

This v0 format is intentionally minimal and Python-first. It stores only what the MVP needs to enumerate residue placements around a fixed ligand:
- a **rigid transform** from ligand/iFG frame → residue backbone stub frame
- a stable id and a small amount of metadata for scoring / audit

## Naming and conventions

- **iFG frame**: a deterministic local frame built from 3 iFG atoms (chain `Y` in parquet), using the same convention as `scripts/02_polar_sites/03_build_ligand_site_frames.py`:
  - origin at atom `a1`
  - x-axis along `a1→a2`
  - y-axis = component of `a1→a0` orthogonal to x
  - z-axis = x × y
- **Stub frame**: a deterministic backbone frame computed from residue `N/CA/C`, matching rifdock `BackboneActor.from_n_ca_c` (see `external/rifdock/schemelib/scheme/actor/BackboneActor.hh`).

Transform definition:

`X_ifg_to_stub` is a 4×4 rigid transform such that:
- given ligand site frame `X_world_ifg`, we can place a residue stub by:
  - `X_world_stub = X_world_ifg * X_ifg_to_stub`

## Frame atom definitions (per cg)

Frame atoms are chosen by `configs/cg_frame_defs.json`.

Important:
- Some cg types have **multiple naming variants** (e.g., `coo` can be ASP-like `CB/CG/OD1` or GLU-like `CG/CD/OE1`).
- The converter tries alternatives in order and uses the first one whose 3 atoms exist in the instance’s chain `Y` rows.

## On-disk layout

For each cg type:

- `processed/03_vdxform/<cg>/parts/<AA>.npz` : one per input parquet file.
- `processed/03_vdxform/<cg>/manifest.json` : conversion manifest for resume/audit.
- Optional packed file:
  - `processed/03_vdxform/<cg>/vdxform_<cg>.npz`
  - `processed/03_vdxform/<cg>/vdxform_<cg>.json`

## Array schema (`*.npz`)

All arrays are aligned by row index `i`:

- `vdm_id_u64` (`uint64`, shape `(n,)`): stable id = `blake2b8(f"{cg}|{rota}|{CG}|{probe_name}|{scrn}|{resname}")`, where `scrn` is the instance residue id (`pdb_segment/pdb_chain/pdb_resnum`, or `UNK/UNK/UNK` if absent in parquet)
- `rota_i32` (`int32`, `(n,)`): parquet `rota` (cluster variant id)
- `CG_i32` (`int32`, `(n,)`): parquet `CG` (iFG naming variant / internal code)
- `aa` (`object`, `(n,)`): 3-letter AA of the placed residue (MVP: the iFG residue in chain `X`)
- `xform_ifg_to_stub_12_f32` (`float32`, `(n,12)`): packed 3×4 row-major `(R|t)` of `X_ifg_to_stub`
- `center_atom_xyz_stub_f32` (`float32`, `(n,m,3)`): center residue atom coordinates in **stub-local** frame, for atom order `m` given by `processed/03_vdxform/<cg>/atom_order.json`; missing atoms are `NaN`
- `cluster_number_i32` (`int32`, `(n,)`): parquet `cluster_number` if present, else `-1`
- `cluster_rank_ABPLE_A_f32` (`float32`, `(n,)`): parquet `cluster_rank_ABPLE_A` if present, else `NaN`
- `C_score_ABPLE_A_f32` (`float32`, `(n,)`): parquet `C_score_ABPLE_A` if present, else `NaN`

Notes:
- This v0 schema deliberately does **not** store full residue atom coords; MVP can place only stubs first.
- As the motif stage matures, we can extend v0→v1 to store minimal sidechain coords and/or an irot mapping.

Update (current repo):
- We now store **center residue atoms in stub-local coords** as `center_atom_xyz_stub_f32` to enable full-atom motif PDB output without Rosetta.
