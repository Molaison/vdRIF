# vdXform → RIF export spec (draft)

Goal: turn Combs-derived vdXform libraries into a rifdock-compatible RIF (or an intermediate store that can be imported into rifdock later), while preserving:
- determinism
- streamability / resume
- the “all ligand polar atoms satisfied” hard constraint at motif-assembly time

This is intentionally written as a **spec** so we can iterate without prematurely committing to a big C++ patch.

## Inputs

1) Ligand bound pose (fixed): `inputs/01_cgmap/<LIG>.pdb`
2) Ligand polar sites: `outputs/02_polar_sites/<LIG>_polar_sites.json`
3) Deterministic ligand site frames: `outputs/02_polar_sites/<LIG>_site_frames.json`
4) vdXform libraries per CG:
   - `processed/03_vdxform/<cg>/vdxform_<cg>.npz`
   - `processed/03_vdxform/<cg>/atom_order.json`

vdXform entry semantics:
- `xform_ifg_to_stub_12_f32`: rigid transform such that `X_world_stub = X_world_ifg * X_ifg_to_stub`
- `center_atom_xyz_stub_f32`: center residue atoms in stub-local coordinates (NaN = missing)
- `cluster_number_i32`, `cluster_rank_ABPLE_A_f32`, `C_score_ABPLE_A_f32`: priors for ranking

## Output options

### Option 1: Intermediate “vdM-RIF table” (Python-first, mmap-friendly)

Store a sparse 6D hash table with:
- `bin_key`: discretized stub xform (cart + orientation bins; deterministic)
- `entries`: top‑K records, each record contains:
  - `irot_id` (bounded, ≤512/1024)
  - `score` (packed or float)
  - `sat1`, `sat2` (indices of satisfied ligand polar atoms; optional)
  - `repr_id` (optional: points back to a vdXform representative id for debugging)

File formats (choose one):
- Parquet (easy to inspect; slower random access)
- custom binary (mmap; closest to rifdock runtime needs)

### Option 2: rifdock `.rif.gz` (C++-compatible)

Export exactly rifdock’s `RifBase` format using a supported `rif_type` (e.g. `Rot10Score6Sat16`).
This requires:
- choosing an irot library ≤1024
- mapping “satisfaction” to rifdock sat group indices (not a bitmask)
- using rifdock’s XformHash discretization (or matching it byte-for-byte)

## Discretization (6D binning)

We should match rifdock’s:
- `hash_cart_resl`
- `hash_angle_resl`
- `cart_bound`

For MVP we can:
- implement the same discretization in Python (for deterministic bin keys)
- validate bin-compatibility against rifdock by spot-checking a small set of xforms

## Scoring / ranking

Within a bin, keep top‑K by:
- primary: vdM prior score (e.g. `C_score_ABPLE_A_f32` with a deterministic transform)
- secondary: number of satisfied ligand polar atoms (prefer multi-satisfying residues)
- tie-break: stable id (e.g. `vdm_id_u64` or a deterministic hash)

## Satisfaction bookkeeping (ligand polar atoms)

Our motif solver uses a `cover_mask_u16` (≤16 polar atoms).
rifdock stores up to `NSat` satisfied target group indices (`sat1`, `sat2`, ...).

Mapping proposal:
- assign each ligand polar atom an integer “sat group id” in `[0, N_polar)`
- for each residue placement, compute which polar atoms are satisfied (geometric cutoffs)
- keep the best 1–2 satisfied polar atoms as `sat1/sat2` for the RIF entry

We keep the full bitmask in the intermediate format even if exporting to rifdock later.

## Resume / determinism requirements

- all enumeration must be chunked and deterministic (sorted file order; stable hashing)
- per-bin top‑K selection must be deterministic (stable tie-break)
- writing must be atomic (tmp + rename)

## Validation

The Python pipeline (`scripts/04_candidates` + `scripts/05_solver`) remains the oracle:
- the exported RIF should be able to reproduce (or upper-bound) the candidate set found by Python enumeration
- the final motif must pass:
  - `scripts/05_solver/03_validate_motif_polar_satisfaction.py`
  - `scripts/05_solver/04_validate_motif_clashes.py`

