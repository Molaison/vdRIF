# Takeaways (notes worth keeping)

## 1) RIF is already designed to compress huge enumeration

In rifdock, a RIF voxel stores a **fixed-capacity top-K** list of rotamer scores (`RotamerScores<N,...>`). This means generation can enumerate enormous candidate placements, but the final memory scales with the number of occupied voxels times `N`, not with total candidates.

From `external/rifdock/apps/rosetta/riflib/RifFactory.cc` (commit `62c3a1a...`):

- `RotScore64`: `RotamerScores<28, RotamerScore<...>>` (K=28)
- `RotScoreSat`: `RotamerScores<14, RotamerScoreSat<uint16_t,9,-4>>` (K=14, 2 sat slots of `uint8_t`)
- `RotScoreSat_1x16`: `RotamerScores<14, RotamerScoreSat<uint16_t,9,-13, SatDatum<uint16_t>, 1>>` (K=14, 1 sat slot of `uint16_t`)
- `RotScoreSat_2x16`: `RotamerScores<19, RotamerScoreSat<uint16_t,9,-13, SatDatum<uint16_t>, 2>>` (K=19, 2 sat slots of `uint16_t`)
- `RotScoreReq`: `RotamerScores<28, RotamerScoreSat<..., SatDatum<uint8_t>, 1>>` (K=28, “requirements” style)

Implication for vdM integration:
- We **can** start from millions of vdM-derived placements if we (a) bin them into the same 6D hash as rifgen and (b) rely on top-K truncation per voxel.
- The bottleneck shifts to **insertion throughput** (hash-map churn) and **front-end candidate filtering**.

## 2) rifgen “polar” enumeration is already very large-scale

The polar RIF generator (`RifGeneratorSimpleHbonds`) does not just loop over a fixed small rotamer set.

Key facts:
- It builds “HB geometry samples” (`RelRotPos`) using `make_hbond_geometries(...)` and caches them on disk (`.rel_rot_pos.gz`).
- It then iterates over every `RelRotPos` and, for each, may apply a **cartesian perturbation grid** via `hbond_cart_sample_hack_range/resl`, i.e. nested `dx/dy/dz` loops.

So the existing rifgen workflow is compatible with “many candidates” in spirit; your vdM library being large is not, by itself, a blocker.

## 3) rifgen UserHotspots is random by default (determinism warning)

`RifGeneratorUserHotspots` seeds `std::mt19937` with `(unsigned int)time(0) + 934875` and generates `NSAMP` random perturbations.

For a deterministic pipeline, do one of:
- Set `hotspot_sample_cart_bound=0`, `hotspot_sample_angle_bound=0`, `hotspot_nsamples=1` (practical workaround; removes degrees of freedom).
- Or implement a dedicated deterministic generator (`RifGeneratorVdM`) that never calls RNG.

## 4) Combs2024 vdM DB schema is rich but not C++-friendly

Combs stores vdMs as gzipped parquet files per chemical group and amino acid:

- Layout: `<db>/<CG>/<AA>.parquet.gzip`
- A vdM instance is keyed by `['rota','CG','probe_name']` (see `external/Combs2024/combs2/design/superpose_ligand.py`).
- Coordinates are spread across `coords_cols` such as `c_x/c_y/c_z`, `c_D_*`, `c_H*_*`, `c_A*_*` (see `external/Combs2024/combs2/design/constants.py`).
- Cluster/enrichment metadata include `cluster_number`, `cluster_rank_*`, `C_score_*`, `centroid*`, `maxdist_to_centroid`, etc.
- The database is distributed as a tarball (~7.5GB at the time of writing) per Combs2024 README.

Implication:
- For rifgen integration, we should **not** try to read parquet in C++.
- Instead: a one-time conversion into a compact binary “vdXform” library with exactly what rifgen needs (rigid transforms + residue identity + scores).

Local quantitative snapshot (your current DB at `~/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs`):
- Total size: ~5.2GB
- Total parquet files: 354 (`18` CG folders; mostly `20` AA files per CG, `isopropyl` has 14)
- Largest CG folders by disk: `bb_cco` (~534MB), `isopropyl` (~509MB), `ph` (~460MB), `phenol` (~435MB), `coo` (~405MB)

## 5) “All polar atoms must be satisfied” implies a hard-constraint solver

RIFDock’s satisfaction fields (`sat1/sat2`) are useful, but:
- a residue placement can cover more than 2 polar atoms,
- and the final motif must satisfy **all** ligand polar atoms with **8–15** residues.

So we should plan a deterministic combinatorial selection step (ILP/CP-SAT) on top of candidate placements:
- Hard constraints: cover every polar atom; no clashes; size 8–15.
- Objective: maximize vdM statistical score + geometric quality + (optional) buriedness.

Practical env note:
- This workspace’s system `python3` is 3.6.8, and the system RDKit is broken (missing FreeINCHI symbols).
- Use the local `uv` environment for anything that needs `pyarrow` / `ortools` / RDKit.
- Known pitfall: `rdkit-pypi==2022.9.5` requires `numpy<2` (wheels are compiled against NumPy 1.x).
- Working commands:
  - `uv sync -p 3.11`
  - `uv sync -p 3.11 --extra rdkit`

## 6) cg_atommap may miss some polar donors; add deterministic fallback

For MTX (`inputs/01_cgmap/MTX.pdb`), Combs CG mapping covers most polar atoms but missed `NA2` and `NA4` as donors in our RDKit-derived polar list.\n
\n
To avoid “uncovered polar atoms” blocking the hard-coverage constraint, `scripts/02_polar_sites/03_build_ligand_site_frames.py` adds an optional deterministic fallback:\n
- For each uncovered donor atom `N*`, create a `bb_cnh` site frame using ligand atoms `(CA,N,H)=(heavy_neighbor_of_N, N, smallest_H_attached_to_N)`.\n
- This leverages the Combs `bb_cnh` vdM library (backbone N–H) to potentially satisfy that donor.\n
\n
This is an MVP hack to keep the pipeline moving; later we should prefer expanding CG typing so donors/acceptors map to the intended sidechain CG types, not just backbone.

## 7) CP-SAT must be made deterministic and must not “fill to max residues”

When many candidates have identical (or near-identical) scores, a naive maximize objective will happily pick extra residues (up to the `max_res` bound) because any positive tie-break weight makes “select more” better.\n
\n
In `scripts/05_solver/01_solve_motif.py`, the objective is encoded lexicographically as a *single* linear maximize:\n
- primary: maximize total score\n
- secondary: minimize residue count\n
- tertiary: deterministic tie-break (candidate rank)\n
\n
Also, CP-SAT parallelism can introduce nondeterminism; for reproducibility we run with `--num-workers 1` and a fixed seed.

## 8) “Coverage” must mean geometric satisfaction (and NaNs can silently break it)

Early iterations mistakenly treated “covered” as “this polar atom belongs to this site frame”, which allowed the solver to pick a motif that *looked* complete but did not satisfy polar atoms in 3D.

Two concrete MVP fixes that made the MTX debug run actually pass validation:
- In `scripts/04_candidates/01_generate_candidates.py`, compute `cover_mask_u16` from coarse **distance-based satisfaction** using the placed residue’s full-atom coordinates (not from site-frame membership).
- When atoms are missing in `center_atom_xyz_stub_f32` they are stored as `NaN`; distance code must treat those as “not present”. Using `np.nansum` (or `nanmin` without guarding all-NaN) can create false positives.

## 9) Symmetric ligand CGs need “swap frames” (OD1/OD2, OE1/OE2)

Combs vdM iFG frames are defined with *labeled* atoms (e.g., carboxylate uses `.../OD1` in `configs/cg_frame_defs.json`). For amino acids this labeling is consistent, but ligand PDB naming like `O1/O2` is arbitrary.

If we align ligand site frames using only one oxygen (e.g., `OD1→O1`), we can accidentally mirror the iFG frame and misplace vdMs.

MVP fix:
- In `scripts/02_polar_sites/03_build_ligand_site_frames.py`, for `cg=="coo"` we emit an extra `cgmap:<key>:swap_*` site frame that swaps `OD1↔OD2` (and `OE1↔OE2`) when those labels exist in the `correspond_names`.

## 10) “All residues look like C-C(-N)-C” is usually backbone-only output

If candidates are missing `center_atom_xyz_stub_f32` (or meta is missing `atom_order`), the solver cannot reconstruct full-atom sidechains and the motif will look like a generic 4-atom backbone stub in many viewers.

Current behavior:
- `scripts/99_harness/01_mtx_determinism_regression.py` treats this as an error unless `--allow-backbone-only` is passed.
- `scripts/05_solver/01_solve_motif.py` now raises if full-atom reconstruction inputs are missing (prevents silent “everything looks like Ala” outputs).

## 11) Sidechain-only satisfaction reduces “hard-to-realize” backbone hbonds

If we allow backbone `N` (donor) / backbone `O` (acceptor) to satisfy ligand polar atoms, then **almost any residue** can satisfy constraints via its backbone, often producing motifs whose sidechains do not point into the pocket.

MVP workaround:
- `scripts/04_candidates/01_generate_candidates.py` defaults to *sidechain-only* satisfaction by excluding backbone `N/O` from donor/acceptor atom sets.
- You can re-enable backbone hbonds with `--allow-backbone-hbonds` if desired.

## 12) Excluding amino acids early (PRO/CYS) is easiest at candidate generation

For ligand-centric polar motifs, some amino acids are often undesirable for practical reasons:
- `PRO` is structurally constrained and often hard to accommodate in a pocket as a designed motif residue.
- `CYS` can complicate expression/handling and downstream design choices.

MVP:
- `scripts/04_candidates/01_generate_candidates.py` adds `--exclude-aa3` (default `PRO,CYS`).
