# Design plan — fuse Combs vdMs into rifgen-style RIF for deterministic ligand motifs

## Primary literature (anchor references)

- vdM / “van der Mer” concept and its use for small-molecule-binding site design:
  - Polizzi, DeGrado, et al. *Science* (2020) “A defined structural unit enables de novo design of small-molecule–binding proteins.” DOI: `10.1126/science.abb8330`.
- RIF / RifGen / RifDock concept (6D hashed rotamer interaction field; hierarchical/branch-and-bound search):
  - Cao, Coventry, et al. *Nature* (2022) “Design of protein-binding proteins from the target structure alone.” DOI: `10.1038/s41586-022-04654-9`.
- RIFDock used in a small-molecule-enabled protein engineering context (motivating the “motif then downstream design” workflow):
  - Langan, et al. *Nat. Biotechnol.* (2019) “De novo design of bioactive protein switches.” DOI: `10.1038/s41587-019-0260-5`.

## 0) Problem statement (frozen requirements)

Input:
- A **ligand bound pose** (protein–ligand complex structure); ligand pose treated as fixed.

Output:
- A **deterministic** motif of **8–15 residues** placed around the ligand such that:
  - **Every ligand polar atom is satisfied** (hard constraint).
  - “Satisfied” means at least one geometrically valid hydrogen bond / salt bridge partner placement (details below).

Out of scope for MVP:
- Downstream scaffold search / backbone generation / full sequence design.

Working test ligand (current session):
- `inputs/01_cgmap/MTX.pdb` (copied from user-provided path; treated as fixed pose)

Primary upstream repos:
- Combs2 / Combs2024 (`npolizzi/Combs2024`)
- rifdock / rifgen (`rifdock/rifdock`)

## 1) What “integrating the core ideas” means (operationally)

### Combs idea we keep
- Use **PDB-derived statistics**: vdM clusters represent recurring interaction geometries between a residue and a chemical group (CG).
- Use cluster enrichment / rank (e.g., `C_score_*`, `cluster_rank_*`) as a **prior** on interaction plausibility.

### RIFGen idea we keep
- Use a **6D discretized field** (position + orientation) with **fixed-capacity top-K** per voxel to compress huge candidate sets.
- Keep **satisfaction labeling** compatible with rifdock (sat1/sat2) so the generated field can be used later.

### New synthesis
Treat each vdM cluster representative as a **data-driven inverse rotamer placement**:

1) Define a local coordinate frame for each ligand polar site / functional group (CG frame).
2) For each vdM cluster rep compatible with that CG, apply a precomputed rigid transform:
   `X_world_residue_stub = X_world_ligand_CG * X_ligandCG_to_residueStub(vdM)`
3) Score and filter; then insert into a RIF voxel keyed by the stub transform.

## 1.5) Concrete integration points (where the two codebases “meet”)

rifdock / rifgen (where a future `RifGeneratorVdM` would plug in):
- `external/rifdock/apps/rosetta/rifgen.cc`: orchestrates rif generation (`generators[igen]->generate_rif(...)`, then `rif_accum->condense()`).
- `external/rifdock/apps/rosetta/riflib/rif/RifGenerator.hh`: `RifGenerator` interface + `RifAccumulator::insert(xform, score, rot, sat1, sat2, ...)`.
- `external/rifdock/apps/rosetta/riflib/rif/RifAccumulators.hh`: how inserts are buffered and merged; `condense()` merges per-thread maps into the final 6D hash map.
- `external/rifdock/apps/rosetta/riflib/RifFactory.cc` and `external/rifdock/schemelib/scheme/objective/storage/RotamerScores.hh`: RIF entry types, per-voxel top-K, and `RotamerBits` bit-packing constraints.
- `external/rifdock/apps/rosetta/riflib/RotamerGenerator.cc` and `external/rifdock/schemelib/scheme/chemical/RotamerIndex.hh`: how irot libraries are built from Rosetta rotamers (today) and where a vdM-derived irot library would eventually need to attach.

Combs2024 (where vdMs are loaded/superposed today):
- `external/Combs2024/combs2/design/superpose_ligand.py`: groups vdMs and superposes ligands onto vdM chemical-group atoms.
- `external/Combs2024/combs2/design/_sample.py`: chunked low-mem loading of vdMs + ligands; shows intended streaming patterns.
- `external/Combs2024/combs2/design/constants.py`: coordinate column conventions and other schema assumptions.

Minimal MVP stance:
- MVP does **not** require a rifdock C++ `RifGeneratorVdM` yet; we can reproduce the same conceptual flow in Python and only later port the tight loop into rifdock if needed.

## 2) Quantitative constraints (from code inspection)

### 2.1 RIF per-voxel top-K is small and fixed

From rifdock `RifFactory.cc`:
- `RotScore64`: K=28
- `RotScoreSat`: K=14
- `RotScoreSat_2x16`: K=19

So even if we enumerate 10^7 placements, the final stored field is bounded by:
`(#occupied 6D voxels) × K × (bytes per entry)`

### 2.2 rifgen already enumerates huge candidate sets

The polar generator `RifGeneratorSimpleHbonds`:
- enumerates `RelRotPos` via `make_hbond_geometries(...)`,
- and then further expands via a cartesian perturbation grid (`dx/dy/dz`) controlled by:
  - `hbond_cart_sample_hack_range`
  - `hbond_cart_sample_hack_resl`

Conclusion:
- The framework is compatible with “very many candidates”; our main problem is **front-end filtering and IO**, not the conceptual mismatch.

### 2.3 rifgen UserHotspots is random unless we remove DOFs

`RifGeneratorUserHotspots` generates `NSAMP` random perturbations seeded by `time(0)`.

For determinism:
- For MVP: enforce `hotspot_sample_cart_bound=0`, `hotspot_sample_angle_bound=0`, `hotspot_nsamples=1`.
- For production: prefer a custom deterministic generator (no RNG).

### 2.4 rifdock RIF rotamer-id bitwidth is a hard cap

In rifdock, the stored “rotamer id” inside a RIF entry is bit-packed. For the default RIF factories:
- `RotScore*` and most `RotScoreSat*`: `RotamerBits=9` ⇒ max 512 distinct ids.
- `Rot10Score6Sat16`: `RotamerBits=10` ⇒ max 1024 distinct ids.

Implication for vdMs:
- We **cannot** store 10^5–10^6 unique vdM ids directly inside a rifdock RIF without changing the data type / bit allocation.
- Any “vdM→RIF” plan needs an intermediate compression step:
  - **cluster** vdM placements into a bounded “irot” library (≤512/1024), or
  - **patch rifdock** to widen `RotamerBits` / storage types (more invasive; higher memory).

Quantitative reality check (local Combs2024 DB, MTX-relevant CGs; see `processed/06_stats/mtx_vdm_db_stats.json`):
- The existing Combs `cluster_number` already exceeds rifdock’s `RotamerBits` caps for every relevant CG:
  - `coo`: 4,182 clusters
  - `conh2`: 4,097 clusters
  - `ccn`: 2,062 clusters
  - `ph`: 4,143 clusters
  - `bb_cnh`: 1,211 clusters
So even “use cluster_number as irot id” is not viable without either (A) additional compression to ≤1024 (or ≤512), or (B) changing rifdock’s packed storage.

Practical nuance:
- This is *independent* of the per-voxel top-K: even if you only store 14–28 entries per voxel, the irot id must still fit in `RotamerBits`.

## 3) Why Combs parquet is a poor direct input for rifgen

Combs vdM DB is stored as gzipped parquet:
- `<db>/<CG>/<AA>.parquet.gzip`
- atom coordinates stored across `coords_cols` (e.g. `c_x/c_y/c_z`, `c_D_*`, `c_H*_`, `c_A*_`) and many metadata columns.
- a vdM “instance” is grouped by `['rota','CG','probe_name']` (see `combs2/design/superpose_ligand.py`), and contains:
  - rows with `chain == 'Y'`: the iFG/CG atoms used for superposition (defines `X_world_ifg`)
  - rows with `chain == 'X'`: the **interacting residue segment** (typically resnum 9/10/11) where the **interaction residue is `resnum==10`**
  - **MVP choice (this repo, 2026-01-10):** export/place **`chain=='X' & resnum==10`** for motif generation; flanks (9/11) can be added later if we decide to build backbone-contiguous motifs.
  - residue identifiers are stored as (`pdb_segment`, `pdb_chain`, `pdb_resnum`) in the parquet; Combs code often normalizes this into a single “seg/chain/resnum” key in-memory

Local database facts (user-provided path):
- Base: `~/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs`
- ~5.2 GB, 354 parquet files (18 CG dirs).

For MTX, CG typing detected 8 CG occurrences across 4 CG types (`outputs/01_cgmap/MTX_cg_atommap.json`):
- `ph`: 3 groups
- `coo`: 2 groups
- `ccn`: 2 groups
- `conh2`: 1 group

Those CG types in the local DB are already **million-scale** when counted as Combs “instances” (unique `(rota, CG, probe_name)` keys), scanned from parquet with `scripts/06_stats/01_run_mtx_vdm_db_stats.sh`:
- `coo`: 348,142 instances
- `ccn`: 145,248 instances
- `conh2`: 303,932 instances
- `ph`: 417,170 instances

Naively combining “all instances × all ligand CG occurrences” is already ~2.54M placed residues for MTX, before any pruning:
- `coo`: 348,142 × 2 ≈ 696k
- `ccn`: 145,248 × 2 ≈ 290k
- `conh2`: 303,932 × 1 ≈ 304k
- `ph`: 417,170 × 3 ≈ 1,251k

So: any MVP must be **streaming**, aggressively filtered, and resumable.

Reading parquet from C++ in rifdock is a bad fit (dependency weight + IO patterns).

Decision:

### 3.x C2 path: widen rifdock RotamerBits (chosen)

Given Combs cluster counts (e.g., `coo`/`ph` ~4k clusters), we choose the C2 approach:
- add a new rifdock `rif_type` with `RotamerBits=12` (4096 rotamer ids) and a larger packed entry type.
- current implementation in vendored rifdock: `external/rifdock/apps/rosetta/riflib/RifFactory.cc` (`Rot12ScoreSat96`).

This makes “cluster_number as irot id” feasible (up to 4096) and shifts the bottleneck to memory/performance rather than bitwidth.

Implementation note:
- In this run-log repo, `external/rifdock` is kept as an upstream submodule pointer; the C2 code change is tracked as a patch file: `external/patches/rifdock/0001-add-Rot12ScoreSat96.patch`.
- Keep parquet for **offline analytics**.
- Build a compact **vdXform** library for rifgen/rifdock consumption.

## 4) Proposed data format: vdXform (target-agnostic)

### 4.1 What rifgen really needs from vdMs
For each candidate placement we need:
- Residue identity (AA, and optionally a canonical rotamer id).
- A rigid transform placing the residue backbone stub (N/CA/C frame) in world coords.
- Which ligand polar atom(s) it satisfies (sat labels; 1–2 for compatibility).
- A score prior (vdM enrichment/rank + geometric quality).

So vdXform should store, per CG type:
- **vdXform-iFG-res (MVP)**: per-instance placement of the iFG residue (the residue owning the CG).
  - Fields: `cg_type`, `vdm_id`, `aa`, `xform_ifg_to_stub`, `score_prior`, optional cached residue coords for output/filters.
- **vdXform-3mer (future)**: rigid 3-residue microenvironment for packing/secondary filters.

### 4.2 On-disk layout
Recommended:
- One file per `cg_type` (or per cg_type × role) to keep sequential reads.
- Header includes version, endianness, and counts.
- Arrays stored as SoA (structure-of-arrays) for cache-friendly iteration.

Rationale:
- rifgen generation is a tight loop; SoA avoids per-entry parsing overhead.
- float32 is enough; RIF discretization is ~Å / degrees.

Implementation status (current repo):
- MVP v0 spec and scripts live in `docs/vdxform_spec.md` and `scripts/03_vdxform/`.
- Frame-atom definitions shared between ligand-site frames and vdXform conversion are in `configs/cg_frame_defs.json`.

## 5) Candidate generation (deterministic) — algorithm sketch

Inputs:
- Ligand pose with heavy atoms and hydrogens (or inferred).
- List of ligand polar atoms + their role (HBA/HBD/+/−).
- vdXform library.

Steps:
1) Identify ligand polar sites (hard list):
   - acceptors: O, N, S with valence rules.
   - donors: N/O with attached H.
   - ionic centers: formal charge.
2) For each polar site, map to one or more **Combs CG types** and build a local CG frame:
   - frame must be deterministic and well-defined from ligand atoms (e.g., acceptor atom + two bonded atoms).
   - practical implementation (already available locally): `~/practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py`
     - it detects Combs CGs via SMARTS patterns and enforces **connectivity-based atom ordering** so ligand atom order matches `cg_dicts.txt`
     - it can be restricted to only CGs present in a given vdM DB via `--vdm-database <.../vdMs>`
3) Enumerate candidate placements (streaming):
   - for each ligand CG occurrence, stream vdXform entries for that CG type.
   - place either:
     - a whole **3mer**: `X_world_3mer = X_world_cg * xform_cg_to_3mer`, or
     - a single residue: `X_world_stub = X_world_cg * xform_cg_to_stub`.
4) Fast filters (cheap, deterministic):
   - distance window to ligand (sanity)
   - internal residue geometry sanity (optional)
   - quick clash with ligand heavy atoms (vdW radii, no protein context in MVP)
5) Assign satisfaction set:
   - which polar atom(s) are satisfied by this placement (usually 1; sometimes 2 for bidentate).
6) Score:
   - `score = w_stat * score_prior + w_geom * hb_geom_score + w_penalty * clash_penalty`
7) Keep only top-M per polar atom (M fixed; e.g., 2000) to bound later ILP size.

Optional:
- Insert into a RIF (RotScoreSat_2x16) as a compressed cache for later rifdock usage.

## 6) Motif selection solver (deterministic; hard polar coverage)

We solve a set cover with conflicts:
- Universe U: ligand polar atoms.
- Items: candidate residue placements; each covers a subset of U.
- Conflicts: residue–residue steric clashes (precompute pairwise).
- Cardinality: 8–15 residues.
- Objective: maximize sum(scores) + optional “coverage redundancy” bonus.

Recommended solver:
- OR-Tools CP-SAT (deterministic with fixed parameters + tie-break rules).

Determinism rules:
- stable candidate ordering
- stable tie-break (cluster_id then transform hash)
- fixed solver seed/parameters if applicable

Outputs:
- `motif.pdb`: ligand + disembodied residues.
- `motif.json`: full audit trail (inputs, thresholds, solver objective, selected ids).

## 7) MVP path (minimal rifdock changes)

Because rifdock is C++/Rosetta-heavy, MVP focuses on:
- Python: read Combs parquet, build vdXform, enumerate candidates around ligand, run deterministic motif solver.
- Optional: emit a rifdock-compatible “pseudo-RIF” later once candidate generation is stable.

To still “use rifgen flow” early:
- We can generate rifgen hotspot inputs with **zero perturbation** (see determinism warning), but this is mainly a stopgap for file-compatibility testing.

## 8) Milestones / definition of done

M0 (research locked):
- Document rifgen constraints + vdM schema + vdXform design (this file + bd issues).

M1 (data layer):
- Converter: parquet → vdXform.
- Unit test: round-trip apply transform reproduces expected geometry for a tiny synthetic example.

M2 (candidate layer):
- Given a test ligand pose, produce per-polar-atom candidate sets with bounded size and stable hashes.

M3 (motif layer):
- Solve 8–15 residue motif with 100% polar coverage on 1–3 test ligands.
- Determinism test: run twice → identical outputs (byte-identical JSON; PDB stable ordering).

M4 (optional RIF export):
- Serialize candidates to a rifdock RIF type for later docking compatibility.
  - Requires resolving the `RotamerBits` cap via (a) clustering to ≤512/1024 “irots”, or (b) widening rifdock storage types.
