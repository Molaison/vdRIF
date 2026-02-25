#!/usr/bin/env bash
set -euo pipefail

# Run MTX pocket-knob sweep on dg with full vdxform data, then pull summary back.
#
# Default sweep matches the historical 6-point grid:
#   ANGLE_LIST=60
#   WEIGHT_LIST=0.05,0.10,0.15
#   MIN_SC_LIST=1,2
#
# Example:
#   bash scripts/07_dg/02_run_mtx_pocket_knob_sweep_on_dg.sh
#
# Optional:
#   DG_HOST=dg
#   DG_WORKDIR=vdRIF_dg_pocket
#   DG_ENV_PREFIX=~/mambaenvs/vdrif-pocket
#   TAG=dg_grid6_$(date +%Y%m%d-%H%M%S)

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
if [[ -n "${COMMON_GIT_DIR}" ]]; then
  COMMON_ROOT="$(dirname "${COMMON_GIT_DIR}")"
else
  COMMON_ROOT="${ROOT}"
fi

DG_HOST="${DG_HOST:-dg}"
DG_WORKDIR="${DG_WORKDIR:-vdRIF_dg_pocket}"
DG_ENV_PREFIX="${DG_ENV_PREFIX:-\$HOME/mambaenvs/vdrif-pocket}"

TAG="${TAG:-dg_grid6_$(date +%Y%m%d-%H%M%S)}"
TOP_PER_SITE="${TOP_PER_SITE:-200}"
TOP_PER_SITE_PER_ATOM="${TOP_PER_SITE_PER_ATOM:-50}"
TIME_LIMIT_S="${TIME_LIMIT_S:-120}"
ANGLE_LIST="${ANGLE_LIST:-60}"
WEIGHT_LIST="${WEIGHT_LIST:-0.05,0.10,0.15}"
MIN_SC_LIST="${MIN_SC_LIST:-1,2}"
SOLVER="${SOLVER:-greedy}"

SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=12 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null)

LOCAL_LIGAND="${COMMON_ROOT}/inputs/01_cgmap/MTX.pdb"
LOCAL_POLAR="${COMMON_ROOT}/outputs/02_polar_sites/MTX_polar_sites.json"
LOCAL_FRAMES="${COMMON_ROOT}/outputs/02_polar_sites/MTX_site_frames.json"
LOCAL_VDXFORM="${COMMON_ROOT}/processed/03_vdxform_full"

for f in "${LOCAL_LIGAND}" "${LOCAL_POLAR}" "${LOCAL_FRAMES}"; do
  [[ -f "${f}" ]] || { echo "Missing required file: ${f}" 1>&2; exit 2; }
done
[[ -d "${LOCAL_VDXFORM}" ]] || { echo "Missing required dir: ${LOCAL_VDXFORM}" 1>&2; exit 2; }

echo "[dg] host=${DG_HOST}"
echo "[dg] workdir=~/${DG_WORKDIR}"
echo "[dg] env_prefix=${DG_ENV_PREFIX}"
echo "[run] tag=${TAG}"
echo "[run] angle_list=${ANGLE_LIST}"
echo "[run] weight_list=${WEIGHT_LIST}"
echo "[run] min_sc_list=${MIN_SC_LIST}"

echo "[sync] code -> dg"
rsync -az --delete -e "ssh ${SSH_OPTS[*]}" \
  --exclude ".git" \
  --exclude ".git/" \
  --exclude ".worktrees/" \
  --exclude ".venv/" \
  --exclude "external/" \
  --exclude "outputs/" \
  --exclude "processed/" \
  --exclude "logs/" \
  --exclude "__pycache__/" \
  --exclude ".pytest_cache/" \
  "${ROOT}/" "${DG_HOST}:~/${DG_WORKDIR}/"

echo "[sync] data -> dg"
ssh "${SSH_OPTS[@]}" "${DG_HOST}" "mkdir -p '${DG_WORKDIR}/inputs/01_cgmap' '${DG_WORKDIR}/outputs/02_polar_sites' '${DG_WORKDIR}/processed/03_vdxform_full'"
rsync -az -e "ssh ${SSH_OPTS[*]}" "${LOCAL_LIGAND}" "${DG_HOST}:~/${DG_WORKDIR}/inputs/01_cgmap/MTX.pdb"
rsync -az -e "ssh ${SSH_OPTS[*]}" "${LOCAL_POLAR}" "${LOCAL_FRAMES}" "${DG_HOST}:~/${DG_WORKDIR}/outputs/02_polar_sites/"
rsync -az -e "ssh ${SSH_OPTS[*]}" "${LOCAL_VDXFORM}/" "${DG_HOST}:~/${DG_WORKDIR}/processed/03_vdxform_full/"

echo "[dg] setup env + run sweep"
ssh "${SSH_OPTS[@]}" "${DG_HOST}" bash -lc "'
set -euo pipefail
source /data/apps/miniforge/25.3.1/etc/profile.d/conda.sh
ENV_PREFIX=\"${DG_ENV_PREFIX}\"
WD=\"${DG_WORKDIR}\"
TAG=\"${TAG}\"

if [[ ! -d \"\${ENV_PREFIX}\" ]]; then
  mamba create -y -p \"\${ENV_PREFIX}\" -c conda-forge \
    python=3.11 numpy=1.26 pandas pyarrow networkx rdkit pip
fi
conda run -p \"\${ENV_PREFIX}\" python -m pip install --quiet ortools==9.10.4067

cd \"\${WD}\"
PY_BIN=\"\${ENV_PREFIX}/bin/python\"
\"\${PY_BIN}\" scripts/99_harness/05_mtx_pocket_knob_sweep.py \
  --repo-root \"\$PWD\" \
  --python-bin \"\${PY_BIN}\" \
  --ligand-pdb \"\$PWD/inputs/01_cgmap/MTX.pdb\" \
  --polar-sites \"\$PWD/outputs/02_polar_sites/MTX_polar_sites.json\" \
  --site-frames \"\$PWD/outputs/02_polar_sites/MTX_site_frames.json\" \
  --vdxform-dir \"\$PWD/processed/03_vdxform_full\" \
  --solver \"${SOLVER}\" \
  --top-per-site \"${TOP_PER_SITE}\" \
  --top-per-site-per-atom \"${TOP_PER_SITE_PER_ATOM}\" \
  --time-limit-s \"${TIME_LIMIT_S}\" \
  --min-lig-donor-angle-list \"${ANGLE_LIST}\" \
  --pocket-contact-weight-list \"${WEIGHT_LIST}\" \
  --min-pocket-sidechain-contacts-list \"${MIN_SC_LIST}\" \
  --tag \"\${TAG}\"

echo \"[dg] outdir=\$PWD/processed/99_harness/mtx_pocket_knob_sweep_\${TAG}\"
'"

LOCAL_OUTDIR="${ROOT}/processed/99_harness/mtx_pocket_knob_sweep_${TAG}"
mkdir -p "${LOCAL_OUTDIR}"
echo "[sync] summary <- dg"
rsync -az -e "ssh ${SSH_OPTS[*]}" \
  "${DG_HOST}:~/${DG_WORKDIR}/processed/99_harness/mtx_pocket_knob_sweep_${TAG}/summary.json" \
  "${DG_HOST}:~/${DG_WORKDIR}/processed/99_harness/mtx_pocket_knob_sweep_${TAG}/summary.md" \
  "${LOCAL_OUTDIR}/"

echo "[done] local summary: ${LOCAL_OUTDIR}/summary.json"
