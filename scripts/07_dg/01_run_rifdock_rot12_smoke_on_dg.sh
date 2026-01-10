#!/usr/bin/env bash
set -euo pipefail

# Runs the rifdock build + rifgen smoke test on the `dg` server.
#
# Requirements on dg:
# - a Rosetta checkout exists (you have a license) and you know its path to `rosetta/main`.
#
# Usage:
#   ROSETTA_MAIN=/path/to/rosetta/main bash scripts/07_dg/01_run_rifdock_rot12_smoke_on_dg.sh
#
# Optional:
#   DG_HOST=dg
#   DG_WORKDIR=vdRIF_dg_smoke
#   REPO_URL=git@github.com:Molaison/vdRIF.git
#   DG_USE_RSYNC=1   # sync local repo to dg instead of git clone/pull (useful if dg can't reach GitHub)

DG_HOST="${DG_HOST:-dg}"
DG_WORKDIR="${DG_WORKDIR:-vdRIF_dg_smoke}"
REPO_URL="${REPO_URL:-git@github.com:Molaison/vdRIF.git}"
ROSETTA_MAIN="${ROSETTA_MAIN:-}"
DG_USE_RSYNC="${DG_USE_RSYNC:-0}"
DG_SUBMODULE_RECURSIVE="${DG_SUBMODULE_RECURSIVE:-0}"
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=10)

if [[ -z "${ROSETTA_MAIN}" ]]; then
  cat <<'EOF' 1>&2
ROSETTA_MAIN is required (path to rosetta/main on dg).
Example:
  ROSETTA_MAIN=~/rosetta/main bash scripts/07_dg/01_run_rifdock_rot12_smoke_on_dg.sh
EOF
  exit 2
fi

echo "[dg] host=${DG_HOST}"
echo "[dg] workdir=${DG_WORKDIR} (remote $HOME/${DG_WORKDIR})"
echo "[dg] repo_url=${REPO_URL}"
echo "[dg] rosetta_main=${ROSETTA_MAIN}"
echo "[dg] use_rsync=${DG_USE_RSYNC}"
echo "[dg] submodule_recursive=${DG_SUBMODULE_RECURSIVE}"

if [[ "${DG_USE_RSYNC}" == "1" ]]; then
  echo "[dg] rsync local repo -> ${DG_HOST}:~/${DG_WORKDIR}"
  ssh "${SSH_OPTS[@]}" "${DG_HOST}" "mkdir -p '${DG_WORKDIR}'"
  rsync -az --delete \
    --exclude 'outputs/' \
    --exclude 'processed/' \
    --exclude 'logs/' \
    --exclude '.venv/' \
    --exclude '__pycache__/' \
    --exclude '.pytest_cache/' \
    --exclude '.mypy_cache/' \
    --exclude 'external/rifdock/build/' \
    ./ "${DG_HOST}:~/${DG_WORKDIR}/"
fi

ssh "${SSH_OPTS[@]}" "${DG_HOST}" bash -lc "'
set -euo pipefail

WD=\"${DG_WORKDIR}\"
ROSETTA_MAIN=\"${ROSETTA_MAIN}\"
REPO_URL=\"${REPO_URL}\"

echo \"[dg] user=\$(whoami) host=\$(hostname) cwd=\$(pwd)\"

if [[ ! -d \"\$ROSETTA_MAIN\" ]]; then
  echo \"ROSETTA_MAIN does not exist on dg: \$ROSETTA_MAIN\" 1>&2
  exit 10
fi
if [[ ! -d \"\$ROSETTA_MAIN/database\" ]]; then
  echo \"Rosetta database not found at: \$ROSETTA_MAIN/database\" 1>&2
  exit 11
fi

if [[ -d \"\$WD/.git\" ]]; then
  echo \"[dg] using existing checkout: \$WD\"
  cd \"\$WD\"
  if [[ \"${DG_USE_RSYNC}\" != \"1\" ]]; then
    echo \"[dg] git fetch/pull (DG_USE_RSYNC!=1)\"
    git fetch --all --prune
    git checkout main || git checkout master
    git pull --rebase
  else
    echo \"[dg] skipping git fetch/pull (DG_USE_RSYNC=1)\"
    git status --porcelain || true
  fi
else
  echo \"[dg] cloning: \$REPO_URL -> \$WD\"
  rm -rf \"\$WD\"
  git clone \"\$REPO_URL\" \"\$WD\"
  cd \"\$WD\"
fi

# Rewrite https GitHub URLs (submodules) to SSH on dg.
git config url.\"git@github.com:\".insteadOf https://github.com/

if [[ "${DG_SUBMODULE_RECURSIVE}" == "1" ]]; then
  git submodule update --init --recursive
else
  git submodule update --init
fi

PATCH=\"external/patches/rifdock/0001-add-Rot12ScoreSat96.patch\"
echo \"[dg] applying rifdock patch: \$PATCH\"
git -C external/rifdock apply --check \"../patches/rifdock/0001-add-Rot12ScoreSat96.patch\" 2>/dev/null && \\
  git -C external/rifdock apply \"../patches/rifdock/0001-add-Rot12ScoreSat96.patch\" || \\
  echo \"[dg] patch already applied (or conflicts); continuing\"

mkdir -p external/rifdock/build
cd external/rifdock/build

export CMAKE_ROSETTA_PATH=\"\$ROSETTA_MAIN\"
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j \$(nproc) rifgen

cd \"\$WD\"
mkdir -p outputs/07_rifgen_smoke

RIFGEN=\"external/rifdock/build/apps/rosetta/rifgen\"
if [[ ! -x \"\$RIFGEN\" ]]; then
  echo \"Missing rifgen binary at \$RIFGEN\" 1>&2
  exit 12
fi

echo \"[dg] running rifgen smoke (Rot12ScoreSat96)\"
\"\$RIFGEN\" -database \"\$ROSETTA_MAIN/database\" @configs/rifgen_smoke_rot12.flags

OUT=\"outputs/07_rifgen_smoke/rot12_smoke.rif.gz\"
if [[ ! -s \"\$OUT\" ]]; then
  echo \"Smoke test FAILED: expected non-empty \$OUT\" 1>&2
  exit 13
fi
echo \"[dg] rifgen wrote: \$OUT\"
ls -lh \"\$OUT\"

echo \"[dg] append-mode load/save check\"
\"\$RIFGEN\" -database \"\$ROSETTA_MAIN/database\" @configs/rifgen_smoke_rot12.flags -rifgen:rif_append_mode true -rifgen:append_mode_clear_sats true

echo \"[dg] OK\"
'"
