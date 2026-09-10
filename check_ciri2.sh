#!/bin/bash
# ---------------------------------------------------------------------------
# Locate CIRI2.pl before submitting run_ciri2_circRNA.sbatch.
#
# Run interactively on a login node (no sbatch needed):
#     bash check_ciri2.sh
#
# On success it prints the export line to paste before `sbatch`, e.g.
#     export CIRI2_PL=/.../envs/circrna/bin/CIRI2.pl
# If CIRI2.pl is not found it lists what to install / where else to look.
# ---------------------------------------------------------------------------

set -uo pipefail

module purge 2>/dev/null || true
module load miniforge 2>/dev/null || true
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate circrna

echo "==> conda env : ${CONDA_DEFAULT_ENV:-<none>}"
echo "==> perl      : $(command -v perl || echo 'MISSING') ($(perl -e 'print $^V' 2>/dev/null))"
echo

# 1. On PATH (bioconda 'ciri2' package installs it here)
CIRI2_PL="$(command -v CIRI2.pl || true)"

# 2. Anywhere under the active conda env
if [[ -z "${CIRI2_PL}" && -n "${CONDA_PREFIX:-}" ]]; then
    CIRI2_PL="$(find "${CONDA_PREFIX}" -name 'CIRI2.pl' -type f 2>/dev/null | head -n1)"
fi

# 3. Bundled inside a CIRIquant install
if [[ -z "${CIRI2_PL}" ]]; then
    CIRIQUANT_BIN="$(command -v CIRIquant || true)"
    if [[ -n "${CIRIQUANT_BIN}" ]]; then
        CIRI2_PL="$(find "$(dirname "${CIRIQUANT_BIN}")/.." -name 'CIRI2.pl' -type f 2>/dev/null | head -n1)"
    fi
fi

echo "-------------------------------------------------------------"
if [[ -z "${CIRI2_PL}" ]]; then
    echo "RESULT: CIRI2.pl NOT found."
    echo
    echo "Fix (in the circrna env):"
    echo "    conda install -n circrna -c bioconda ciri2"
    echo "or point CIRI2_PL at a manual copy / the one inside CIRIquant:"
    echo "    find / -name CIRI2.pl 2>/dev/null"
    exit 1
fi

echo "RESULT: found CIRI2.pl"
echo "    path : ${CIRI2_PL}"
perl "${CIRI2_PL}" 2>&1 | grep -i -m1 'version\|CIRI2' || true

# Confirm perl can actually load the script (catches missing Perl deps)
if perl -c "${CIRI2_PL}" >/dev/null 2>&1; then
    echo "    perl -c: OK (script compiles)"
else
    echo "    perl -c: WARNING -- script does not compile cleanly:"
    perl -c "${CIRI2_PL}" 2>&1 | sed 's/^/        /'
fi

echo
echo "Use it with:"
echo "    export CIRI2_PL=\"${CIRI2_PL}\""
echo "    sbatch run_ciri2_circRNA.sbatch"
