#!/bin/bash
# ---------------------------------------------------------------------------
# Confirm every circRNA-seq CIRIquant task (run_ciriquant_circRNA_py2.sbatch,
# step 5 / uncorrected) finished successfully, before trusting it as the
# --RNaseR input to run_ciriquant_RNaseR_corrected_py2.sbatch.
#
# Usage (run on a login node, no sbatch needed):
#   bash check_ciriquant_circRNA.sh [JOBID]
#     JOBID (optional) -- the array job id from `sbatch run_ciriquant_circRNA_py2.sbatch`.
#     If given, also checks Slurm's own exit-code record via sacct.
#
# Always runs a file-based check regardless of JOBID: for every expected
# circRNA-seq sample (master_samplesheet.csv), confirms
#     ${OUTDIR}/<sample>/<sample>.gtf
# exists, is non-empty, and has at least one circRNA quantification row.
# This catches anything sacct alone would not (e.g. a task Slurm reports as
# COMPLETED that nonetheless wrote a truncated or empty gtf).
# ---------------------------------------------------------------------------

set -uo pipefail

SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
OUTDIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/circRNAseq/ciriquant_noRNaseRcorrection_py2"
JOBID="${1:-}"

###########
# 1. Slurm accounting (if a job id was given)
###########
if [[ -n "${JOBID}" ]]; then
    echo "== sacct for job ${JOBID} =="
    sacct -j "${JOBID}" --format=JobID,State,ExitCode -X
    echo

    N_BAD=$(sacct -j "${JOBID}" --format=State -X -n | grep -vc '^COMPLETED')
    if [[ "${N_BAD}" -gt 0 ]]; then
        echo "WARNING: ${N_BAD} array task(s) did not report COMPLETED -- see above." >&2
    else
        echo "sacct: all tasks COMPLETED"
    fi
    echo
else
    echo "== no JOBID given, skipping sacct check (pass it as: bash check_ciriquant_circRNA.sh <jobid>) =="
    echo
fi

###########
# 2. File-based check: one non-empty, real gtf per expected sample
###########
mapfile -t SAMPLES < <(awk -F',' \
    'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
    "${SAMPLESHEET}")
echo "== checking ${#SAMPLES[@]} expected circRNA-seq samples in ${OUTDIR} =="

N_OK=0; N_MISSING=0; N_EMPTY=0
MISSING=(); EMPTY=()
for SM in "${SAMPLES[@]}"; do
    GTF="${OUTDIR}/${SM}/${SM}.gtf"
    if [[ ! -s "${GTF}" ]]; then
        N_MISSING=$((N_MISSING+1))
        MISSING+=("${SM}")
        continue
    fi
    N_CIRC=$(awk '$3=="circRNA"' "${GTF}" | wc -l)
    if [[ "${N_CIRC}" -eq 0 ]]; then
        N_EMPTY=$((N_EMPTY+1))
        EMPTY+=("${SM}")
    else
        N_OK=$((N_OK+1))
    fi
done

echo "  OK (gtf present, >=1 circRNA row) : ${N_OK}"
echo "  MISSING (no/empty gtf file)       : ${N_MISSING}"
echo "  EMPTY (gtf exists, 0 circRNA rows): ${N_EMPTY}"

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo
    echo "Missing samples:"
    printf '  %s\n' "${MISSING[@]}"
fi
if [[ ${#EMPTY[@]} -gt 0 ]]; then
    echo
    echo "Samples with 0 circRNA rows (check their .log in ${OUTDIR}/<sample>/):"
    printf '  %s\n' "${EMPTY[@]}"
fi

echo
if [[ "${N_MISSING}" -eq 0 && "${N_EMPTY}" -eq 0 ]]; then
    echo "RESULT: all ${#SAMPLES[@]} samples look complete -> safe to run run_ciriquant_RNaseR_corrected_py2.sbatch"
    exit 0
else
    echo "RESULT: $((N_MISSING+N_EMPTY)) of ${#SAMPLES[@]} samples are incomplete -- resolve before running the" >&2
    echo "        RNase R correction (it silently SKIPS any rnaseq sample paired to a missing gtf, exit 0, no error)." >&2
    exit 1
fi
