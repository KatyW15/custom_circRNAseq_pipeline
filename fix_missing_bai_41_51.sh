#!/bin/bash
# ---------------------------------------------------------------------------
# Fix for array tasks 41 (TP00949_tumor_rnaseq) and 51 (TP00973-02_tumor_rnaseq)
# from run_ciriquant_RNaseR_corrected_py2.sbatch (job 32657652): both died with
#   OSError: [Errno 2] No such file or directory: '.../align/<sample>.sorted.bam.bai'
# CIRIquant found a leftover .sorted.bam from an earlier interrupted attempt,
# skipped re-alignment, then choked looking for the .bai that was never built.
# Fix: remove the stale per-sample output dir so CIRIquant redoes alignment +
# indexing cleanly, then resubmit just these two array indices.
#
# Run on the Alpine login node (not locally -- /scratch/alpine is not mounted
# on this machine).
# ---------------------------------------------------------------------------

set -euo pipefail

OUTDIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/ciriquant_py2"
SAMPLES=("TP00949_tumor_rnaseq" "TP00973-02_tumor_rnaseq")

for SM in "${SAMPLES[@]}"; do
    DIR="${OUTDIR}/${SM}"
    if [[ -d "${DIR}" ]]; then
        echo "removing stale output dir: ${DIR}"
        rm -rf "${DIR}"
    else
        echo "no existing output dir for ${SM} (nothing to remove)"
    fi
done

echo "resubmitting array indices 41,51"
sbatch --array=41,51 run_ciriquant_RNaseR_corrected_py2.sbatch
