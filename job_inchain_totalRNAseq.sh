#!/bin/bash
#
# Submit the total RNA-seq pre-processing chain as dependent Slurm jobs:
#
#     trim --> bwa --> ciri2 --> fastQC          (fastQC:  afterany ciri2)
#              \--> flagstat                      (flagstat: afterany bwa)
#
# The array stages (trim, bwa, ciri2) are chained with --dependency=aftercorr:
# task N of each waits only on task N of the previous, so one bad sample drops
# out of the rest of the chain instead of blocking every sample. The two QC jobs
# are single jobs that run once their input stage finishes, each with afterany so
# a failed sample does not block them:
#   fastQC   -- raw-read QC, afterany on ciri2 (the last array stage)
#   flagstat -- run_flagstat_totalRNAseq.sbatch: samtools flagstat over the BWA
#               SAMs + a mapping-rate summary table; afterany on bwa, so the
#               summary is ready as soon as mapping finishes (not held behind
#               the slow ciri2 step).
# CIRIquant / RNase R correction (run_ciriquant_totalRNAseq.sbatch) and the
# strandedness check are NOT in this chain -- run those separately.
#
# This script just calls sbatch a few times and exits in seconds. Run it on a
# login node; once it returns the queued jobs are independent of your session,
# so you can log off.
#
#     bash job_inchain_totalRNAseq.sh
#
# Based on my_circRNAseq_pipeline/job_inchain.sh.
#

set -euo pipefail

# --export=NONE: child jobs start from a clean environment (each *_totalRNAseq
# .sbatch does its own module load / conda activate) -- keeps the submitting
# shell's env from leaking in.
SUBMIT="sbatch --export=NONE"

# Work from the dir holding the *_totalRNAseq.sbatch scripts (this file's dir).
CHAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${CHAIN_DIR}"
echo "chain dir: ${CHAIN_DIR}"

for s in run_fastp_totalRNAseq.sbatch run_bwa_totalRNAseq.sbatch \
         run_ciri2_totalRNAseq.sbatch run_fastqc_totalRNAseq.sbatch \
         run_flagstat_totalRNAseq.sbatch; do
    [[ -f "${s}" ]] || { echo "ERROR: ${s} not found in ${CHAIN_DIR}" >&2; exit 1; }
done

SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"

# One array task per total RNA-seq sample (modality "rnaseq" with real FASTQ paths)
N=$(awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
    "${SAMPLESHEET}" | wc -l)
echo "total RNA-seq samples: ${N}"
[[ "${N}" -gt 0 ]] || { echo "ERROR: no rnaseq samples found in ${SAMPLESHEET}" >&2; exit 1; }

# Phase 1 -- trim
tmp1=$(${SUBMIT} --array=1-${N} run_fastp_totalRNAseq.sbatch)
jid1=${tmp1##* }
echo "trim   : ${jid1}"

# Phase 2 -- BWA-MEM  (task N after trim task N)
tmp2=$(${SUBMIT} --array=1-${N} --dependency=aftercorr:${jid1} run_bwa_totalRNAseq.sbatch)
jid2=${tmp2##* }
echo "bwa    : ${jid2}  (aftercorr ${jid1})"

# Phase 3 -- CIRI2  (task N after bwa task N)
tmp3=$(${SUBMIT} --array=1-${N} --dependency=aftercorr:${jid2} run_ciri2_totalRNAseq.sbatch)
jid3=${tmp3##* }
echo "ciri2  : ${jid3}  (aftercorr ${jid2})"

# Phase 4 -- raw-read FastQC (single job, independent; not blocked by a failed sample)
tmp4=$(${SUBMIT} --dependency=afterany:${jid3} run_fastqc_totalRNAseq.sbatch)
jid4=${tmp4##* }
echo "fastqc : ${jid4}  (afterany ${jid3})"

# Phase 5 -- mapping QC: samtools flagstat over the BWA SAMs + summary table.
# Only needs the SAMs (bwa), so depends on that; afterany -> profiles whatever mapped.
tmp5=$(${SUBMIT} --dependency=afterany:${jid2} run_flagstat_totalRNAseq.sbatch)
jid5=${tmp5##* }
echo "flagstat: ${jid5}  (afterany ${jid2})"

echo
echo "chain submitted. watch with:  squeue -u \$USER -o '%.18i %.20j %.10T %.10M %R'"
echo "with aftercorr, a failed sample only drops its own downstream tasks."
echo "to stop the whole chain:  scancel ${jid1} ${jid2} ${jid3} ${jid4} ${jid5}"
