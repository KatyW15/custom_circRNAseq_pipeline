#!/bin/bash
#SBATCH -J bwa_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bwa_totalRNAseq_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bwa_totalRNAseq_%A_%a.err
#SBATCH -t 08:00:00    # per-sample wall time (D-HH:MM:SS)
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 16      # 16 cores -> bwa mem -t 16
#SBATCH --mem=32GB
#SBATCH --array=1-65


## 65 = total RNA-seq samples (modality "rnaseq" with real paired FASTQ paths).
## Recount if the samplesheet changes:
#   awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l

# ---------------------------------------------------------------------------
# VALIDATION branch, map step -- see circRNA_analysis_pipeline_outline.txt STEP 6.
#
# Independent circRNA detection on the RNase R(-) total RNA-seq libraries so we
# can ask which circRNA-seq consensus BSJs are ALSO seen in non-enriched data.
# Same aligner/params as the circRNA branch (bwa mem -T 19, UNSORTED SAM) so
# CIRI2 detection is comparable between the two data sets.
#
# NOTE: this is only for the orthogonal validation. The RNase R *correction*
# (run_ciriquant_totalRNAseq.sbatch) does its own alignment internally and does
# NOT need this SAM.
#
# Input : trimmed FASTQs from run_fastp_totalRNAseq.sbatch
#           ${TRIMDIR}/${SM}_R1.fq.gz , ${SM}_R2.fq.gz
# Output: ${ALIGNDIR}/${SM}.sam  (name-ordered; consumed by run_ciri2_totalRNAseq.sbatch)
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
module load samtools
#activate my circrna env with all dependancies
conda activate circrna

echo ____________________________________
which bwa
bwa 2>&1 | grep -i version || true
which samtools
samtools --version 2>&1 | head -n 1 || true
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
MODALITY_DIR="${ANALYSIS_DIR}/total_RNAseq"         # trimmed/bwa_mapped/ciri2 all nest here
TRIMDIR="${MODALITY_DIR}/trimmed"     # trimmed FASTQs from run_fastp_totalRNAseq.sbatch
ALIGNDIR="${MODALITY_DIR}/bwa_mapped" # output SAMs
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# BWA index prefix -- the .amb/.ann/.bwt/.pac/.sa files must already exist
# alongside this FASTA (built once by run_bwa_index.sbatch)
REF="/projects/$USER/HNSCC_analysis/reference_genomes/hg38/GRCh38.primary_assembly.genome.fa"

mkdir -p "${ALIGNDIR}"

###########
# Pull this task's sample from the master samplesheet
###########
# Same filter as run_fastp_totalRNAseq.sbatch, so array index N maps to the
# same sample across every total RNA step.
#   6 sample_id   7 fastq_1   8 fastq_2
line_number=$SLURM_ARRAY_TASK_ID
echo "array number ${SLURM_ARRAY_TASK_ID}"

SM=$(awk -F',' \
    'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
    "${SAMPLESHEET}" | sed -n "${line_number}p")

if [[ -z "${SM}" ]]; then
    echo "ERROR: no rnaseq row ${line_number} in ${SAMPLESHEET}" >&2
    exit 1
fi
echo "sample ${SM}"

R1="${TRIMDIR}/${SM}_R1.fq.gz"
R2="${TRIMDIR}/${SM}_R2.fq.gz"
[[ -s "${R1}" ]] || { echo "ERROR: missing trimmed R1: ${R1}" >&2; exit 1; }
[[ -s "${R2}" ]] || { echo "ERROR: missing trimmed R2: ${R2}" >&2; exit 1; }
[[ -s "${REF}.bwt" ]] || { echo "ERROR: BWA index not found for ${REF} (run run_bwa_index.sbatch)" >&2; exit 1; }

###########
# Map with BWA-MEM
###########
SAM="${ALIGNDIR}/${SM}.sam"
if [[ -s "${SAM}" ]]; then
    echo "${SM}: ${SAM} already exists, skipping bwa mem"
else
    bwa mem -T 19 -t "${THREADS}" "${REF}" \
        "${R1}" "${R2}" \
        > "${SAM}"
    echo "${SM}: bwa mem done"
fi

###########
# Alignment QC -- overall mapping rate
###########
STATSDIR="${ALIGNDIR}/flagstat"
mkdir -p "${STATSDIR}"
FLAGSTAT="${STATSDIR}/${SM}.flagstat.txt"
if [[ ! -s "${FLAGSTAT}" ]]; then
    samtools flagstat -@ "${THREADS}" "${SAM}" > "${FLAGSTAT}"
fi
echo "${SM}: flagstat ->"
grep -E 'primary mapped|properly paired' "${FLAGSTAT}" || cat "${FLAGSTAT}"

echo "${SM}: done"
