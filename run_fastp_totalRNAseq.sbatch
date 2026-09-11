#!/bin/bash
#SBATCH -J fastp_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastp_totalRNAseq_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastp_totalRNAseq_%A_%a.err
#SBATCH -t 06:00:00    # per-sample wall time (D-HH:MM:SS)
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> fastp -w 8
#SBATCH --mem=16GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=YOUR_EMAIL@example.com
#SBATCH --array=1-65


## 65 = total RNA-seq samples (modality "rnaseq" with real paired FASTQ paths).
## Recount if the samplesheet changes:
#   awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l

# ---------------------------------------------------------------------------
# STEP 6 total RNA-seq branch -- trim.
# Same fastp settings as the circRNA branch (run_fastp_circRNA.sbatch) so the
# two data sets are processed identically before detection / quantification.
# Outputs -> trimmed/total_RNAseq/ and QC/total_RNAseq/trimmed/.
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge

#activate my circrna env with all dependancies
conda activate circrna

echo ____________________________________
which fastp
fastp --version
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
TRIMDIR="${ANALYSIS_DIR}/trimmed/total_RNAseq"      # trimmed FASTQs
QCDIR="${ANALYSIS_DIR}/QC/total_RNAseq/trimmed"     # fastp json/html reports
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${TRIMDIR}" "${QCDIR}"

###########
# Pull this task's sample straight from the master samplesheet
###########
# Same approach as run_fastp_circRNA.sbatch, filter changed to modality "rnaseq".
#   6 sample_id   7 fastq_1   8 fastq_2
line_number=$SLURM_ARRAY_TASK_ID
echo "array number ${SLURM_ARRAY_TASK_ID}"

ROW=$(awk -F',' -v OFS='\t' \
    'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6, $7, $8}' \
    "${SAMPLESHEET}" | sed -n "${line_number}p")

if [[ -z "${ROW}" ]]; then
    echo "ERROR: no rnaseq row ${line_number} in ${SAMPLESHEET}" >&2
    exit 1
fi
IFS=$'\t' read -r SM R1 R2 <<< "${ROW}"
echo "sample ${SM}"
echo "  R1: ${R1}"
echo "  R2: ${R2}"

###########
# Trim with fastp
###########
if [[ -s "${TRIMDIR}/${SM}_R1.fq.gz" && -s "${TRIMDIR}/${SM}_R2.fq.gz" ]]; then
    echo "${SM}: trimmed FASTQs already exist, skipping fastp"
else
    fastp \
        -i "${R1}" -I "${R2}" \
        -o "${TRIMDIR}/${SM}_R1.fq.gz" -O "${TRIMDIR}/${SM}_R2.fq.gz" \
        --detect_adapter_for_pe \
        -q 20 -u 30 \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
        -l 50 \
        -w "${THREADS}" \
        -j "${QCDIR}/${SM}.fastp.json" -h "${QCDIR}/${SM}.fastp.html"
fi

echo "${SM}: done"
