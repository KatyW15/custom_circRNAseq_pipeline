#!/bin/bash
#SBATCH -J fastp_circRNA
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastp_circRNA_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastp_circRNA_%A_%a.err
#SBATCH -t 06:00:00    # per-sample wall time (D-HH:MM:SS)
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> fastp -w 8
#SBATCH --mem=16GB
#SBATCH --array=1-57

## hard coded array above, if want to define array in command line see below
## Submit one array task per circRNA sample. Get the count (= array range) with:
#   awk -F',' 'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l
#   sbatch --array=1-57 run_fastp_circRNA.sbatch

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
MODALITY_DIR="${ANALYSIS_DIR}/circRNAseq"       # trimmed/bwa_mapped/ciri2 all nest here
TRIMDIR="${MODALITY_DIR}/trimmed"      # trimmed FASTQs
QCDIR="/projects/$USER/HNSCC_analysis/my_circRNAseq_analysis/QC/circRNAseq/trimmed"   # fastp json/html reports
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${TRIMDIR}" "${QCDIR}"

###########
# Pull this task's sample straight from the master samplesheet
###########
# master_samplesheet.csv columns:
#   1 patient_id  2 batch_id  3 tissue  4 hpv_status  5 modality
#   6 sample_id   7 fastq_1   8 fastq_2  9 survival_time_days  10 vital_status  11 notes
# Filter to circRNA rows with real paired FASTQ paths 
# (drops rows with placeholder text like "NA")
# then the array number picks the Nth such row -- no intermediate file needed.
line_number=$SLURM_ARRAY_TASK_ID
echo "array number ${SLURM_ARRAY_TASK_ID}"

ROW=$(awk -F',' -v OFS='\t' \
    'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6, $7, $8}' \
    "${SAMPLESHEET}" | sed -n "${line_number}p")

if [[ -z "${ROW}" ]]; then
    echo "ERROR: no circRNA row ${line_number} in ${SAMPLESHEET}" >&2
    exit 1
fi
IFS=$'\t' read -r SM R1 R2 <<< "${ROW}"
echo "sample ${SM}"
echo "  R1: ${R1}"
echo "  R2: ${R2}"

###########
# Trim with fastp
###########
fastp \
    -i "${R1}" -I "${R2}" \
    -o "${TRIMDIR}/${SM}_R1.fq.gz" -O "${TRIMDIR}/${SM}_R2.fq.gz" \
    --detect_adapter_for_pe \
    -q 20 -u 30 \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
    -l 50 \
    -w "${THREADS}" \
    -j "${QCDIR}/${SM}.fastp.json" -h "${QCDIR}/${SM}.fastp.html"

echo "${SM}: done"
