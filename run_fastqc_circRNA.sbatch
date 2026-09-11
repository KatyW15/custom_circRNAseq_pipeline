#!/bin/bash
#SBATCH -J fastqc_circRNA
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastqc_circRNA_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/fastqc_circRNA_%j.err
#SBATCH -t 1-00:00:00   # D-HH:MM:SS
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> FastQC -t 8
#SBATCH --mem=16GB

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge

#activate my circrna env with all dependancies
conda activate circrna


echo ____________________________________
which fastqc
fastqc --version
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
OUTDIR="/projects/$USER/HNSCC_analysis/my_circRNAseq_analysis/QC/circRNAseq/raw_FastQC"
MODALITY="circrnaseq"          # value in the 'modality' column to keep
THREADS="${SLURM_CPUS_PER_TASK:-8}"

FILELIST="/projects/$USER/HNSCC_analysis/my_circRNAseq_analysis/circRNA_fastq_files.txt"

mkdir -p "${OUTDIR}"

###########
# Build the list of FASTQ files for circRNA samples only
###########
# Columns in master_samplesheet.csv:
#   1 patient_id  2 batch_id  3 tissue  4 hpv_status  5 modality
#   6 sample_id   7 fastq_1   8 fastq_2  9 survival_time_days
#  10 vital_status 11 notes
# fastq_1 / fastq_2 are columns 7 and 8 and contain no commas, so a plain
# comma split is safe even though 'notes' (col 11) may contain commas.
# Some circRNA rows have placeholder text ("no/no corresponding fastq") in the
# fastq columns instead of a real path — the regex keeps only real .fastq.gz /
# .fq.gz paths and drops those.
awk -F',' 'NR > 1 && $5 == "'"${MODALITY}"'" {
    if ($7 ~ /\.(fastq|fq)\.gz$/) print $7
    if ($8 ~ /\.(fastq|fq)\.gz$/) print $8
}' "${SAMPLESHEET}" | sort -u > "${FILELIST}"

N=$(wc -l < "${FILELIST}")
echo "Found ${N} FASTQ files for modality='${MODALITY}'"
if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: no files matched — check SAMPLESHEET path and MODALITY value." >&2
    exit 1
fi

# Fail early if any listed file is missing
MISSING=0
while IFS= read -r f; do
    if [[ ! -s "${f}" ]]; then
        echo "MISSING: ${f}" >&2
        MISSING=1
    fi
done < "${FILELIST}"
[[ "${MISSING}" -eq 0 ]] || { echo "ERROR: some FASTQ files are missing (see above)." >&2; exit 1; }

###########
# Run FastQC
###########
echo "Running FastQC with ${THREADS} threads -> ${OUTDIR}"
mapfile -t FILES < "${FILELIST}"
fastqc --threads "${THREADS}" --outdir "${OUTDIR}" "${FILES[@]}"

###########
# Aggregate with MultiQC (optional; skipped if module not present)
###########
if command -v multiqc >/dev/null 2>&1; then
    echo "Running MultiQC"
    multiqc --force --outdir "${OUTDIR}" --filename multiqc_circRNA_fastqc "${OUTDIR}"
fi

echo "Done. Reports in ${OUTDIR}"
