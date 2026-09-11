#!/bin/bash
#SBATCH -J ciri2_totalRNAseq
#SBATCH --partition amem
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri2_totalRNAseq_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri2_totalRNAseq_%A_%a.err
#SBATCH -t 03:00:00   # per-sample wall time (D-HH:MM:SS); CIRI2 is slow on deep libraries
#SBATCH --qos=mem-normal
#SBATCH -n 1 -c 8       # 8 cores -> CIRI2 -T 8
#SBATCH --mem=256GB     # CIRI2 holds all candidate junction reads in RAM during
                        # scan 1; needs at least 100GB, running with minimum amem mem request
                        # (256GB) because amem node wait times are less.
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=YOUR_EMAIL@example.com
#SBATCH --array=1-65


## 65 = total RNA-seq samples (modality "rnaseq" with real paired FASTQ paths).
## Recount if the samplesheet changes:
#   awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l

# ---------------------------------------------------------------------------
# VALIDATION branch, detect step -- see circRNA_analysis_pipeline_outline.txt STEP 6.
#
# CIRI2 on the RNase R(-) total RNA-seq SAM. Identical command/filters to
# run_ciri2_circRNA.sbatch so the BSJ call sets are directly comparable.
# Output goes to ciri2/total_RNAseq/ (separate from the circRNA-seq ciri2/ dir);
# build_bsj_validation.sbatch intersects the two.
#
# Input : name-ordered BWA-MEM SAM from run_bwa_totalRNAseq.sbatch
#           ${ALIGNDIR}/${SM}.sam   (NOT coordinate-sorted)
# Output: ${OUTDIR}/${SM}.ciri          raw CIRI2 BSJ table (1-based coords)
#         ${OUTDIR}/${SM}.CIRI2.log     CIRI2 run log
#         ${OUTDIR}/${SM}.bsj2.ciri     first-pass filter: junction_reads >= 2
#         ${OUTDIR}/${SM}.bsj2.bed      BED6 of the >=2-read BSJ loci (0-based
#                                       half-open) -> build_bsj_validation.sbatch
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
#activate my circrna env with all dependancies
conda activate circrna

# CIRI2 ships as a standalone Perl script. Same resolution logic as
# run_ciri2_circRNA.sbatch (bioconda 'ciri2' puts CIRI2.pl on PATH; otherwise
# fall back to the CIRIquant-bundled copy). Override with CIRI2_PL=/path.
CIRI2_PL="${CIRI2_PL:-$(command -v CIRI2.pl || true)}"
if [[ -z "${CIRI2_PL}" && -n "${CONDA_PREFIX:-}" ]]; then
    CIRI2_PL="$(find "${CONDA_PREFIX}" -name 'CIRI2.pl' -type f 2>/dev/null | head -n1)"
fi
[[ -n "${CIRI2_PL}" ]] || { echo "ERROR: CIRI2.pl not found on PATH or under ${CONDA_PREFIX:-<no conda env>} (set CIRI2_PL=/path/to/CIRI2.pl)" >&2; exit 1; }

echo ____________________________________
which bwa perl
perl "${CIRI2_PL}" --help 2>&1 | grep -i -m1 version || echo "CIRI2.pl: ${CIRI2_PL}"
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
ALIGNDIR="${ANALYSIS_DIR}/bwa_mapped/total_RNAseq"   # input SAMs from run_bwa_totalRNAseq.sbatch
OUTDIR="${ANALYSIS_DIR}/ciri2/total_RNAseq"          # separate from the circRNA-seq ciri2/ dir
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Same references as the rest of the HNSCC build
REF="/projects/$USER/HNSCC_analysis/reference_genomes/hg38/GRCh38.primary_assembly.genome.fa"
GTF="/projects/$USER/HNSCC_analysis/reference_genomes/gencode/gencode.v50.annotation.gtf"

# Set to 1 to remove the (large) BWA SAM once CIRI2 succeeds.
DELETE_SAM="${DELETE_SAM:-0}"

mkdir -p "${OUTDIR}"

###########
# Pull this task's sample from the master samplesheet
###########
# Identical filter to run_fastp/bwa_totalRNAseq.sbatch, so array index N maps to
# the same sample across every total RNA step.
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

SAM="${ALIGNDIR}/${SM}.sam"
OUT="${OUTDIR}/${SM}.ciri"
LOG="${OUTDIR}/${SM}.CIRI2.log"

[[ -s "${SAM}" ]] || { echo "ERROR: missing BWA SAM: ${SAM} (run run_bwa_totalRNAseq.sbatch)" >&2; exit 1; }
[[ -s "${REF}" ]] || { echo "ERROR: missing reference FASTA: ${REF}" >&2; exit 1; }
[[ -s "${GTF}" ]] || { echo "ERROR: missing annotation GTF: ${GTF}" >&2; exit 1; }

###########
# Detect back-splice junctions with CIRI2
###########
if [[ -s "${OUT}" ]]; then
    echo "${SM}: ${OUT} already exists, skipping CIRI2"
else
    perl "${CIRI2_PL}" \
        -I "${SAM}" \
        -O "${OUT}" \
        -F "${REF}" \
        -A "${GTF}" \
        -G "${LOG}" \
        -T "${THREADS}"
fi

###########
# First-pass filter (#junction_reads >= 2) + BED of BSJ loci for the validation step
###########
# CIRI2 columns: 1 circRNA_ID  2 chr  3 circRNA_start(1-based)  4 circRNA_end
#                5 #junction_reads ... 11 strand
BSJ2="${OUTDIR}/${SM}.bsj2.ciri"
BED="${OUTDIR}/${SM}.bsj2.bed"

awk 'BEGIN{FS=OFS="\t"} NR==1 || $5+0 >= 2' "${OUT}" > "${BSJ2}"

awk 'BEGIN{FS=OFS="\t"} NR>1 && $5+0 >= 2 {print $2, $3-1, $4, $1, $5, $11}' "${OUT}" \
    | sort -k1,1 -k2,2n > "${BED}"

N_ALL=$(($(wc -l < "${OUT}") - 1))
N_BSJ2=$(wc -l < "${BED}")
echo "${SM}: ${N_ALL} circRNAs called, ${N_BSJ2} with >=2 junction reads"

###########
# Optionally reclaim scratch
###########
if [[ "${DELETE_SAM}" == "1" ]]; then
    echo "${SM}: removing ${SAM}"
    rm -f "${SAM}"
fi

echo "${SM}: done"
