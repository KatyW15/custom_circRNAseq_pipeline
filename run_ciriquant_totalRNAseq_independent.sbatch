#!/bin/bash
#SBATCH -J ciriquant_totalRNAseq_independent
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_totalRNAseq_independent_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_totalRNAseq_independent_%A_%a.err
#SBATCH -t 16:00:00    # per-sample wall time (D-HH:MM:SS); BWA + HISAT2 + StringTie
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> CIRIquant -t 8
#SBATCH --mem=48GB
#SBATCH --array=1-65

## 65 = total RNA-seq samples (modality "rnaseq" with real paired FASTQ paths).
## No circRNA-seq pairing needed for this run, so all 65 get quantified. Recount
## if the samplesheet changes:
#   awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l

# ---------------------------------------------------------------------------
# Independent CIRIquant quantification of the total RNA-seq libraries.
#
# Plain CIRIquant on the total RNA-seq reads -- no --RNaseR. The RNase R(+)/
# RNase R(-) pairing/correction is still being decided (see
# run_ciriquant_totalRNAseq.sbatch); this script does not depend on it and does
# not need a matched circRNA-seq sample. It quantifies the same consensus locus
# set as the circRNA-seq run so the three planned analyses share the same rows:
#   1. circRNA-seq, uncorrected     run_ciriquant_circRNA.sbatch
#   2. RNase R corrected (TBD)      run_ciriquant_totalRNAseq.sbatch
#   3. total RNA-seq, independent   this script
#
# Unlike the RNase R(+) circRNA-seq libraries, junction ratio here IS meaningful
# -- total RNA-seq is not RNase R-digested, so FSJ reads are real.
#
# Input : trimmed FASTQs from run_fastp_totalRNAseq.sbatch
#           ${TRIMDIR}/${SM}_R1.fq.gz , ${SM}_R2.fq.gz
#         consensus BED from build_bsj_consensus.sbatch (circRNA-seq branch)
#           ${CONSENSUS_DIR}/cohort.consensus.union.bed
# Output: ${OUTDIR}/${SM}/${SM}.gtf   per circRNA: BSJ, FSJ, junction ratio, CPM
#
# !!! LIBRARY_TYPE below is NOT confirmed for the rnaseq libraries -- different
# prep/batch from the circRNA-seq set. Confirm with run_infer_strandedness_test.sbatch
# pointed at a trimmed rnaseq sample before trusting the junction-ratio/CPM values.
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
#activate my ciriquant env (bwa=0.7.17, hisat2=2.2.0, stringtie=2.1.1, samtools>=1.10, ciriquant, python=3.9)
#separate from circrna env because CIRIquant breaks on python>=3.12 (distutils removed)
conda activate ciriquant

echo ____________________________________
which CIRIquant bwa hisat2 stringtie samtools
CIRIquant --version 2>&1 | head -1 || true
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
MODALITY_DIR="${ANALYSIS_DIR}/total_RNAseq"           # trimmed/bwa_mapped/ciri2 all nest here
TRIMDIR="${MODALITY_DIR}/trimmed"                     # trimmed FASTQs from run_fastp_totalRNAseq.sbatch
CONSENSUS_DIR="${ANALYSIS_DIR}/circRNAseq/consensus"  # from build_bsj_consensus.sbatch (circRNA-seq branch)
OUTDIR="${MODALITY_DIR}/ciriquant_noRNaseRcorrection" # independent quantification (no --RNaseR)
THREADS="${SLURM_CPUS_PER_TASK:-8}"

LIBRARY_TYPE="${LIBRARY_TYPE:-2}"   # !!! CONFIRM for rnaseq (see header). 0 unstr / 1 fwd / 2 rev

# Same references as the rest of the HNSCC build.
REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
FASTA="${REF_DIR}/hg38/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode/gencode.v50.annotation.gtf"
BWA_INDEX="${FASTA}"
HISAT_INDEX="${REF_DIR}/hg38/GRCh38.primary_assembly.genome_hisat2"

mkdir -p "${OUTDIR}"

###########
# Pull this task's sample from the master samplesheet
###########
# Identical filter to run_fastp/bwa/ciri2_totalRNAseq.sbatch.
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
CONSENSUS_BED="${CONSENSUS_DIR}/cohort.consensus.union.bed"

###########
# Preflight
###########
[[ -s "${R1}" ]]                || { echo "ERROR: missing trimmed R1: ${R1} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }
[[ -s "${R2}" ]]                || { echo "ERROR: missing trimmed R2: ${R2} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }
[[ -s "${CONSENSUS_BED}" ]]     || { echo "ERROR: missing consensus BED: ${CONSENSUS_BED} (run build_bsj_consensus.sbatch)" >&2; exit 1; }
[[ -s "${FASTA}" ]]             || { echo "ERROR: missing reference FASTA: ${FASTA}" >&2; exit 1; }
[[ -s "${FASTA}.fai" ]]         || { echo "ERROR: missing ${FASTA}.fai (samtools faidx)" >&2; exit 1; }
[[ -s "${GTF}" ]]               || { echo "ERROR: missing annotation GTF: ${GTF}" >&2; exit 1; }
[[ -s "${BWA_INDEX}.bwt" ]]     || { echo "ERROR: BWA index not found for ${BWA_INDEX} (run run_bwa_index.sbatch)" >&2; exit 1; }
[[ -s "${HISAT_INDEX}.1.ht2" ]] || { echo "ERROR: HISAT2 index not found at ${HISAT_INDEX} (built in step 0)" >&2; exit 1; }

###########
# Working CIRIquant config (absolute tool paths from this env); one per task
###########
CONFIG="$(mktemp "${OUTDIR}/config.${SM}.XXXXXX.yaml")"
trap 'rm -f "${CONFIG}"' EXIT
cat > "${CONFIG}" <<YAML
name: hg38_gencode_v50
tools:
  bwa: $(command -v bwa)
  hisat2: $(command -v hisat2)
  stringtie: $(command -v stringtie)
  samtools: $(command -v samtools)
reference:
  fasta: ${FASTA}
  gtf: ${GTF}
  bwa_index: ${BWA_INDEX}
  hisat_index: ${HISAT_INDEX}
YAML
echo "config:"; sed 's/^/    /' "${CONFIG}"
echo "bed   : ${CONSENSUS_BED} ($(wc -l < "${CONSENSUS_BED}") loci)"
echo "library-type: ${LIBRARY_TYPE}"

###########
# Quantify (no correction)
###########
GTF_OUT="${OUTDIR}/${SM}/${SM}.gtf"
if [[ -s "${GTF_OUT}" ]]; then
    echo "${SM}: ${GTF_OUT} already exists, skipping CIRIquant"
else
    CIRIquant \
        -t "${THREADS}" \
        -1 "${R1}" -2 "${R2}" \
        --config "${CONFIG}" \
        --library-type "${LIBRARY_TYPE}" \
        --bed "${CONSENSUS_BED}" \
        -o "${OUTDIR}/${SM}" \
        -p "${SM}"
fi

N_CIRC=$(awk '$3=="circRNA"' "${GTF_OUT}" 2>/dev/null | wc -l || echo 0)
echo "${SM}: ${N_CIRC} circRNAs quantified -> ${GTF_OUT}"
echo "${SM}: done"
