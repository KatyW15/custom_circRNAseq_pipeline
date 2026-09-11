#!/bin/bash
#SBATCH -J ciriquant_circRNA
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_circRNA_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_circRNA_%A_%a.err
#SBATCH -t 16:00:00    # per-sample wall time (D-HH:MM:SS); BWA + HISAT2 + StringTie
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> CIRIquant -t 8
#SBATCH --mem=48GB
#SBATCH --array=1-57

## hard coded array above, if want to define array in command line see below
## Submit one array task per circRNA sample. Get the count (= array range) with:
#   awk -F',' 'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l
#   sbatch --array=1-57 run_ciriquant_circRNA.sbatch

# ---------------------------------------------------------------------------
# STEP 5 (CIRIquant) of the circRNA-seq pipeline -- see
# circRNA_analysis_pipeline_outline.txt
#
# Quantifies a FIXED consensus locus set in every sample (prior-guided mode,
# --bed) so the downstream matrix has consistent rows. CIRIquant internally:
# HISAT2 + StringTie (linear model) and BWA + pseudo-circular realignment.
#
# Input : trimmed FASTQs from run_fastp_circRNA.sbatch
#           ${TRIMDIR}/${SM}_R1.fq.gz , ${SM}_R2.fq.gz
#         consensus BED from build_bsj_consensus.sbatch
#           ${CONSENSUS_DIR}/cohort.consensus.union.bed   (default; whole cohort)
#           -- or per-sample ${CONSENSUS_DIR}/${SM}.consensus.bed (set CONSENSUS_BED)
# Output: ${OUTDIR}/${SM}/${SM}.gtf   per circRNA: BSJ, FSJ, junction ratio, CPM
#         ${OUTDIR}/${SM}/${SM}.log
#         ${OUTDIR}/${SM}/{align,circ,gene}/  intermediates
#
# RNase R caveat (outline STEP 6): these libraries are RNase R treated, so the
# junction ratio / CPM here are detection-scaled, NOT corrected for enrichment.
# Accurate cross-sample expression needs a matched RNase R(-) mock per sample
# + `CIRIquant --RNaseR` -- a later step. Report BSJ counts until then.
#
# --library-type 2  == reverse / fr-firststrand, confirmed by RSeQC in step 0.
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
#activate my circrna env with all dependancies (bwa, hisat2, stringtie, samtools, ciriquant)
conda activate circrna

echo ____________________________________
which CIRIquant bwa hisat2 stringtie samtools
CIRIquant --version 2>&1 | head -1 || true
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
TRIMDIR="${ANALYSIS_DIR}/circRNAseq/trimmed"   # trimmed FASTQs from run_fastp_circRNA.sbatch
CONSENSUS_DIR="${ANALYSIS_DIR}/consensus"  # from build_bsj_consensus.sbatch (cohort-wide, not modality-nested)
OUTDIR="${ANALYSIS_DIR}/ciriquant"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
LIBRARY_TYPE=2                             # reverse (fr-firststrand)

# Same references as the rest of the HNSCC build. bwa_index prefix == the FASTA
# path (bwa index was run on the .fa); hisat_index prefix has no extension.
REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
FASTA="${REF_DIR}/hg38/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode/gencode.v50.annotation.gtf"
BWA_INDEX="${FASTA}"
HISAT_INDEX="${REF_DIR}/hg38/GRCh38.primary_assembly.genome_hisat2"

# Locus set to quantify. Default: the cohort-wide union (consistent matrix rows
# across all samples). For per-sample loci instead:
#   sbatch --export=ALL,CONSENSUS_BED_MODE=persample run_ciriquant_circRNA.sbatch
CONSENSUS_BED_MODE="${CONSENSUS_BED_MODE:-union}"

mkdir -p "${OUTDIR}"

###########
# Pull this task's sample from the master samplesheet
###########
# Identical filter to run_fastp/bwa/ciri2, so array index N is the same sample.
#   6 sample_id   7 fastq_1   8 fastq_2
line_number=$SLURM_ARRAY_TASK_ID
echo "array number ${SLURM_ARRAY_TASK_ID}"

SM=$(awk -F',' \
    'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
    "${SAMPLESHEET}" | sed -n "${line_number}p")

if [[ -z "${SM}" ]]; then
    echo "ERROR: no circRNA row ${line_number} in ${SAMPLESHEET}" >&2
    exit 1
fi
echo "sample ${SM}"

R1="${TRIMDIR}/${SM}_R1.fq.gz"
R2="${TRIMDIR}/${SM}_R2.fq.gz"

if [[ "${CONSENSUS_BED_MODE}" == "persample" ]]; then
    CONSENSUS_BED="${CONSENSUS_DIR}/${SM}.consensus.bed"
else
    CONSENSUS_BED="${CONSENSUS_DIR}/cohort.consensus.union.bed"
fi

###########
# Preflight
###########
[[ -s "${R1}" ]]            || { echo "ERROR: missing trimmed R1: ${R1}" >&2; exit 1; }
[[ -s "${R2}" ]]            || { echo "ERROR: missing trimmed R2: ${R2}" >&2; exit 1; }
[[ -s "${CONSENSUS_BED}" ]] || { echo "ERROR: missing consensus BED: ${CONSENSUS_BED} (run build_bsj_consensus.sbatch)" >&2; exit 1; }
[[ -s "${FASTA}" ]]         || { echo "ERROR: missing reference FASTA: ${FASTA}" >&2; exit 1; }
[[ -s "${FASTA}.fai" ]]    || { echo "ERROR: missing ${FASTA}.fai (samtools faidx)" >&2; exit 1; }
[[ -s "${GTF}" ]]          || { echo "ERROR: missing annotation GTF: ${GTF}" >&2; exit 1; }
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
echo "bed   : ${CONSENSUS_BED} ($(wc -l < "${CONSENSUS_BED}") loci, mode=${CONSENSUS_BED_MODE})"

###########
# Quantify
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
