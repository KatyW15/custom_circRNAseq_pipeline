#!/bin/bash
#SBATCH -J ciriquant_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_totalRNAseq_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_totalRNAseq_%A_%a.err
#SBATCH -t 16:00:00    # per-sample wall time (D-HH:MM:SS); BWA + HISAT2 + StringTie
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> CIRIquant -t 8
#SBATCH --mem=48GB
#SBATCH --array=1-65

#
# NOT part of job_inchain_totalRNAseq.sh -- run this separately, and only after
# run_ciriquant_circRNA.sbatch has finished for the matched samples.

## 65 = total RNA-seq samples (modality "rnaseq" with real paired FASTQ paths).
## ~8 have no circRNA-seq counterpart -> those tasks exit 0 (skip). Recount with:
#   awk -F',' 'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l

# ---------------------------------------------------------------------------
# STEP 6a (RNase R effect correction) -- see circRNA_analysis_pipeline_outline.txt
#
# Runs CIRIquant on the matched RNase R(-) total RNA-seq library with --RNaseR
# pointing at that sample's RNase R(+) CIRIquant result from
# run_ciriquant_circRNA.sbatch. CIRIquant detects circRNAs in this mock library,
# takes the overlap with the RNase R(+) circRNA set, fits a Gaussian mixture
# model to the per-circRNA enrichment ratios, and writes enrichment-corrected
# (pre-treatment) expression. The header of the output gtf carries the estimated
# RNase R efficiency.
#
# KEEPS BOTH:
#   uncorrected  = ${ANALYSIS_DIR}/ciriquant/<circRNA sample>/...       (untouched)
#   corrected    = ${OUTDIR}/<rnaseq sample>/<rnaseq sample>.gtf        (this script)
#
# Pairing (sample_id = <patient>_<tissue>_<modality>):
#   TP#####_tumor_rnaseq   <->   TP#####_tumor_circrnaseq
# rnaseq samples with no circRNA-seq counterpart cannot be corrected -> skip (exit 0).
#
# NOT passed:
#   --bed      : the circRNA locus set is inherited from the --RNaseR gtf;
#                CIRIquant still detects de novo in the mock for the GMM overlap.
#   --no-gene  : the mock has real linear signal -- keep gene quantification on
#                (CIRIquant uses it for normalisation; gene table is a byproduct).
#
# !!! LIBRARY_TYPE below is NOT confirmed for the rnaseq libraries. They are a
# different prep/batch from the circRNA-seq set -- do NOT assume they are also
# reverse. Confirm first: point run_infer_strandedness_test.sbatch at one
# trimmed rnaseq sample (edit SAMPLE + the _circrnaseq_ infix) and set the value.
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
TRIMDIR="${ANALYSIS_DIR}/total_RNAseq/trimmed"          # trimmed FASTQs from run_fastp_totalRNAseq.sbatch
CIRC_CIRIQUANT_DIR="${ANALYSIS_DIR}/circRNAseq/ciriquant_noRNaseRcorrection"   # RNase R(+) results from run_ciriquant_circRNA.sbatch
OUTDIR="${ANALYSIS_DIR}/ciriquant/total_RNAseq"         # enrichment-corrected results
THREADS="${SLURM_CPUS_PER_TASK:-8}"

LIBRARY_TYPE="${LIBRARY_TYPE:-2}"             # !!! CONFIRM for rnaseq (see header). 0 unstr / 1 fwd / 2 rev
CIRC_MODALITY_SUFFIX="_circrnaseq"            # how the circRNA-seq sample_id ends
RNA_MODALITY_SUFFIX="_rnaseq"                 # how the rnaseq sample_id ends

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

# Matched RNase R(+) circRNA-seq CIRIquant result
CIRC_SM="${SM%${RNA_MODALITY_SUFFIX}}${CIRC_MODALITY_SUFFIX}"
RNASER_GTF="${CIRC_CIRIQUANT_DIR}/${CIRC_SM}/${CIRC_SM}.gtf"
echo "matched circRNA-seq sample: ${CIRC_SM}"
echo "RNase R(+) gtf            : ${RNASER_GTF}"

if [[ ! -s "${RNASER_GTF}" ]]; then
    echo "${SM}: no matched circRNA-seq CIRIquant result -- cannot RNase R-correct this sample. Skipping." >&2
    exit 0
fi

###########
# Preflight
###########
[[ -s "${R1}" ]]                || { echo "ERROR: missing trimmed R1: ${R1} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }
[[ -s "${R2}" ]]                || { echo "ERROR: missing trimmed R2: ${R2} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }
[[ -s "${FASTA}" ]]            || { echo "ERROR: missing reference FASTA: ${FASTA}" >&2; exit 1; }
[[ -s "${FASTA}.fai" ]]       || { echo "ERROR: missing ${FASTA}.fai (samtools faidx)" >&2; exit 1; }
[[ -s "${GTF}" ]]             || { echo "ERROR: missing annotation GTF: ${GTF}" >&2; exit 1; }
[[ -s "${BWA_INDEX}.bwt" ]]     || { echo "ERROR: BWA index not found for ${BWA_INDEX} (run run_bwa_index.sbatch)" >&2; exit 1; }
[[ -s "${HISAT_INDEX}.1.ht2" ]] || { echo "ERROR: HISAT2 index not found at ${HISAT_INDEX} (built in step 0)" >&2; exit 1; }

###########
# Working CIRIquant config (absolute tool paths from this env); one per task
###########
# Same template as run_ciriquant_circRNA.sbatch / ciriquant_config.yaml.
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
echo "library-type: ${LIBRARY_TYPE}"

###########
# Quantify + RNase R correction
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
        -o "${OUTDIR}/${SM}" \
        -p "${SM}" \
        --RNaseR "${RNASER_GTF}"
fi

###########
# Report the estimated RNase R efficiency from the output header
###########
if [[ -s "${GTF_OUT}" ]]; then
    echo "${SM}: RNase R correction header ->"
    grep -i -m5 -E 'RNaseR|efficiency|enrich' "${GTF_OUT}" || head -n 15 "${GTF_OUT}" | grep '^#' || true
    N_CIRC=$(awk '$3=="circRNA"' "${GTF_OUT}" | wc -l || echo 0)
    echo "${SM}: ${N_CIRC} circRNAs in corrected output -> ${GTF_OUT}"
fi

echo "${SM}: done"
