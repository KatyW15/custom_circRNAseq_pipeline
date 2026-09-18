#!/bin/bash
#SBATCH -J ciriquant_circRNA_py2
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_circRNA_py2_%A_%a.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciriquant_circRNA_py2_%A_%a.err
#SBATCH -t 16:00:00    # per-sample wall time (D-HH:MM:SS); BWA + HISAT2 + StringTie
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> CIRIquant -t 8
#SBATCH --mem=48GB
#SBATCH --array=1-57

## hard coded array above, if want to define array in command line see below
## Submit one array task per circRNA sample. Get the count (= array range) with:
#   awk -F',' 'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/' \
#       /projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv | wc -l
#   sbatch --array=1-57 run_ciriquant_circRNA_py2.sbatch
#
## First time using this env: sanity-check on ONE sample before committing all 57 --
#   sbatch --array=1 run_ciriquant_circRNA_py2.sbatch

# ---------------------------------------------------------------------------
# STEP 5 (CIRIquant).
#
# Runs in the `CIRIquant_github_env` env: the official CIRIquant v1.1.3
# release installed per the developers' own documented method, on its
# natively supported Python 2 interpreter -- see
# https://ciri-cookbook.readthedocs.io/en/latest/CIRIquant_1_installation.html
#
# This is the primary/only CIRIquant path for this pipeline. An earlier
# attempt ran bioconda's CIRIquant under a patched Python 3.11 env instead;
# it kept surfacing new upstream Python-2-only bugs one at a time across
# multi-hour runs and was still producing bad output on some samples after
# five separate fixes, while this native Python 2 install ran the full
# cohort cleanly on the first attempt with zero patches needed. That env and
# its patch scripts are retired -- see software_requirements.md "Retired:
# ciriquant" for the history.
#
# Output: ${ANALYSIS_DIR}/ciriquant_noRNaseRcorrection_py2/<sample>/<sample>.gtf
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
# CIRIquant_github_env: official CIRIquant v1.1.3 install (CIRI-cookbook
# Method 2 pinned environment.yml -- python=2.7.15, pysam==0.15.2, etc.),
# moved into the standard envs dir after being built elsewhere.
#
# set +u around activation only: this env's compiler-toolchain packages
# (binutils_linux-64 and friends) ship activate.d hook scripts that reference
# variables like $ADDR2LINE without a default -- fine under normal bash, but
# fatal ("unbound variable") under our set -u. Not a problem with the env;
# just an under-defensive activation script. Restore -u immediately after.
set +u
conda activate CIRIquant_github_env
set -u

echo ____________________________________
which CIRIquant bwa hisat2 stringtie samtools
CIRIquant --version 2>&1 | head -1 || true
python --version 2>&1 || true
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/circRNAseq"
TRIMDIR="${ANALYSIS_DIR}/trimmed"   # trimmed FASTQs from run_fastp_circRNA.sbatch (shared, unchanged)
CONSENSUS_DIR="${ANALYSIS_DIR}/consensus"  # from build_bsj_consensus.sbatch (shared, unchanged)
OUTDIR="${ANALYSIS_DIR}/ciriquant_noRNaseRcorrection_py2"   # SEPARATE from the patched-build OUTDIR
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
# across all samples, and directly comparable to the patched-build run). For
# per-sample loci instead:
#   sbatch --export=ALL,CONSENSUS_BED_MODE=persample run_ciriquant_circRNA_py2.sbatch
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
# check variables
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
name: hg38_gencode_v50_py2
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
GTF_IS_COMPLETE=0
if [[ -s "${GTF_OUT}" ]] && [[ "$(awk '$3=="circRNA"' "${GTF_OUT}" | wc -l)" -gt 0 ]]; then
    GTF_IS_COMPLETE=1
fi

if [[ "${GTF_IS_COMPLETE}" -eq 1 ]]; then
    echo "${SM}: ${GTF_OUT} already exists with circRNA rows, skipping CIRIquant"
else
    if [[ -s "${GTF_OUT}" ]]; then
        echo "${SM}: ${GTF_OUT} exists but has 0 circRNA rows -- stale/incomplete output from an earlier failed run; removing and re-running"
        rm -rf "${OUTDIR}/${SM}"
    fi
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
