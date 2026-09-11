#!/bin/bash
#SBATCH -J rrna_check_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/rrna_check_totalRNAseq_%j.out  # Output file with the job ID
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/rrna_check_totalRNAseq_%j.err  # Error file with the job ID
#SBATCH -t 10:00:00     # whole batch in one job (D-HH:MM:SS); ~65 samples, subsampled reads
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> sortmerna -a 8
#SBATCH --mem=24GB


## Submit:  sbatch check_rrna_totalRNAseq.sbatch
## Re-do every sample from scratch:  FORCE=1 sbatch check_rrna_totalRNAseq.sbatch
## Point at a different trimmed-FASTQ folder:
##   TRIMDIR=/scratch/alpine/$USER/.../total_RNAseq/trimmed sbatch check_rrna_totalRNAseq.sbatch

# ---------------------------------------------------------------------------
# Ribosomal RNA content screen for the ribo-depleted total RNA-seq ("rnaseq"
# modality) libraries -- diagnostic for high FastQC duplication rates.
#
# Rationale: incomplete rRNA depletion is a common, purely bioinformatically-
# checkable driver of inflated FastQC "Sequence Duplication Levels" in total
# RNA-seq that has nothing to do with PCR duplication or UMIs -- a handful of
# rRNA fragments sequenced millions of times will blow up the duplication
# metric on their own. This is independent of, and should be checked before,
# any conclusion about PCR-driven duplication (see run_dupradar_totalRNAseq.sbatch
# for that separate question).
#
# Method: SortMeRNA (field standard for this) against its own default
# SILVA + Rfam rRNA database, run on a seqtk subsample of each sample's
# trimmed reads (same subsampling convention as run_infer_strandedness_test.sbatch).
#
# Requires a NEW conda env not yet in this repo's README (SortMeRNA is not in
# 'circrna' or 'rseqc'):
#   conda create -n sortmerna -c bioconda -c conda-forge sortmerna=4.3.7
# Pinned to 4.3.7 because the database file downloaded below (smr_v4.3_*) is
# versioned to match the SortMeRNA 4.3.x index format. The preflight check
# below fails fast with this same instruction if the env is missing.
#
# Output (resumable, one row per sample, skips samples already done):
#   ${RRNA_DIR}/<sample>.rrna_sortmerna.log   SortMeRNA's own summary log (kept; per-sample workdir is deleted)
#   ${RRNA_DIR}/rrna_summary.tsv              sample  n_subsampled_pairs  pct_reads_rrna
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
conda activate circrna     # seqtk only, for subsampling

echo ____________________________________
which seqtk
echo ____________________________________
echo ">> Preflight: checking 'sortmerna' conda env"
if ! conda run -n sortmerna sortmerna --version >/dev/null 2>&1; then
    echo "ERROR: conda env 'sortmerna' was not found." >&2
    echo "       Create it with:" >&2
    echo "         conda create -n sortmerna -c bioconda -c conda-forge sortmerna=4.3.7" >&2
    exit 1
fi
conda run -n sortmerna sortmerna --version
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
MODALITY="rnaseq"

REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
RRNA_DB_DIR="${REF_DIR}/rrna_sortmerna"
RRNA_DB="${RRNA_DB_DIR}/smr_v4.3_default_db.fasta"
IDX_DIR="${RRNA_DB_DIR}/idx"          # built once on first sample, reused after
DB_URL="https://github.com/sortmerna/sortmerna/releases/download/v4.3.7/smr_v4.3_default_db.fasta.gz"

ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
TRIMDIR="${TRIMDIR:-${ANALYSIS_DIR}/total_RNAseq/trimmed}"
RRNA_DIR="${ANALYSIS_DIR}/total_RNAseq/QC/rrna_screen"
TMPDIR_LOCAL="${RRNA_DIR}/tmp"

N_READS=2000000     # subsampled read pairs per sample, same as the strandedness check
THREADS="${SLURM_CPUS_PER_TASK:-8}"
FORCE="${FORCE:-0}" # FORCE=1 -> ignore existing per-sample reports and redo

mkdir -p "${RRNA_DB_DIR}" "${IDX_DIR}" "${RRNA_DIR}" "${TMPDIR_LOCAL}"

###########
# 1. Fetch SortMeRNA's default rRNA database (once, reused every run)
###########
# smr_v4.3_default_db.fasta = SortMeRNA's bundled SILVA (arc/bac/euk 16S+18S,
# 23S+28S) + Rfam (5S, 5.8S) rRNA reference set -- the field-standard database
# for this check, not a hand-picked pair of accessions.
if [[ ! -s "${RRNA_DB}" ]]; then
    echo ">> Downloading SortMeRNA default rRNA database"
    curl -L -o "${RRNA_DB}.gz" "${DB_URL}"
    gunzip "${RRNA_DB}.gz"
else
    echo ">> rRNA database already present, skipping: ${RRNA_DB}"
fi

###########
# 2. Sample list (modality == rnaseq), same columns as the other totalRNAseq scripts
###########
#   5 modality   6 sample_id   7 fastq_1   8 fastq_2
mapfile -t SAMPLES < <(
    awk -F',' -v m="${MODALITY}" \
        'NR>1 && $5==m && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
        "${SAMPLESHEET}"
)
N="${#SAMPLES[@]}"
echo "Found ${N} samples for modality='${MODALITY}'"
[[ "${N}" -gt 0 ]] || { echo "ERROR: no rnaseq rows in ${SAMPLESHEET}" >&2; exit 1; }

###########
# Helper: does an existing report look finished?
###########
report_is_complete() {
    local f="$1"
    [[ -s "${f}" ]] || return 1
    grep -q 'Total reads passing E-value threshold' "${f}" || return 1
    return 0
}

###########
# 3. Loop: subsample -> SortMeRNA vs. default rRNA DB -> parse summary log
###########
n_done=0; n_skip=0; n_fail=0
for SM in "${SAMPLES[@]}"; do
    R1="${TRIMDIR}/${SM}_R1.fq.gz"
    R2="${TRIMDIR}/${SM}_R2.fq.gz"
    OUT="${RRNA_DIR}/${SM}.rrna_sortmerna.log"

    if [[ "${FORCE}" != "1" ]] && report_is_complete "${OUT}"; then
        echo "SKIP  ${SM} (already done)"
        n_skip=$((n_skip+1))
        continue
    fi

    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "MISSING trimmed FASTQ for ${SM} (${R1} / ${R2}) -- skipping" >&2
        n_fail=$((n_fail+1))
        continue
    fi

    echo "RUN   ${SM}"
    SUB_R1="${TMPDIR_LOCAL}/${SM}.sub_R1.fq.gz"
    SUB_R2="${TMPDIR_LOCAL}/${SM}.sub_R2.fq.gz"
    seqtk sample -s100 "${R1}" "${N_READS}" | gzip > "${SUB_R1}"
    seqtk sample -s100 "${R2}" "${N_READS}" | gzip > "${SUB_R2}"

    WORKDIR="${TMPDIR_LOCAL}/${SM}.smr_work"
    rm -rf "${WORKDIR}"
    mkdir -p "${WORKDIR}"

    # --paired_in: if either mate hits the rRNA DB, count the pair as rRNA
    # (conservative for a contamination screen). No --fastx/--aligned/--other
    # -- we only need the summary log, not the actual rRNA/non-rRNA reads.
    ok=1
    conda run --no-capture-output -n sortmerna sortmerna \
        --ref "${RRNA_DB}" \
        --reads "${SUB_R1}" --reads "${SUB_R2}" \
        --workdir "${WORKDIR}" \
        --idx-dir "${IDX_DIR}" \
        --paired_in \
        -a "${THREADS}" \
        -v || ok=0

    LOG="$(find "${WORKDIR}" -iname '*.log' -print -quit 2>/dev/null || true)"
    if [[ "${ok}" -eq 1 && -n "${LOG}" ]] && report_is_complete "${LOG}"; then
        cp -f "${LOG}" "${OUT}"
        n_done=$((n_done+1))
        grep 'Total reads passing E-value threshold' "${OUT}" | sed "s/^/      ${SM}: /"
    else
        echo "FAIL  ${SM}" >&2
        n_fail=$((n_fail+1))
    fi

    rm -f "${SUB_R1}" "${SUB_R2}"
    rm -rf "${WORKDIR}"   # keep only the copied-out log; kvdb/readb caches are per-sample and large
done

echo "----------------------------------------------------------------------"
echo "rRNA screen: ${n_done} generated, ${n_skip} skipped, ${n_fail} failed"

###########
# 4. Build the cohort summary table
###########
SUMMARY="${RRNA_DIR}/rrna_summary.tsv"
{
    printf 'sample\tn_subsampled_pairs\tpct_reads_rrna\n'
    for r in "${RRNA_DIR}"/*.rrna_sortmerna.log; do
        [[ -e "${r}" ]] || continue
        sm="$(basename "${r}")"; sm="${sm%.rrna_sortmerna.log}"
        pct=$(grep 'Total reads passing E-value threshold' "${r}" | sed -E 's/.*\(([0-9.]+)\).*/\1/')
        printf '%s\t%s\t%s\n' "${sm}" "${N_READS}" "${pct:-NA}"
    done
} > "${SUMMARY}"

echo "wrote ${SUMMARY}"
column -t -s $'\t' "${SUMMARY}" || cat "${SUMMARY}"

cat <<'EOF'

============================================================
 HOW TO READ THIS
============================================================
 pct_reads_rrna = % of a 2M-read-pair subsample that SortMeRNA assigns to its
 default rRNA database (SILVA 16S/18S/23S/28S + Rfam 5S/5.8S) at the default
 E-value threshold:

   roughly <2-3%  ..... ribo-depletion looks effective; rRNA is unlikely to
                        be a major driver of the FastQC duplication numbers.
   roughly 5-20%  ..... some residual rRNA; worth noting as a contributor,
                        probably not the whole story.
   >20-30%+ ........... incomplete ribo-depletion is likely a real, separate
                        contributor to the high duplication rates -- this has
                        nothing to do with PCR duplication or UMIs and is not
                        fixable retroactively in the data you already have.

 Also worth a quick look regardless of these numbers: FastQC's
 "Overrepresented sequences" table for a couple of samples -- if the top
 entries look rRNA-like, that corroborates this screen.
============================================================
EOF

echo "done"
