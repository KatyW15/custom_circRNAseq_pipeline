#!/bin/bash
# ---------------------------------------------------------------------------
# Verify every sample made it through the total RNA-seq validation chain
# (job_inchain_totalRNAseq.sh: fastp -> bwa -> ciri2 -> fastqc / flagstat)
# without opening each .out/.err file by hand.
#
# Run interactively on a login node (no sbatch needed):
#     bash check_totalRNAseq_chain.sh
#
# 1. sacct summary   -- per stage, count of array tasks by exit state
#                        (last 7 days; widen with SACCT_SINCE=now-30days)
# 2. per-sample table -- for each of the N samples in master_samplesheet.csv,
#                        checks the real output file for each stage, not just
#                        "the job finished":
#     fastp    : trimmed R1/R2 FASTQs + fastp.json report exist, non-empty
#     bwa      : SAM exists + non-empty, per-sample flagstat report present
#     ciri2    : .ciri + .bsj2.bed exist, CIRI2 log has no error markers
#     fastqc   : raw_FastQC has a *_fastqc.zip for both mates      (batch job)
#     flagstat : sample has a row in mapping_summary.tsv           (batch job)
# 3. error grep      -- slurm_out/*_totalRNAseq_*.err for real failure markers
#                        (OOM, time-limit kill, traceback, segfault), skipping
#                        the normal `which`/version chatter these scripts print
#
# Exits non-zero (and lists the offending samples) if anything is missing.
# ---------------------------------------------------------------------------

set -uo pipefail

SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
MODALITY_DIR="${ANALYSIS_DIR}/total_RNAseq"    # trimmed/bwa_mapped/ciri2 all nest here now
QC_DIR="/projects/$USER/HNSCC_analysis/my_circRNAseq_analysis/QC/totalRNAseq"
SLURM_OUT="${ANALYSIS_DIR}/slurm_out"          # not part of the modality reorg

TRIMDIR="${MODALITY_DIR}/trimmed"
FASTP_QCDIR="${QC_DIR}/trimmed"
ALIGNDIR="${MODALITY_DIR}/bwa_mapped"
FLAGSTATDIR="${ALIGNDIR}/flagstat"
CIRI2DIR="${MODALITY_DIR}/ciri2"
FASTQCDIR="${QC_DIR}/raw_FastQC"

# Pre-reorg location, checked as a fallback in case historical fastp reports
# from before the migration were never physically moved to the path above.
# (Lived on /projects, not /scratch, pre-reorg.)
FASTP_QCDIR_LEGACY="/projects/$USER/HNSCC_analysis/my_circRNAseq_analysis/QC/total_RNAseq/trimmed"

SACCT_SINCE="${SACCT_SINCE:-now-7days}"

[[ -f "${SAMPLESHEET}" ]] || { echo "ERROR: samplesheet not found: ${SAMPLESHEET}" >&2; exit 1; }

# Same filter as job_inchain_totalRNAseq.sh / each *_totalRNAseq.sbatch, so
# array index N maps to the same sample here as it did in the chain. Also pull
# the raw fastq_1/fastq_2 paths (cols 7,8): run_fastqc_totalRNAseq.sbatch runs
# FastQC directly on those raw paths (never renamed to the sample_id), so its
# output zip names come from the RAW file's own basename, not ${SM}.
mapfile -t SAMPLE_ROWS < <(awk -F',' -v OFS=$'\t' \
    'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6,$7,$8}' \
    "${SAMPLESHEET}")
N=${#SAMPLE_ROWS[@]}
echo "samplesheet: ${N} total RNA-seq samples"
[[ "${N}" -gt 0 ]] || { echo "ERROR: no rnaseq samples found in ${SAMPLESHEET}" >&2; exit 1; }
echo

# FastQC strips a recognized fastq extension off the INPUT file's own basename
# and appends _fastqc.zip -- reproduce that so the check looks for the right name.
fastqc_zip_name() {
    local b
    b="$(basename "$1")"
    b="${b%.fastq.gz}"; b="${b%.fq.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    printf '%s_fastqc.zip' "${b}"
}

###########
# 1. Slurm accounting -- did every array task actually finish?
###########
echo "== sacct summary since ${SACCT_SINCE} (jobname,state -> task count) =="
if command -v sacct >/dev/null 2>&1; then
    sacct -u "$USER" -S "${SACCT_SINCE}" \
        --format=JobID,JobName%22,State,ExitCode -n -P \
        | awk -F'|' '$2 ~ /_totalRNAseq$/ && $1 !~ /\./ {state[$2","$3]++} END {for (k in state) print k, state[k]}' \
        | sort
else
    echo "sacct not found on PATH -- skip (are you on a login node?)"
fi
echo

###########
# 2. Per-sample file completeness
###########
printf '%-20s %-6s %-6s %-6s %-6s %-6s\n' "sample" "fastp" "bwa" "ciri2" "fastqc" "flagst"
n_fail=0
n_fastp_legacy=0
FAILED_SAMPLES=()
for ROW in "${SAMPLE_ROWS[@]}"; do
    IFS=$'\t' read -r SM RAW_R1 RAW_R2 <<< "${ROW}"
    ok_fastp="FAIL"; ok_bwa="FAIL"; ok_ciri2="FAIL"; ok_fastqc="FAIL"; ok_flagstat="FAIL"

    if [[ -s "${TRIMDIR}/${SM}_R1.fq.gz" && -s "${TRIMDIR}/${SM}_R2.fq.gz" ]]; then
        if [[ -s "${FASTP_QCDIR}/${SM}.fastp.json" ]]; then
            ok_fastp="ok"
        elif [[ -s "${FASTP_QCDIR_LEGACY}/${SM}.fastp.json" ]]; then
            ok_fastp="ok"; n_fastp_legacy=$((n_fastp_legacy+1))
        fi
    fi

    if [[ -s "${ALIGNDIR}/${SM}.sam" && -s "${FLAGSTATDIR}/${SM}.flagstat.txt" ]] \
        && grep -q 'primary mapped' "${FLAGSTATDIR}/${SM}.flagstat.txt" 2>/dev/null; then
        ok_bwa="ok"
    fi

    if [[ -s "${CIRI2DIR}/${SM}.ciri" && -s "${CIRI2DIR}/${SM}.bsj2.bed" ]]; then
        if [[ -s "${CIRI2DIR}/${SM}.CIRI2.log" ]] \
            && grep -qiE 'error|died|can.t (locate|open)' "${CIRI2DIR}/${SM}.CIRI2.log"; then
            ok_ciri2="FAIL"
        else
            ok_ciri2="ok"
        fi
    fi

    zip1="$(fastqc_zip_name "${RAW_R1}")"
    zip2="$(fastqc_zip_name "${RAW_R2}")"
    [[ -s "${FASTQCDIR}/${zip1}" && -s "${FASTQCDIR}/${zip2}" ]] && ok_fastqc="ok"

    if [[ -s "${FLAGSTATDIR}/mapping_summary.tsv" ]] \
        && awk -F'\t' -v sm="${SM}" '$1==sm{f=1} END{exit !f}' "${FLAGSTATDIR}/mapping_summary.tsv"; then
        ok_flagstat="ok"
    fi

    printf '%-20s %-6s %-6s %-6s %-6s %-6s\n' "${SM}" "${ok_fastp}" "${ok_bwa}" "${ok_ciri2}" "${ok_fastqc}" "${ok_flagstat}"

    # fastqc/flagstat are single batch jobs over everything, not per-sample
    # array tasks -- a gap there is more likely "job hasn't run yet" than a
    # per-sample failure, so only fastp/bwa/ciri2 count toward pass/fail here.
    if [[ "${ok_fastp}" == FAIL || "${ok_bwa}" == FAIL || "${ok_ciri2}" == FAIL ]]; then
        n_fail=$((n_fail+1))
        FAILED_SAMPLES+=("${SM}")
    fi
done
if [[ "${n_fastp_legacy}" -gt 0 ]]; then
    echo
    echo "NOTE: found fastp reports only at the pre-reorg location for" \
         "${n_fastp_legacy} sample(s) -- consider moving them into ${QC_DIR}/trimmed to finish the migration."
fi
echo

###########
# 3. Scan stderr logs for real errors (skip normal stderr chatter, e.g. `which`/version prints)
###########
echo "== error markers in ${SLURM_OUT}/*_totalRNAseq_*.err =="
if [[ -d "${SLURM_OUT}" ]]; then
    grep -lEi 'error|traceback|core dumped|cancelled|oom|out of memory|due to time limit|segmentation fault' \
        "${SLURM_OUT}"/*_totalRNAseq_*.err 2>/dev/null | sort
    echo "(files listed above have at least one hit -- open just those, not all logs)"
else
    echo "slurm_out dir not found: ${SLURM_OUT}"
fi
echo

###########
# Summary
###########
echo "----------------------------------------------------------------------"
if [[ "${n_fail}" -eq 0 ]]; then
    echo "RESULT: all ${N} samples completed fastp -> bwa -> ciri2."
    echo "        Check the fastqc/flagstat columns above for the two batch QC jobs."
    exit 0
else
    echo "RESULT: ${n_fail}/${N} sample(s) missing output somewhere in fastp/bwa/ciri2:"
    printf '  %s\n' "${FAILED_SAMPLES[@]}"
    exit 1
fi
