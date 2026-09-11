#!/bin/bash
#SBATCH -J flagstat_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -t 02:00:00     # whole batch in one job (D-HH:MM:SS); flagstat is a single sequential pass per file
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8       # 8 cores -> samtools flagstat -@ 8
#SBATCH --mem=8GB       # flagstat streams the file, needs almost no RAM
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/mapping_stats_totalRNAseq_%j.out  # Output file with the job ID
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/mapping_stats_totalRNAseq_%j.err  # Error file with the job ID


## Submit:  sbatch run_flagstat_totalRNAseq.sbatch
## Re-do every report from scratch:  FORCE=1 sbatch run_flagstat_totalRNAseq.sbatch
## Point at a different folder of alignments:
##   SAMDIR=/scratch/alpine/$USER/.../bwa_mapped sbatch run_flagstat_totalRNAseq.sbatch

# ---------------------------------------------------------------------------
# Alignment QC for the total RNA-seq validation branch: run `samtools flagstat`
# on every SAM (or BAM/CRAM) in SAMDIR and collect the mapping stats.
# Twin of run_flagstat_circRNA.sbatch; only SAMDIR and the slurm_out names differ.
# run_bwa_totalRNAseq.sbatch already writes a per-sample report inline, so on a
# clean re-run this job mostly just (re)builds the summary table.
#
#   - one report per sample under ${SAMDIR}/flagstat/<sample>.flagstat.txt
#   - resumable, samples whose report already exists AND looks complete are skipped
#   - a summary table is (re)built from ALL reports every run, one row per
#     sample: ${SAMDIR}/flagstat/mapping_summary.tsv
#         sample  total_reads  mapped_reads  mapping_pct
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load samtools

###########
# Variables
###########
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"

# Folder holding the alignments to profile. Default: the BWA-MEM SAMs from
# run_bwa_totalRNAseq.sbatch. Override on the command line with SAMDIR=... sbatch run_flagstat_totalRNAseq.sbatch
SAMDIR="${SAMDIR:-${ANALYSIS_DIR}/bwa_mapped/total_RNAseq}"

# File extension to look for (sam | bam | cram).
EXT="${EXT:-sam}"

# FORCE=1 -> ignore existing reports and regenerate all of them.
FORCE="${FORCE:-0}" #skip if report already exsists

THREADS="${SLURM_CPUS_PER_TASK:-8}"

STATSDIR="${SAMDIR}/flagstat"
SUMMARY="${STATSDIR}/mapping_summary.tsv"

[[ -d "${SAMDIR}" ]] || { echo "ERROR: SAMDIR not found: ${SAMDIR}" >&2; exit 1; }
mkdir -p "${STATSDIR}"

###########
# Helper: does an existing report look like a finished flagstat run?
###########
# `samtools flagstat` computes everything, then prints ~13-16 lines in one go.
# A job killed by a timeout can leave an empty or truncated file, so require
# both the "primary mapped" line and the final "different chr" line.
report_is_complete() {
    local f="$1"
    [[ -s "${f}" ]] || return 1
    grep -q 'primary mapped'  "${f}" || return 1
    grep -q 'mate mapped to a different chr' "${f}" || return 1
    return 0
}

###########
# Run flagstat over every alignment file
###########
shopt -s nullglob
FILES=( "${SAMDIR}"/*."${EXT}" )
shopt -u nullglob

(( ${#FILES[@]} > 0 )) || { echo "ERROR: no *.${EXT} files in ${SAMDIR}" >&2; exit 1; }
echo "found ${#FILES[@]} *.${EXT} file(s) in ${SAMDIR}"

n_done=0; n_skip=0; n_fail=0
for f in "${FILES[@]}"; do
    sm="$(basename "${f}")"; sm="${sm%.${EXT}}"
    out="${STATSDIR}/${sm}.flagstat.txt"

    if [[ "${FORCE}" != "1" ]] && report_is_complete "${out}"; then
        echo "SKIP  ${sm} (report already complete)"
        n_skip=$((n_skip+1))
        continue
    fi

    echo "RUN   ${sm}"
    # write to a temp file so an interrupted run never leaves a half report
    # that would later be mistaken for complete
    tmp="${out}.partial.$$"
    if samtools flagstat -@ "${THREADS}" "${f}" > "${tmp}" && report_is_complete "${tmp}"; then
        mv -f "${tmp}" "${out}"
        n_done=$((n_done+1))
        grep -E 'primary mapped|properly paired' "${out}" | sed "s/^/      ${sm}: /" || true
    else
        rm -f "${tmp}"
        echo "FAIL  ${sm} (flagstat did not finish)" >&2
        n_fail=$((n_fail+1))
    fi
done

echo "----------------------------------------------------------------------"
echo "flagstat: ${n_done} generated, ${n_skip} skipped, ${n_fail} failed"

###########
# Build the combined mapping table from ALL reports -- one row per sample
###########
# Columns: sample  total_reads  mapped_reads  mapping_pct   (primary reads only)
{
    printf 'sample\ttotal_reads\tmapped_reads\tmapping_pct\n'
    for r in "${STATSDIR}"/*.flagstat.txt; do
        [[ -e "${r}" ]] || continue
        sm="$(basename "${r}")"; sm="${sm%.flagstat.txt}"
        awk -v sm="${sm}" '
            /[0-9]+ \+ [0-9]+ primary$/          { tot = $1 }
            /[0-9]+ \+ [0-9]+ primary mapped \(/ { mapped = $1; pct = $6; gsub(/[(%]/,"",pct) }
            END {
                if (tot == "")    tot = "NA"
                if (mapped == "") { mapped = "NA"; pct = "NA" }
                printf "%s\t%s\t%s\t%s\n", sm, tot, mapped, pct
            }' "${r}"
    done
} > "${SUMMARY}"

echo "wrote ${SUMMARY}"
column -t -s $'\t' "${SUMMARY}" || cat "${SUMMARY}"

# Optional: fold every report into one MultiQC section
#   module load multiqc 2>/dev/null || conda activate circrna
#   multiqc -f -o "${STATSDIR}" "${STATSDIR}"

echo "done"
