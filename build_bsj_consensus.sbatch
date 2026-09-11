#!/bin/bash
#SBATCH -J bsj_consensus
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bsj_consensus_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bsj_consensus_%j.err
#SBATCH -t 00:30:00    # whole cohort, seconds of real work
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 1
#SBATCH --mem=4GB


# ---------------------------------------------------------------------------
# STEP 4 (consensus) of the circRNA-seq pipeline -- see
# circRNA_analysis_pipeline_outline.txt
#
# Merge the per-tool BSJ call sets for each sample into one consensus locus
# set, tag every locus with how many (and which) tools called it, and keep
# those passing MIN_TOOLS. Then slop-merge every sample's consensus into one
# cohort-wide union BED -- the fixed locus set that run_ciriquant_circRNA.sbatch
# quantifies in every sample (so the DE matrix has consistent rows).
#
# NOW: CIRI2 is the only caller wired in, so with MIN_TOOLS=1 this just
# collapses CIRI2's own >=2-read loci (merging duplicates within SLOP bp) and
# passes them through. When CIRCexplorer2 / find_circ are added:
#   1. uncomment their entries in TOOL_BED below
#   2. set MIN_TOOLS=2  -- the >=2-of-3 precision filter from the outline
#
# Input  (per sample, BED6: chr start[0-based] end name bsj_reads strand):
#   ${CIRI2DIR}/${SM}.bsj2.bed        from run_ciri2_circRNA.sbatch (circRNAseq/ciri2/)
# Output (${OUTDIR} = ${ANALYSIS_DIR}/consensus):
#   ${SM}.consensus.bed              per-sample consensus, BED6+2
#                                    (extra cols: n_tools, tools_csv)
#   cohort.consensus.union.bed       BED6 union across all samples -> CIRIquant --bed
#   consensus_summary.tsv            per-sample locus counts
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/circRNAseq"
CIRI2DIR="${ANALYSIS_DIR}/circRNAseq/ciri2"
OUTDIR="${ANALYSIS_DIR}/consensus"

MIN_TOOLS="${MIN_TOOLS:-1}"   # bump to 2 once a second caller is wired in
SLOP="${SLOP:-5}"             # bp margin for matching BSJs across tools / samples

# Per-tool normalised BSJ BEDs. "__SM__" is replaced with the sample id.
# Every file must be BED6 (chr, start 0-based, end, name, bsj_reads, strand).
declare -A TOOL_BED=(
    [ciri2]="${CIRI2DIR}/__SM__.bsj2.bed"
    # [circexplorer2]="${ANALYSIS_DIR}/circexplorer2/__SM__.circexplorer2.bsj.bed"
    # [find_circ]="${ANALYSIS_DIR}/find_circ/__SM__.find_circ.bsj.bed"
)

mkdir -p "${OUTDIR}"

[[ -s "${SAMPLESHEET}" ]] || { echo "ERROR: missing samplesheet: ${SAMPLESHEET}" >&2; exit 1; }

###########
# Sample list -- identical filter to every other step, so ids line up
###########
mapfile -t SAMPLES < <(awk -F',' \
    'NR>1 && $5=="circrnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
    "${SAMPLESHEET}")
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no circRNA samples in ${SAMPLESHEET}" >&2; exit 1; }
echo "consensus: ${#SAMPLES[@]} samples, MIN_TOOLS=${MIN_TOOLS}, SLOP=${SLOP}"

###########
# Stranded single-linkage merge (sorted input: chr, strand, start)
#   in : chr start end name reads strand tool
#   out: chr start end chr:start-end max_reads strand n_tools tools_csv
###########
merge_awk='
BEGIN {
    FS = OFS = "\t"
    NPRI = split("ciri2 circexplorer2 find_circ", PRI, " ")
    for (i = 1; i <= NPRI; i++) prilist[PRI[i]] = 1
}
function flush(   i, n, list, t) {
    if (!have) return
    n = 0; list = ""
    for (i = 1; i <= NPRI; i++)
        if (PRI[i] in seen) { n++; list = list (list == "" ? "" : ",") PRI[i] }
    for (t in seen)
        if (!(t in prilist)) { n++; list = list (list == "" ? "" : ",") t }
    if (n >= mintools)
        print cchr, cstart, cend, cchr ":" cstart "-" cend, cmax, cstrand, n, list
}
{
    chr = $1; s = $2 + 0; e = $3 + 0; reads = $5 + 0; strand = $6; tool = $7
    if (have && chr == cchr && strand == cstrand && s <= cend + slop) {
        if (e > cend)    cend = e
        if (reads > cmax) cmax = reads
        seen[tool] = 1
    } else {
        flush()
        have = 1
        cchr = chr; cstart = s; cend = e; cstrand = strand; cmax = reads
        delete seen; seen[tool] = 1
    }
}
END { flush() }
'

###########
# Per-sample consensus
###########
SUMMARY="${OUTDIR}/consensus_summary.tsv"
printf 'sample\tpooled_intervals\tconsensus_loci\tmin_tools\ttools_available\n' > "${SUMMARY}"

UNION_TMP="${OUTDIR}/.union.tmp"
: > "${UNION_TMP}"

for SM in "${SAMPLES[@]}"; do
    POOL="${OUTDIR}/.${SM}.pooled.tmp"
    : > "${POOL}"
    avail=()
    for tool in "${!TOOL_BED[@]}"; do
        f="${TOOL_BED[$tool]//__SM__/${SM}}"
        if [[ -s "${f}" ]]; then
            awk -v t="${tool}" 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=6 {print $1,$2,$3,$4,$5,$6,t}' \
                "${f}" >> "${POOL}"
            avail+=("${tool}")
        else
            echo "WARN ${SM}: no ${tool} BED (${f})" >&2
        fi
    done

    if [[ ${#avail[@]} -eq 0 ]]; then
        echo "ERROR ${SM}: no input BEDs found, skipping" >&2
        printf '%s\t0\t0\t%s\t-\n' "${SM}" "${MIN_TOOLS}" >> "${SUMMARY}"
        rm -f "${POOL}"
        continue
    fi

    OUT="${OUTDIR}/${SM}.consensus.bed"
    sort -k1,1 -k6,6 -k2,2n "${POOL}" \
        | awk -v slop="${SLOP}" -v mintools="${MIN_TOOLS}" "${merge_awk}" \
        | sort -k1,1 -k2,2n > "${OUT}"

    n_pool=$(wc -l < "${POOL}")
    n_cons=$(wc -l < "${OUT}")
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${SM}" "${n_pool}" "${n_cons}" "${MIN_TOOLS}" "$(IFS=,; echo "${avail[*]}")" >> "${SUMMARY}"
    echo "${SM}: ${n_pool} pooled -> ${n_cons} consensus (>=${MIN_TOOLS} tool/s; ${avail[*]})"

    cut -f1-6 "${OUT}" >> "${UNION_TMP}"
    rm -f "${POOL}"
done

###########
# Cohort union -> the locus set CIRIquant quantifies in every sample
###########
UNION_BED="${OUTDIR}/cohort.consensus.union.bed"
if [[ -s "${UNION_TMP}" ]]; then
    sort -k1,1 -k6,6 -k2,2n "${UNION_TMP}" \
        | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,"union"}' \
        | awk -v slop="${SLOP}" -v mintools=1 "${merge_awk}" \
        | cut -f1-6 \
        | sort -k1,1 -k2,2n > "${UNION_BED}"
else
    : > "${UNION_BED}"
fi
rm -f "${UNION_TMP}"

echo "----"
echo "cohort union: $(wc -l < "${UNION_BED}") loci -> ${UNION_BED}"
echo "summary     : ${SUMMARY}"
echo "next        : sbatch run_ciriquant_circRNA.sbatch"
