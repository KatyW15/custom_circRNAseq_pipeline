#!/bin/bash
#SBATCH -J bsj_validation
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bsj_validation_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bsj_validation_%j.err
#SBATCH -t 00:30:00    # whole cohort, seconds of real work
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 1
#SBATCH --mem=4GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=YOUR_EMAIL@example.com

#
# NOT part of job_inchain_totalRNAseq.sh -- run this after both
# build_bsj_consensus.sbatch (circRNA branch) and run_ciri2_totalRNAseq.sbatch.

# ---------------------------------------------------------------------------
# STEP 6b (orthogonal validation) -- see circRNA_analysis_pipeline_outline.txt
#
# Intersect the circRNA-seq consensus BSJ set with circRNAs detected de novo in
# the matched RNase R(-) total RNA-seq libraries. A consensus circRNA that is
# ALSO found in non-enriched data is much less likely to be an RNase R artifact
# (rolling-circle / template switching) -> per-circRNA confidence flag alongside
# the STEP 6a --RNaseR expression correction.
#
# Matching is a stranded single-linkage merge with +/- SLOP bp on the BSJ
# coordinates (same method as build_bsj_consensus.sbatch).
#
# Input:
#   circRNA-seq per-sample consensus  ${ANALYSIS_DIR}/consensus/<circ sample>.consensus.bed
#                                       (BED6+2: chr start end name reads strand n_tools tools_csv)
#   total RNA-seq CIRI2 loci          ${ANALYSIS_DIR}/ciri2/total_RNAseq/<rnaseq sample>.bsj2.bed
#                                       (BED6: chr start end name reads strand)
#
# Pairing (sample_id = <patient>_<tissue>_<modality>):
#   TP#####_tumor_rnaseq   <->   TP#####_tumor_circrnaseq
#
# Output (${ANALYSIS_DIR}/validation):
#   <circ sample>.validation.bed     per matched pair, one row per consensus locus:
#       chr start end locus circseq_bsj strand found_in_totalRNA(0/1) totalRNA_bsj
#   pair_overview.tsv                per pair: n_consensus, n_validated, pct
#   cohort.circRNA_vs_totalRNA.tsv   per consensus locus across the cohort:
#       chr start end locus strand n_circseq_samples n_totalRNA_samples validated(yes/no)
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
CONSENSUS_DIR="${ANALYSIS_DIR}/consensus"              # from build_bsj_consensus.sbatch
CIRI2_RNA_DIR="${ANALYSIS_DIR}/ciri2/total_RNAseq"     # from run_ciri2_totalRNAseq.sbatch
OUTDIR="${ANALYSIS_DIR}/validation"

SLOP="${SLOP:-5}"                       # bp margin for matching BSJs (match build_bsj_consensus)
RNA_MODALITY_SUFFIX="_rnaseq"
CIRC_MODALITY_SUFFIX="_circrnaseq"

mkdir -p "${OUTDIR}"
[[ -s "${SAMPLESHEET}" ]] || { echo "ERROR: missing samplesheet: ${SAMPLESHEET}" >&2; exit 1; }

###########
# rnaseq sample list -- same filter as run_*_totalRNAseq.sbatch
###########
mapfile -t RNA_SAMPLES < <(awk -F',' \
    'NR>1 && $5=="rnaseq" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6}' \
    "${SAMPLESHEET}")
[[ ${#RNA_SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no rnaseq samples in ${SAMPLESHEET}" >&2; exit 1; }
echo "validation: ${#RNA_SAMPLES[@]} rnaseq samples, SLOP=${SLOP}"

###########
# Per-pair merge: circseq consensus vs total-RNA CIRI2
#   in : chr start end name reads strand cls        (cls = circseq | rnaseq)
#   out: chr start end chr:start-end circseq_bsj strand found_in_totalRNA totalRNA_bsj
###########
pair_awk='
BEGIN { FS = OFS = "\t" }
function flush() {
    if (!have) return
    if (has_c)
        print cchr, cstart, cend, cchr ":" cstart "-" cend, c_reads, cstrand, (has_r ? 1 : 0), r_reads
}
{
    chr = $1; s = $2 + 0; e = $3 + 0; reads = $5 + 0; strand = $6; cls = $7
    if (have && chr == cchr && strand == cstrand && s <= cend + slop) {
        if (e > cend) cend = e
    } else {
        flush()
        have = 1; cchr = chr; cstart = s; cend = e; cstrand = strand
        has_c = 0; has_r = 0; c_reads = 0; r_reads = 0
    }
    if (cls == "circseq") { has_c = 1; if (reads > c_reads) c_reads = reads }
    else                  { has_r = 1; if (reads > r_reads) r_reads = reads }
}
END { flush() }
'

###########
# Cohort merge: every circseq consensus locus + every total-RNA CIRI2 locus,
# tagged by sample, collapsed across the whole cohort.
#   in : chr start end name reads strand cls smp
#   out: chr start end chr:start-end strand n_circseq_samples n_totalRNA_samples validated
###########
cohort_awk='
BEGIN { FS = OFS = "\t" }
function flush(   k, nc, nr) {
    if (!have) return
    nc = 0; for (k in cS) nc++
    nr = 0; for (k in rS) nr++
    if (nc > 0)
        print cchr, cstart, cend, cchr ":" cstart "-" cend, cstrand, nc, nr, (nr > 0 ? "yes" : "no")
}
{
    chr = $1; s = $2 + 0; e = $3 + 0; strand = $6; cls = $7; smp = $8
    if (have && chr == cchr && strand == cstrand && s <= cend + slop) {
        if (e > cend) cend = e
    } else {
        flush()
        have = 1; cchr = chr; cstart = s; cend = e; cstrand = strand
        delete cS; delete rS
    }
    if (cls == "circseq") cS[smp] = 1; else rS[smp] = 1
}
END { flush() }
'

###########
# Per-pair validation
###########
OVERVIEW="${OUTDIR}/pair_overview.tsv"
printf 'circ_sample\trnaseq_sample\tn_consensus\tn_validated\tpct_validated\n' > "${OVERVIEW}"

COHORT_TMP="${OUTDIR}/.cohort.tmp"
: > "${COHORT_TMP}"

for RNA_SM in "${RNA_SAMPLES[@]}"; do
    CIRC_SM="${RNA_SM%${RNA_MODALITY_SUFFIX}}${CIRC_MODALITY_SUFFIX}"
    CONS="${CONSENSUS_DIR}/${CIRC_SM}.consensus.bed"
    RNA_BED="${CIRI2_RNA_DIR}/${RNA_SM}.bsj2.bed"

    if [[ ! -s "${CONS}" ]]; then
        echo "WARN ${RNA_SM}: no circseq consensus (${CONS}), skipping pair" >&2
        continue
    fi
    if [[ ! -s "${RNA_BED}" ]]; then
        echo "WARN ${RNA_SM}: no total-RNA CIRI2 BED (${RNA_BED}), skipping pair" >&2
        continue
    fi

    POOL="${OUTDIR}/.${RNA_SM}.pooled.tmp"
    {
        awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=6 {print $1,$2,$3,$4,$5,$6,"circseq"}' "${CONS}"
        awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=6 {print $1,$2,$3,$4,$5,$6,"rnaseq"}'  "${RNA_BED}"
    } > "${POOL}"

    OUT="${OUTDIR}/${CIRC_SM}.validation.bed"
    sort -k1,1 -k6,6 -k2,2n "${POOL}" \
        | awk -v slop="${SLOP}" "${pair_awk}" \
        | sort -k1,1 -k2,2n > "${OUT}"

    n_cons=$(wc -l < "${OUT}")
    n_val=$(awk 'BEGIN{FS="\t"} $7==1' "${OUT}" | wc -l)
    pct=$(awk -v a="${n_val}" -v b="${n_cons}" 'BEGIN{printf (b>0 ? "%.1f" : "NA"), (b>0 ? 100*a/b : 0)}')
    printf '%s\t%s\t%s\t%s\t%s\n' "${CIRC_SM}" "${RNA_SM}" "${n_cons}" "${n_val}" "${pct}" >> "${OVERVIEW}"
    echo "${CIRC_SM} <-> ${RNA_SM}: ${n_val}/${n_cons} consensus circRNAs also in total RNA (${pct}%)"

    # feed the cohort pool (tag each locus with its sample)
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=6 {print $1,$2,$3,$4,$5,$6,"circseq","'"${CIRC_SM}"'"}' "${CONS}"    >> "${COHORT_TMP}"
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=6 {print $1,$2,$3,$4,$5,$6,"rnaseq","'"${RNA_SM}"'"}'  "${RNA_BED}" >> "${COHORT_TMP}"
    rm -f "${POOL}"
done

###########
# Cohort-wide table
###########
COHORT="${OUTDIR}/cohort.circRNA_vs_totalRNA.tsv"
printf 'chr\tstart\tend\tlocus\tstrand\tn_circseq_samples\tn_totalRNA_samples\tvalidated\n' > "${COHORT}"
if [[ -s "${COHORT_TMP}" ]]; then
    sort -k1,1 -k6,6 -k2,2n "${COHORT_TMP}" \
        | awk -v slop="${SLOP}" "${cohort_awk}" \
        | sort -k1,1 -k2,2n >> "${COHORT}"
fi
rm -f "${COHORT_TMP}"

n_loci=$(($(wc -l < "${COHORT}") - 1))
n_val=$(awk 'BEGIN{FS="\t"} NR>1 && $8=="yes"' "${COHORT}" | wc -l)
echo "----"
echo "cohort: ${n_val}/${n_loci} consensus circRNA loci validated in >=1 total RNA-seq library"
echo "  per-pair beds : ${OUTDIR}/<circ sample>.validation.bed"
echo "  pair overview : ${OVERVIEW}"
echo "  cohort table  : ${COHORT}"
