#!/bin/bash
#SBATCH -J dupradar_totalRNAseq
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/dupradar_totalRNAseq_%j.out  # Output file with the job ID
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/dupradar_totalRNAseq_%j.err  # Error file with the job ID
#SBATCH -t 12:00:00     # whole subset in one job (D-HH:MM:SS); HISAT2 whole-genome align is the slow part
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 16      # 16 cores -> hisat2 -p 16 / samtools -@ 16 / featureCounts (inside dupRadar) -p 16
#SBATCH --mem=48GB

## Submit:  sbatch run_dupradar_totalRNAseq.sbatch
## Re-do every sample from scratch:  FORCE=1 sbatch run_dupradar_totalRNAseq.sbatch
## Run on specific samples instead of the auto-picked subset:
##   SAMPLES="TP00471_tumor_rnaseq TP00471_normal_rnaseq" sbatch run_dupradar_totalRNAseq.sbatch
## Widen/narrow the auto-picked subset (default 3 tumor/normal pairs = 6 samples):
##   N_PATIENTS=5 sbatch run_dupradar_totalRNAseq.sbatch

# ---------------------------------------------------------------------------
# dupRadar diagnostic for the ribo-depleted total RNA-seq ("rnaseq" modality)
# libraries -- distinguishes expression-driven (biological) duplication from
# PCR/technical duplication, which the raw FastQC duplication % cannot do.
#
# Why a representative SUBSET, not all 65 samples: dupRadar needs a real
# splice-aware whole-genome alignment (HISAT2), unlike the BWA-MEM SAM this
# repo already builds for CIRI2 (STEP 6b) -- BWA-MEM is not splice-aware and
# would misassign/over-split intron-spanning reads, corrupting dupRadar's
# per-gene duplication model. That alignment is the expensive part, and the
# duplication-vs-expression *pattern* this diagnostic is after is a property
# of the library prep protocol, not expected to vary meaningfully sample to
# sample -- so a handful of samples spanning tumor/normal is enough to answer
# "is this duplication biological or technical", without aligning the whole
# cohort just for a QC question. Default: auto-picks the first N_PATIENTS
# (default 3) patients in the master samplesheet that have both a tumor and a
# normal rnaseq sample -- 6 samples total. Override with SAMPLES= to name
# exact sample_ids instead.
#
# Do NOT treat this script's output BAMs as a substitute for the real
# quantification alignment (nf-core/rnaseq Stage 1a, STAR+Salmon) -- this is a
# side diagnostic, not part of the production BAM/count path.
#
# Strandedness caveat (see circRNA_analysis_pipeline_outline.txt STEP 6a /
# README "Before you run" #6): the rnaseq-modality libraries' strandedness has
# NOT been confirmed the way the circRNA-seq libraries' was (reverse,
# STEP 0). This script aligns unstranded by default (no --rna-strandness) and
# uses STRANDED=0 for dupRadar below. If you've since run
# run_infer_strandedness_test.sbatch against an rnaseq sample, set STRANDED
# to match (0 unstranded / 1 forward / 2 reverse -- same convention as
# CIRIquant's --library-type) for a more accurate per-gene duplication call.
#
# Requires a NEW conda env not yet in this repo's README (fastqc/fastp/bwa/
# hisat2/stringtie/samtools/ciriquant/seqtk live in 'circrna'; RSeQC lives in
# 'rseqc'): a Bioconductor env for dupRadar itself, e.g.:
#   conda create -n dupradar -c bioconda -c conda-forge bioconductor-dupradar
# The preflight check below fails fast with this same instruction if missing.
#
# Output (resumable, one sample at a time, skips samples already done):
#   ${DUPRADAR_DIR}/<sample>.markdup.bam(.bai)      HISAT2 BAM, duplicates MARKED (not removed)
#   ${DUPRADAR_DIR}/<sample>.dupradar.tsv           per-gene dupRadar table
#   ${DUPRADAR_DIR}/<sample>.dupradar_fit.txt       duprateExpFit() intercept/slope
#   ${DUPRADAR_DIR}/<sample>.dupradar_plot.pdf      expression-vs-duplication-rate plot (the main thing to eyeball)
#   ${DUPRADAR_DIR}/dupradar_summary.tsv            cohort summary: sample, intercept, slope
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
conda activate circrna     # hisat2, samtools

echo ____________________________________
which hisat2 ; hisat2 --version | head -1
which samtools ; samtools --version | head -1
echo ____________________________________
echo ">> Preflight: checking 'dupradar' conda env"
if ! conda run -n dupradar Rscript -e 'suppressMessages(library(dupRadar))' >/dev/null 2>&1; then
    echo "ERROR: conda env 'dupradar' with the Bioconductor dupRadar package was not found." >&2
    echo "       Create it with:" >&2
    echo "         conda create -n dupradar -c bioconda -c conda-forge bioconductor-dupradar" >&2
    exit 1
fi
echo "   dupradar env OK"
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
MODALITY="rnaseq"

REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
GTF="${REF_DIR}/gencode/gencode.v50.annotation.gtf"
HISAT_INDEX="${REF_DIR}/hg38/GRCh38.primary_assembly.genome_hisat2"   # already built (run_infer_strandedness_test.sbatch)

ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
TRIMDIR="${TRIMDIR:-${ANALYSIS_DIR}/total_RNAseq/trimmed}"
DUPRADAR_DIR="${ANALYSIS_DIR}/total_RNAseq/QC/dupradar"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
FORCE="${FORCE:-0}"
N_PATIENTS="${N_PATIENTS:-3}"
STRANDED="${STRANDED:-0}"   # 0 unstranded / 1 forward / 2 reverse -- see caveat above

mkdir -p "${DUPRADAR_DIR}"

[[ -s "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}" >&2; exit 1; }
[[ -s "${HISAT_INDEX}.1.ht2" ]] || { echo "ERROR: HISAT2 index not found: ${HISAT_INDEX} (run run_infer_strandedness_test.sbatch first)" >&2; exit 1; }

###########
# 1. Sample list: explicit SAMPLES= override, else auto-pick tumor/normal pairs
###########
if [[ -n "${SAMPLES:-}" ]]; then
    read -r -a SAMPLE_LIST <<< "${SAMPLES}"
else
    #   1 patient_id  3 tissue  5 modality  6 sample_id  7/8 fastq_1/2
    mapfile -t ROWS < <(
        awk -F',' -v m="${MODALITY}" \
            'NR>1 && $5==m && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $1","$3","$6}' \
            "${SAMPLESHEET}"
    )
    declare -A COUNT TUMOR_ID NORMAL_ID
    ORDER=()
    for row in "${ROWS[@]}"; do
        IFS=',' read -r pid tissue sid <<< "${row}"
        [[ -n "${COUNT[$pid]:-}" ]] || ORDER+=("${pid}")
        COUNT[$pid]=$(( ${COUNT[$pid]:-0} + 1 ))
        [[ "${tissue}" == "tumor"  ]] && TUMOR_ID[$pid]="${sid}"
        [[ "${tissue}" == "normal" ]] && NORMAL_ID[$pid]="${sid}"
    done

    SAMPLE_LIST=()
    picked=0
    for pid in "${ORDER[@]}"; do
        if [[ -n "${TUMOR_ID[$pid]:-}" && -n "${NORMAL_ID[$pid]:-}" ]]; then
            SAMPLE_LIST+=("${TUMOR_ID[$pid]}" "${NORMAL_ID[$pid]}")
            picked=$((picked+1))
            [[ "${picked}" -ge "${N_PATIENTS}" ]] && break
        fi
    done
fi

N="${#SAMPLE_LIST[@]}"
echo "Running dupRadar on ${N} sample(s): ${SAMPLE_LIST[*]}"
[[ "${N}" -gt 0 ]] || { echo "ERROR: no samples selected -- check SAMPLESHEET / SAMPLES / N_PATIENTS" >&2; exit 1; }

###########
# R script used for every sample (written once, reused)
###########
RSCRIPT="${DUPRADAR_DIR}/_run_dupradar.R"
cat > "${RSCRIPT}" <<'RSCRIPT_EOF'
suppressMessages(library(dupRadar))
args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]; gtf <- args[2]; stranded <- as.numeric(args[3])
threads <- as.numeric(args[4]); out_prefix <- args[5]

dm <- analyzeDuprates(bam, gtf, stranded, paired = TRUE, threads = threads)
write.table(dm, file = paste0(out_prefix, ".dupradar.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fit <- duprateExpFit(DupMat = dm)
sink(paste0(out_prefix, ".dupradar_fit.txt"))
cat("intercept (2^, approx. duplication fraction at ~zero expression):",
    2^fit$intercept, "\n")
cat("slope:", fit$slope, "\n")
sink()

pdf(paste0(out_prefix, ".dupradar_plot.pdf"))
duprateExpDensPlot(DupMat = dm)
title("Expression vs. duplication rate")
dev.off()
RSCRIPT_EOF

###########
# 2. Loop: HISAT2 align -> sort -> mark duplicates (samtools, not removed) -> dupRadar
###########
n_done=0; n_skip=0; n_fail=0
for SM in "${SAMPLE_LIST[@]}"; do
    OUT_PREFIX="${DUPRADAR_DIR}/${SM}"
    MARKED_BAM="${OUT_PREFIX}.markdup.bam"
    FIT="${OUT_PREFIX}.dupradar_fit.txt"

    if [[ "${FORCE}" != "1" && -s "${FIT}" ]]; then
        echo "SKIP  ${SM} (already done)"
        n_skip=$((n_skip+1))
        continue
    fi

    R1="${TRIMDIR}/${SM}_R1.fq.gz"
    R2="${TRIMDIR}/${SM}_R2.fq.gz"
    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "MISSING trimmed FASTQ for ${SM} (${R1} / ${R2}) -- skipping" >&2
        n_fail=$((n_fail+1))
        continue
    fi

    echo "RUN   ${SM}"

    ###########
    # HISAT2 align, coordinate-sorted BAM
    ###########
    RAW_BAM="${DUPRADAR_DIR}/${SM}.sorted.bam"
    if [[ ! -s "${RAW_BAM}" ]]; then
        hisat2 -p "${THREADS}" -x "${HISAT_INDEX}" \
               -1 "${R1}" -2 "${R2}" \
            | samtools sort -@ "${THREADS}" -o "${RAW_BAM}" -
    fi

    ###########
    # Mark duplicates with samtools only (no Picard in this repo's envs):
    # name-sort -> fixmate (adds mate-score tags markdup needs) -> coord-sort -> markdup
    ###########
    NAME_SORTED="${DUPRADAR_DIR}/${SM}.namesorted.bam"
    FIXMATE="${DUPRADAR_DIR}/${SM}.fixmate.bam"
    RESORTED="${DUPRADAR_DIR}/${SM}.resorted.bam"

    samtools sort -n -@ "${THREADS}" -o "${NAME_SORTED}" "${RAW_BAM}"
    samtools fixmate -m -@ "${THREADS}" "${NAME_SORTED}" "${FIXMATE}"
    samtools sort -@ "${THREADS}" -o "${RESORTED}" "${FIXMATE}"
    samtools markdup -@ "${THREADS}" "${RESORTED}" "${MARKED_BAM}"   # marks only -- no -r, dupRadar needs the flag set
    samtools index "${MARKED_BAM}"

    rm -f "${NAME_SORTED}" "${FIXMATE}" "${RESORTED}"

    ###########
    # dupRadar
    ###########
    if conda run --no-capture-output -n dupradar Rscript "${RSCRIPT}" \
            "${MARKED_BAM}" "${GTF}" "${STRANDED}" "${THREADS}" "${OUT_PREFIX}"; then
        n_done=$((n_done+1))
        echo "${SM}: $(cat "${FIT}" | tr '\n' ' ')"
    else
        echo "FAIL  ${SM} (dupRadar)" >&2
        n_fail=$((n_fail+1))
    fi
done

echo "----------------------------------------------------------------------"
echo "dupRadar: ${n_done} generated, ${n_skip} skipped, ${n_fail} failed"

###########
# 3. Cohort summary
###########
SUMMARY="${DUPRADAR_DIR}/dupradar_summary.tsv"
{
    printf 'sample\tintercept_frac\tslope\n'
    for f in "${DUPRADAR_DIR}"/*.dupradar_fit.txt; do
        [[ -e "${f}" ]] || continue
        sm="$(basename "${f}")"; sm="${sm%.dupradar_fit.txt}"
        icpt=$(awk -F': ' '/intercept/{print $2}' "${f}")
        slope=$(awk -F': ' '/slope/{print $2}' "${f}")
        printf '%s\t%s\t%s\n' "${sm}" "${icpt:-NA}" "${slope:-NA}"
    done
} > "${SUMMARY}"

echo "wrote ${SUMMARY}"
column -t -s $'\t' "${SUMMARY}" || cat "${SUMMARY}"

cat <<'EOF'

============================================================
 HOW TO READ THIS
============================================================
 The PDF plot per sample is the primary diagnostic -- eyeball it first:
   - Duplication rate rising smoothly with expression, low at low-expression
     genes -> largely BIOLOGICAL duplication (expected for total RNA-seq).
     Do not deduplicate; proceed with quantification/DE as planned.
   - Duplication rate high even for LOW-expression genes (a high, fairly flat
     curve near the y-axis) -> a real PCR/technical duplication signature.
     Without UMIs this isn't fixable in the data already sequenced, but it's
     useful to know for interpreting effect sizes and for future library prep.

 intercept_frac is a rough single-number summary of the same thing: the
 duplication fraction the fitted curve extrapolates to at ~zero expression.
 Small (e.g. <0.05-0.1) supports the biological-duplication read; large points
 toward a technical component. Treat this as a rough numeric anchor for the
 plot, not a hard pass/fail threshold.

 Compare across the tumor/normal pairs run here -- a consistent pattern across
 samples further supports "this is how the total RNA-seq libraries behave",
 rather than a one-off outlier sample.
============================================================
EOF

echo "done"
