#!/bin/bash
#SBATCH -J ciri_de
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri_de_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri_de_%j.err
#SBATCH -t 00:30:00
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 1
#SBATCH --mem=8GB

# ---------------------------------------------------------------------------
# STEP 7b (differential expression) of the circRNA-seq pipeline -- see
# circRNA_analysis_pipeline_outline.txt, STEP 7 Downstream.
#
# Runs CIRI_DE_replicate (CIRIquant's bundled edgeR wrapper), tumor vs normal,
# on both matrices from build_ciriquant_matrix.sbatch (run that first):
#   corrected    RNase R-adjusted expression -- PRIMARY result.
#   uncorrected  raw circRNA-seq expression  -- SENSITIVITY check.
# Grouping (T=tumor / C=normal) and pairing (patient_id) live in sample.lst,
# written by build_ciriquant_matrix.sbatch from master_samplesheet.csv.
# CIRI_DE_replicate fits circRNA BSJ counts jointly with host-gene expression
# (--gene), so a circRNA is only called differential beyond what its host
# gene's own expression change would predict.
#
# R/edgeR live in their OWN env ("ciriDE"), not inside the CIRIquant env --
# even py2's pysam/libcurl/openssl stack conflicts with the ancient
# krb5/openssl chain the CIRIquant-docs-recommended r-base=3.6 pin drags in,
# so a same-env `conda install` fails to solve. It doesn't need to be
# same-env anyway: CIRI_DE_replicate (CIRIquant/replicate.py) just does
# `subprocess.call('Rscript .../CIRI_DE.R ...')` -- it picks up whatever
# Rscript is first on $PATH at call time, so this script prepends the ciriDE
# env's bin/ to PATH right before invoking it (see below), same idea as the
# separate dupradar/sortmerna envs already in this pipeline.
#
# Create it once, if the preflight check below fails (unpinned -- let conda
# pick a mutually-compatible modern set rather than the old CIRIquant-docs
# versions that caused the same-env conflict).
#
# Channel order matters: conda-forge MUST come before bioconda, or the
# solver pulls r-base from bioconda's old GCC4-linked build, which crashes
# at runtime with "error while loading shared libraries: libgfortran.so.3:
# cannot open shared object file" (that ABI doesn't exist on modern systems).
# bioconductor-edger/limma only exist on bioconda, so they still resolve
# from there regardless of order; only r-base's toolchain needs conda-forge
# to win:
#   conda create -n ciriDE -c conda-forge -c bioconda \
#       r-base bioconductor-edger bioconductor-limma r-statmod r-optparse
#
# Input : ${OUTDIR}/<mode>/{library_info.csv, circRNA_bsj.csv, gene_count_matrix.csv}
#           (from build_ciriquant_matrix.sbatch)
# Output: ${OUTDIR}/<mode>/circRNA_de.tsv   per-locus differential BSJ test (logFC, p-value, FDR)
#         ${OUTDIR}/<mode>/gene_de.tsv      per-gene DE, byproduct of the same edgeR fit
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
# CIRIquant_github_env: official CIRIquant v1.1.3 install (python=2.7.15) --
# provides the CIRI_DE_replicate entry point. set +u around activation only --
# see run_ciriquant_circRNA_py2.sbatch header for why.
set +u
conda activate CIRIquant_github_env
set -u

echo ">> Preflight: checking 'ciriDE' conda env (R + edgeR)"
if ! conda run -n ciriDE Rscript -e 'suppressMessages(library(edgeR))' >/dev/null 2>&1; then
    echo "ERROR: conda env 'ciriDE' with R + edgeR was not found or is broken." >&2
    echo "       Create it with:" >&2
    echo "         conda create -n ciriDE -c conda-forge -c bioconda \\" >&2
    echo "             r-base bioconductor-edger bioconductor-limma r-statmod r-optparse" >&2
    exit 1
fi
conda run -n ciriDE Rscript -e 'cat("edgeR", as.character(packageVersion("edgeR")), "\n")'
echo "   ciriDE env OK"

# Make that env's Rscript resolve first on PATH -- CIRI_DE_replicate (still
# run as the CIRIquant_github_env python entry point below) shells out to
# plain `Rscript`, inheriting this process's PATH, so this is enough; no need
# to activate/switch envs for the actual call.
CIRIDE_BASE="$(conda info --base)/envs/ciriDE/bin"
[[ -d "${CIRIDE_BASE}" ]] || { echo "ERROR: ${CIRIDE_BASE} not found" >&2; exit 1; }
export PATH="${CIRIDE_BASE}:${PATH}"
echo ____________________________________
which CIRI_DE_replicate Rscript
echo ____________________________________

###########
# Variables
###########
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
OUTDIR="${ANALYSIS_DIR}/ciri_de"

###########
# CIRI_DE_replicate per mode
###########
for MODE in corrected uncorrected; do
    MODE_DIR="${OUTDIR}/${MODE}"
    LIB="${MODE_DIR}/library_info.csv"
    BSJ="${MODE_DIR}/circRNA_bsj.csv"
    GENE="${MODE_DIR}/gene_count_matrix.csv"

    MISSING=0
    for f in "${LIB}" "${BSJ}" "${GENE}"; do
        if [[ ! -s "${f}" ]]; then
            echo "${MODE}: missing ${f} -- run build_ciriquant_matrix.sbatch first, skipping" >&2
            MISSING=1
        fi
    done
    [[ "${MISSING}" -eq 0 ]] || continue

    echo "=== ${MODE}: CIRI_DE_replicate ==="
    CIRI_DE_replicate \
        --lib  "${LIB}" \
        --bsj  "${BSJ}" \
        --gene "${GENE}" \
        --out  "${MODE_DIR}/circRNA_de.tsv" \
        --out2 "${MODE_DIR}/gene_de.tsv"

    N_LOCI=$(($(wc -l < "${MODE_DIR}/circRNA_de.tsv") - 1))
    echo "${MODE}: done -> ${MODE_DIR}/circRNA_de.tsv (${N_LOCI} circRNA loci tested)"
    echo "  columns: $(head -1 "${MODE_DIR}/circRNA_de.tsv")"
    echo "  check the FDR/adjusted-p column above for significance -- not assumed by this script"
done
