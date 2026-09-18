#!/bin/bash
#SBATCH -J ciri_de_prep
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri_de_prep_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/ciri_de_prep_%j.err
#SBATCH -t 00:30:00
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 1
#SBATCH --mem=8GB

# ---------------------------------------------------------------------------
# STEP 7a (matrix build) of the circRNA-seq pipeline -- see
# circRNA_analysis_pipeline_outline.txt, STEP 7 Downstream.
#
# Single aggregation job (no array, cheap): merges every sample's already-
# computed CIRIquant per-locus GTF into cohort-wide matrices via prep_CIRIquant
# + prepDE.py, in two modes:
#   corrected    RNase R-adjusted expression (STEP 6a,
#                run_ciriquant_RNaseR_corrected_py2.sbatch). PRIMARY input for DE.
#                Only patients with a matched total RNA-seq mock are included;
#                circRNA-seq-only patients are skipped and logged.
#   uncorrected  raw circRNA-seq CIRIquant output (STEP 5,
#                run_ciriquant_circRNA_py2.sbatch), no RNase R correction.
#                SENSITIVITY check; includes every circRNA-seq sample.
# Tumor vs normal grouping (T/C) and per-patient pairing (patient_id) come
# straight from master_samplesheet.csv, so CIRI_DE_replicate
# (run_ciri_de.sbatch, next step) can run a paired edgeR test.
#
# Runs in CIRIquant_github_env (official py2 CIRIquant install) -- the patched
# bioconda/Python 3.11 `ciriquant` env is retired (kept repeatedly hitting new
# upstream Python-2-only bugs; py2 has been running cleanly). prep_CIRIquant
# and prepDE.py are unaffected either way (no Python 2/3 issues found in
# either during the proactive sweep), but every input path below now points
# at the py2 run's outputs, which is what actually matters here.
#
# Input : master_samplesheet.csv
#         STEP 5  circRNAseq/ciriquant_noRNaseRcorrection_py2/<circ sample>/{<circ sample>.gtf, gene/<circ sample>_out.gtf}
#         STEP 6a ciriquant_py2/<rnaseq sample>/{<rnaseq sample>.gtf, gene/<rnaseq sample>_out.gtf}
# Output: ${OUTDIR}/<mode>/sample.lst              name group(T/C) patient_id gtf_path   (CIRI_DE_replicate input)
#         ${OUTDIR}/<mode>/sample_gene.lst         name gtf_path                          (prepDE.py input)
#         ${OUTDIR}/<mode>/library_info.csv        prep_CIRIquant output; library sizes for DE
#         ${OUTDIR}/<mode>/circRNA_info.csv        prep_CIRIquant output; per-locus annotation
#         ${OUTDIR}/<mode>/circRNA_bsj.csv         prep_CIRIquant output; BSJ count matrix
#         ${OUTDIR}/<mode>/circRNA_ratio.csv       prep_CIRIquant output; junction-ratio matrix
#         ${OUTDIR}/<mode>/gene_count_matrix.csv   prepDE.py output; gene-level count matrix (CIRI_DE_replicate --gene)
#         ${OUTDIR}/<mode>/transcript_count_matrix.csv  prepDE.py output (unused downstream, kept)
# ---------------------------------------------------------------------------

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge
# CIRIquant_github_env: official CIRIquant v1.1.3 install (python=2.7.15).
# set +u around activation only -- see run_ciriquant_circRNA_py2.sbatch header
# for why (binutils_linux-64's activate.d script references an unset var).
set +u
conda activate CIRIquant_github_env
set -u

echo ____________________________________
which prep_CIRIquant prepDE.py
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"
ANALYSIS_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis"
UNCORR_DIR="${ANALYSIS_DIR}/circRNAseq/ciriquant_noRNaseRcorrection_py2"   # STEP 5
CORR_DIR="${ANALYSIS_DIR}/ciriquant_py2"                 # STEP 6a
OUTDIR="${ANALYSIS_DIR}/ciri_de"
RNA_MODALITY_SUFFIX="_rnaseq"                 # how the matched rnaseq sample_id ends

mkdir -p "${OUTDIR}/corrected" "${OUTDIR}/uncorrected"

[[ -s "${SAMPLESHEET}" ]] || { echo "ERROR: missing samplesheet: ${SAMPLESHEET}" >&2; exit 1; }

###########
# Build sample.lst + sample_gene.lst for both modes from the samplesheet.
# Matrix column names use the clean <patient>_<tissue> label (no modality
# suffix) so corrected/uncorrected matrices line up on the same names.
###########
python3 - "${SAMPLESHEET}" "${UNCORR_DIR}" "${CORR_DIR}" "${OUTDIR}" "${RNA_MODALITY_SUFFIX}" <<'PY'
import csv, os, sys

samplesheet, uncorr_dir, corr_dir, outdir, rna_suf = sys.argv[1:6]

group_code = {"tumor": "T", "normal": "C"}

with open(samplesheet, newline="") as f:
    rows = [r for r in csv.DictReader(f) if r["modality"] == "circrnaseq"]

uncorr_lines, uncorr_gene_lines = [], []
corr_lines, corr_gene_lines = [], []
n_corrected_skipped = 0

for row in rows:
    patient, tissue, circ_sm = row["patient_id"], row["tissue"], row["sample_id"]
    if tissue not in group_code:
        print(f"WARNING: {circ_sm} has unexpected tissue '{tissue}', skipping", file=sys.stderr)
        continue
    label = f"{patient}_{tissue}"          # clean matrix column name
    grp = group_code[tissue]

    # -- uncorrected (sensitivity): every circRNA-seq sample, STEP 5 output
    unc_gtf = os.path.join(uncorr_dir, circ_sm, f"{circ_sm}.gtf")
    unc_gene_gtf = os.path.join(uncorr_dir, circ_sm, "gene", f"{circ_sm}_out.gtf")
    if os.path.exists(unc_gtf) and os.path.getsize(unc_gtf) > 0:
        uncorr_lines.append(f"{label} {unc_gtf} {grp} {patient}")
        uncorr_gene_lines.append(f"{label} {unc_gene_gtf}")
    else:
        print(f"WARNING: missing uncorrected gtf for {circ_sm}: {unc_gtf} -- excluded from uncorrected matrix", file=sys.stderr)

    # -- corrected (primary): needs the matched rnaseq sample's STEP 6a output
    rna_sm = f"{patient}_{tissue}{rna_suf}"
    cor_gtf = os.path.join(corr_dir, rna_sm, f"{rna_sm}.gtf")
    cor_gene_gtf = os.path.join(corr_dir, rna_sm, "gene", f"{rna_sm}_out.gtf")
    if os.path.exists(cor_gtf) and os.path.getsize(cor_gtf) > 0:
        corr_lines.append(f"{label} {cor_gtf} {grp} {patient}")
        corr_gene_lines.append(f"{label} {cor_gene_gtf}")
    else:
        n_corrected_skipped += 1
        print(f"INFO: no RNase R-corrected result for {circ_sm} (expected {cor_gtf}) -- excluded from corrected/primary matrix", file=sys.stderr)

def write(path, lines):
    with open(path, "w") as f:
        if lines:
            f.write("\n".join(lines) + "\n")

write(os.path.join(outdir, "uncorrected", "sample.lst"), uncorr_lines)
write(os.path.join(outdir, "uncorrected", "sample_gene.lst"), uncorr_gene_lines)
write(os.path.join(outdir, "corrected", "sample.lst"), corr_lines)
write(os.path.join(outdir, "corrected", "sample_gene.lst"), corr_gene_lines)

n_t_unc = sum(1 for l in uncorr_lines if l.split(" ")[2] == "T")
n_c_unc = len(uncorr_lines) - n_t_unc
n_t_cor = sum(1 for l in corr_lines if l.split(" ")[2] == "T")
n_c_cor = len(corr_lines) - n_t_cor
print(f"uncorrected: {n_t_unc} tumor / {n_c_unc} normal ({len(uncorr_lines)} total)")
print(f"corrected  : {n_t_cor} tumor / {n_c_cor} normal ({len(corr_lines)} total); "
      f"{n_corrected_skipped} circRNA-seq samples had no matched RNase R correction and were excluded")
PY

###########
# prep_CIRIquant + prepDE.py per mode
###########
for MODE in corrected uncorrected; do
    MODE_DIR="${OUTDIR}/${MODE}"

    if [[ ! -s "${MODE_DIR}/sample.lst" ]]; then
        echo "${MODE}: sample.lst is empty -- nothing to quantify, skipping" >&2
        continue
    fi

    echo "=== ${MODE}: prep_CIRIquant ($(wc -l < "${MODE_DIR}/sample.lst") samples) ==="
    prep_CIRIquant \
        -i "${MODE_DIR}/sample.lst" \
        --lib "${MODE_DIR}/library_info.csv" \
        --circ "${MODE_DIR}/circRNA_info.csv" \
        --bsj "${MODE_DIR}/circRNA_bsj.csv" \
        --ratio "${MODE_DIR}/circRNA_ratio.csv"

    echo "=== ${MODE}: prepDE.py (gene-level counts, for CIRI_DE_replicate --gene) ==="
    ( cd "${MODE_DIR}" && prepDE.py -i sample_gene.lst )
    [[ -s "${MODE_DIR}/gene_count_matrix.csv" ]] || echo "WARNING: ${MODE_DIR}/gene_count_matrix.csv not produced by prepDE.py" >&2

    echo "${MODE}: done -> ${MODE_DIR}"
done
