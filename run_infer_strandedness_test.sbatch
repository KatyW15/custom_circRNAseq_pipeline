#!/bin/bash
#SBATCH -J infer_strand_TEST
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/infer_strand_TEST_%j.out  # Output file with the job ID
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/infer_strand_TEST_%j.err  # Error file with the job ID
#SBATCH -t 06:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 8  # ask for number of nodes/cores
#SBATCH --mem=32GB  # Specify memory allocation


set -euo pipefail

#
# Determine the strandedness of the total RNA-seq libraries, checking ONE
# SAMPLE FROM EACH BATCH (pilot, batch1 -- master_samplesheet.csv batch_id) to
# confirm both batches share the same strandedness before applying one
# CIRIquant --library-type (and nf-core/circrna `strandedness` column) to all
# samples. If the two batches disagree, --library-type must be set per batch
# instead of globally.
#
# Method: for each batch, auto-pick the first rnaseq sample_id with real
# trimmed FASTQs, align a small random subset of its reads with HISAT2
# (splice-aware), and run RSeQC infer_experiment.py against a GENCODE gene
# model. Compare the two reports with the "HOW TO READ THE RESULT" block
# printed at the end.
#
# The HISAT2 index and BED12 built here are shared/reused across both batches:
# the index is a required CIRIquant config input (hisat_index); the BED12 is
# generic RSeQC QC.
#

###########
# Software setup
###########
module purge
module load miniforge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate circrna
# RSeQC is not py3.14-compatible in circrna, so it lives in its own env and
# the infer_experiment.py step calls it via `conda run`.
# see .yaml fies for env configurations

echo ____________________________________
which hisat2 ; hisat2 --version | head -1
echo ____________________________________
conda run -n rseqc which infer_experiment.py #makes sure rseqc env works, calling env in line
echo ____________________________________

###########
# Variables
###########
SAMPLESHEET="/projects/$USER/HNSCC_analysis/scripts/master_samplesheet.csv"

REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
FASTA="${REF_DIR}/hg38/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode/gencode.v50.annotation.gtf"
HISAT_INDEX="${REF_DIR}/hg38/GRCh38.primary_assembly.genome_hisat2"
BED12="${REF_DIR}/gencode/gencode.v50.annotation.bed12"

#FASTQ_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/circRNAseq/trimmed" ##circRNAseq
FASTQ_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/total_RNAseq/trimmed" ##total RNA-seq

WORKDIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/strand_check"
N_READS=2000000

mkdir -p "${WORKDIR}"

###########
# Pick one rnaseq sample per batch (first real-FASTQ row for that batch)
###########
# master_samplesheet.csv columns: 1 patient_id 2 batch_id 3 tissue 4 hpv_status
# 5 modality 6 sample_id 7 fastq_1 8 fastq_2
PILOT_SAMPLE=$(awk -F',' 'NR>1 && $5=="rnaseq" && $2=="pilot"  && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6; exit}' "${SAMPLESHEET}")
BATCH1_SAMPLE=$(awk -F',' 'NR>1 && $5=="rnaseq" && $2=="batch1" && $7 ~ /\.f(ast)?q\.gz$/ && $8 ~ /\.f(ast)?q\.gz$/ {print $6; exit}' "${SAMPLESHEET}")

[[ -n "${PILOT_SAMPLE}"  ]] || { echo "ERROR: no rnaseq sample found for batch_id=pilot in ${SAMPLESHEET}"  >&2; exit 1; }
[[ -n "${BATCH1_SAMPLE}" ]] || { echo "ERROR: no rnaseq sample found for batch_id=batch1 in ${SAMPLESHEET}" >&2; exit 1; }

declare -A BATCH_SAMPLE=( [pilot]="${PILOT_SAMPLE}" [batch1]="${BATCH1_SAMPLE}" )
echo "pilot  sample: ${PILOT_SAMPLE}"
echo "batch1 sample: ${BATCH1_SAMPLE}"

###########
# 1. GENCODE gene model -> BED12  (built once, then reused)
###########
if [ ! -s "${BED12}" ]; then
  echo ">> [1/2] Converting GTF to BED12"
  gtfToGenePred -genePredExt "${GTF}" "${WORKDIR}/tmp.genePred"
  genePredToBed "${WORKDIR}/tmp.genePred" "${BED12}"
  rm -f "${WORKDIR}/tmp.genePred"
else
  echo ">> [1/2] BED12 already present, skipping: ${BED12}"
fi

###########
# 2. HISAT2 genome index  (built once, ~30-60 min; also a CIRIquant input)
###########
if [ ! -s "${HISAT_INDEX}.1.ht2" ]; then
  echo ">> [2/2] Building HISAT2 index"
  hisat2-build -p 8 "${FASTA}" "${HISAT_INDEX}"
else
  echo ">> [2/2] HISAT2 index already present, skipping: ${HISAT_INDEX}"
fi

###########
# 3. Per batch: subsample reads, align (splice-aware), infer strandedness
###########
for BATCH in pilot batch1; do
  SAMPLE="${BATCH_SAMPLE[${BATCH}]}"
  R1="${FASTQ_DIR}/${SAMPLE}_R1.fq.gz"
  R2="${FASTQ_DIR}/${SAMPLE}_R2.fq.gz"
  [[ -s "${R1}" ]] || { echo "ERROR: missing trimmed R1 for ${BATCH} sample ${SAMPLE}: ${R1} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }
  [[ -s "${R2}" ]] || { echo "ERROR: missing trimmed R2 for ${BATCH} sample ${SAMPLE}: ${R2} (run run_fastp_totalRNAseq.sbatch)" >&2; exit 1; }

  echo ">> [${BATCH}] Subsampling ${N_READS} read pairs from ${SAMPLE} and aligning"
  seqtk sample -s100 "${R1}" ${N_READS} | gzip > "${WORKDIR}/${SAMPLE}.sub_R1.fq.gz"
  seqtk sample -s100 "${R2}" ${N_READS} | gzip > "${WORKDIR}/${SAMPLE}.sub_R2.fq.gz"

  hisat2 -p 8 -x "${HISAT_INDEX}" \
         -1 "${WORKDIR}/${SAMPLE}.sub_R1.fq.gz" -2 "${WORKDIR}/${SAMPLE}.sub_R2.fq.gz" \
         --no-unal \
    | samtools sort -@ 8 -o "${WORKDIR}/${SAMPLE}.sub.bam" -
  samtools index "${WORKDIR}/${SAMPLE}.sub.bam"

  echo ">> [${BATCH}] Running RSeQC infer_experiment.py on ${SAMPLE} (from the dedicated 'rseqc' env)"
  conda run --no-capture-output -n rseqc infer_experiment.py \
      -i "${WORKDIR}/${SAMPLE}.sub.bam" \
      -r "${BED12}" \
      -s 400000 \
    | tee "${WORKDIR}/${SAMPLE}.infer_strandedness.txt"
  echo
done

###########
# HOW TO READ THE RESULT
###########
cat <<EOF

============================================================
 infer_experiment.py output -> CIRIquant --library-type
============================================================
 Two reports were generated -- one per batch:
   pilot  (${PILOT_SAMPLE}) : ${WORKDIR}/${PILOT_SAMPLE}.infer_strandedness.txt
   batch1 (${BATCH1_SAMPLE}): ${WORKDIR}/${BATCH1_SAMPLE}.infer_strandedness.txt

 Each report gives two informative fractions:

   "1++,1--,2+-,2-+"   = forward / sense     (FR, fr-secondstrand)
   "1+-,1-+,2++,2--"   = reverse / antisense (RF, fr-firststrand)

 Choose from whichever dominates, per report:

   both ~0.5 .................. unstranded  -> --library-type 0   (strandedness: unstranded)
   forward fraction >= ~0.80 .. forward     -> --library-type 1   (strandedness: forward)
   reverse fraction >= ~0.80 .. reverse     -> --library-type 2   (strandedness: reverse)

 Then compare the two:
   SAME category in both reports  -> one --library-type is safe for every total
                                      RNA-seq sample regardless of batch.
   DIFFERENT categories           -> the two batches were prepped differently;
                                      set --library-type per batch instead (e.g.
                                      filter run_ciriquant_totalRNAseq*.sbatch's
                                      sample selection by batch_id, or add a
                                      per-batch LIBRARY_TYPE lookup).
============================================================
EOF
