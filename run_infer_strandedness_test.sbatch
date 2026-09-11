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
# Determine the strandedness of the RNase R circRNA libraries so the correct
# CIRIquant --library-type (and nf-core/circrna `strandedness` column) can be
# set. Strandedness is fixed by library prep, so this only needs to run once
# on one representative sample (SAMPLE01_tumor).
#
# Method: align a small random subset of reads with HISAT2 (splice-aware) and
# run RSeQC infer_experiment.py against a GENCODE gene model. Read the result
# with the "HOW TO READ THE RESULT" block printed at the end of the job.
#
# The HISAT2 index and the BED12 built here are reused later: the index is a
# required CIRIquant config input (hisat_index); the BED12 is generic RSeQC QC.
#

###########
# Software setup
###########
module purge
module load miniforge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate circrna   
# RSeQC is not py3.14-compatible in circrna, so it lives in its own env and
# step 4 calls it via `conda run`. 
# see .yaml fies for env configurations

echo ____________________________________
which hisat2 ; hisat2 --version | head -1
echo ____________________________________
conda run -n rseqc which infer_experiment.py #makes sure rseqc env works, calling env in line
echo ____________________________________

###########
# Variables
###########
REF_DIR="/projects/$USER/HNSCC_analysis/reference_genomes"
FASTA="${REF_DIR}/hg38/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode/gencode.v50.annotation.gtf"
HISAT_INDEX="${REF_DIR}/hg38/GRCh38.primary_assembly.genome_hisat2"
BED12="${REF_DIR}/gencode/gencode.v50.annotation.bed12"

#FASTQ_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/trimmed" ##circRNAseq
FASTQ_DIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/total_RNAseq/trimmed" ##total RNA-seq
SAMPLE="SAMPLE01_tumor"
R1="${FASTQ_DIR}/${SAMPLE}*_R1.fq.gz"
R2="${FASTQ_DIR}/${SAMPLE}*_R2.fq.gz"

WORKDIR="/scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/strand_check"
N_READS=2000000

mkdir -p "${WORKDIR}"

###########
# 1. GENCODE gene model -> BED12  (built once, then reused)
###########
if [ ! -s "${BED12}" ]; then
  echo ">> [1/4] Converting GTF to BED12"
  gtfToGenePred -genePredExt "${GTF}" "${WORKDIR}/tmp.genePred"
  genePredToBed "${WORKDIR}/tmp.genePred" "${BED12}"
  rm -f "${WORKDIR}/tmp.genePred"
else
  echo ">> [1/4] BED12 already present, skipping: ${BED12}"
fi

###########
# 2. HISAT2 genome index  (built once, ~30-60 min; also a CIRIquant input)
###########
if [ ! -s "${HISAT_INDEX}.1.ht2" ]; then
  echo ">> [2/4] Building HISAT2 index"
  hisat2-build -p 8 "${FASTA}" "${HISAT_INDEX}"
else
  echo ">> [2/4] HISAT2 index already present, skipping: ${HISAT_INDEX}"
fi

###########
# 3. Subsample reads and align (splice-aware, default unstranded mode)
###########
echo ">> [3/4] Subsampling ${N_READS} read pairs and aligning"
seqtk sample -s100 "${R1}" ${N_READS} | gzip > "${WORKDIR}/sub_R1.fq.gz"
seqtk sample -s100 "${R2}" ${N_READS} | gzip > "${WORKDIR}/sub_R2.fq.gz"

hisat2 -p 8 -x "${HISAT_INDEX}" \
       -1 "${WORKDIR}/sub_R1.fq.gz" -2 "${WORKDIR}/sub_R2.fq.gz" \
       --no-unal \
  | samtools sort -@ 8 -o "${WORKDIR}/${SAMPLE}.sub.bam" -
samtools index "${WORKDIR}/${SAMPLE}.sub.bam"

###########
# 4. Infer strandedness
###########
echo ">> [4/4] Running RSeQC infer_experiment.py (from the dedicated 'rseqc' env)"
conda run --no-capture-output -n rseqc infer_experiment.py \
    -i "${WORKDIR}/${SAMPLE}.sub.bam" \
    -r "${BED12}" \
    -s 400000 \
  | tee "${WORKDIR}/${SAMPLE}.infer_strandedness.txt"

###########
# HOW TO READ THE RESULT
###########
cat <<'EOF'

============================================================
 infer_experiment.py output -> CIRIquant --library-type
============================================================
 The report gives two informative fractions:

   "1++,1--,2+-,2-+"   = forward / sense     (FR, fr-secondstrand)
   "1+-,1-+,2++,2--"   = reverse / antisense (RF, fr-firststrand)

 Choose from whichever dominates:

   both ~0.5 .................. unstranded  -> --library-type 0   (strandedness: unstranded)
   forward fraction >= ~0.80 .. forward     -> --library-type 1   (strandedness: forward)
   reverse fraction >= ~0.80 .. reverse     -> --library-type 2   (strandedness: reverse)

 Update samplesheet `strandedness` and the CIRIquant command to match.
============================================================
EOF
