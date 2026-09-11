#!/bin/bash
#SBATCH -J bwa_index
#SBATCH --partition acpu
#SBATCH --account=amc-general
#SBATCH -o /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bwa_index_%j.out
#SBATCH -e /scratch/alpine/$USER/HNSCC_analysis/my_circRNAseq_analysis/slurm_out/bwa_index_%j.err
#SBATCH -t 03:00:00    # D-HH:MM:SS  (human genome ~1-1.5 h)
#SBATCH --qos=cpu-normal
#SBATCH -n 1 -c 1      # bwa index is single-threaded
#SBATCH --mem=8GB

set -euo pipefail

###########
# Software setup
###########
module purge
module load miniforge

#activate my circrna env with all dependancies
conda activate circrna

echo ____________________________________
which bwa
bwa 2>&1 | grep -i version || true
echo ____________________________________

###########
# Variables
###########
REF="/projects/$USER/HNSCC_analysis/reference_genomes/hg38/GRCh38.primary_assembly.genome.fa"

###########
# Build the index (writes REF.amb/.ann/.bwt/.pac/.sa next to the FASTA)
###########
if [[ -s "${REF}.bwt" ]]; then
    echo "Index already present (${REF}.bwt) -- nothing to do."
    exit 0
fi

echo "Indexing ${REF}"
bwa index "${REF}"
echo "done"
