#!/bin/bash

# Phase 1
# Setting up sbatch to avoid dependency error
alias sbatch='sbatch --export=NONE'
SRC_DIR="/path/to/fastq_dir/"
total_files=`ls "$SRC_DIR"/*.fastq.gz | wc -l`

tmp1=$(sbatch --array=1-$total_files slurm_script.sh)
jid1=`echo ${tmp1##* }`

# Phase 2
tmp2=$(sbatch --array=1-$total_files --dependency=afterok:$jid1 Alpine_BLAST.slurm)

jid2=`echo ${tmp2##* }`




