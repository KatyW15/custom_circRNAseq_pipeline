# custom_circRNAseq_pipeline

A hand-built circRNA-seq analysis pipeline (**CIRI2 + CIRIquant**) for HNSCC
tumor/normal samples: RNase R–treated, ribo-depleted total RNA, paired-end.
It is a set of SLURM (`sbatch`) scripts for a HPC cluster (developed on CU
Boulder's Alpine), plus a written protocol. A parallel
[nf-core/circrna](https://nf-co.re/circrna) run is used to cross-check results.

> **Note:** these scripts are shared as a reference/template. Cluster paths have
> been replaced with `$USER` and emails with a placeholder — see
> [Before you run](#before-you-run). No sequencing data or patient information is
> included.

## Approach

Two accepted strategies exist for circRNA analysis:

1. CIRIquant end-to-end (CIRI2 runs inside it), single tool.
2. **≥2-tool consensus detection, then CIRIquant for quantification.** ← this pipeline

Option 2 is more conservative and matches the nf-core/circrna cross-check
(`--tools circexplorer2,ciri,find_circ --min_tools 2 --bsj_reads 2`). CIRCexplorer2
and find_circ callers are outlined but not yet implemented here.

## Pipeline steps

| Step | Purpose | Script |
|------|---------|--------|
| 0 | Strandedness check (once, one sample): HISAT2 + RSeQC `infer_experiment.py`. Result for this dataset: **reverse → `--library-type 2`** | [run_infer_strandedness_test.sbatch](run_infer_strandedness_test.sbatch) |
| 1 | Raw-read QC: FastQC (+ optional MultiQC) | [run_fastqc_circRNA.sbatch](run_fastqc_circRNA.sbatch) |
| 2 | Adapter + light quality trim: fastp, paired-end, min length 50 bp | [run_fastp_circRNA.sbatch](run_fastp_circRNA.sbatch) |
| 3 | Map with BWA-MEM (`-T 19`, **unsorted SAM** — CIRI2 needs name order + SA tags) | [run_bwa_index.sbatch](run_bwa_index.sbatch), [run_bwa_circRNA.sbatch](run_bwa_circRNA.sbatch) |
| 4 | BSJ detection with CIRI2 from the BWA SAM; first-pass filter `#junction_reads >= 2`; emit BED of loci | [run_ciri2_circRNA.sbatch](run_ciri2_circRNA.sbatch) |
| 4b | Locate `CIRI2.pl` before submitting step 4 | [check_ciri2.sh](check_ciri2.sh) |
| 5 | Quantification: CIRIquant (HISAT2 + StringTie linear model, pseudo-circular realignment) | *outlined* |
| 6 | RNase R effect correction (needs a matched RNase R(−) mock library per sample) | *outlined* |
| 7 | Downstream: `prep_CIRIquant` matrices, `CIRI_DE_replicate` (edgeR) tumor vs normal | *outlined* |

Full rationale for each step (parameter choices, RNase R caveats, coordinate
conventions) is in
[circRNA_analysis_pipeline_outline.txt](circRNA_analysis_pipeline_outline.txt).

## References / build

- Genome: `GRCh38.primary_assembly.genome.fa` (GENCODE-distributed)
- Annotation: `gencode.v50.annotation.gtf`
- Indexes built once and reused: `bwa index -a bwtsw`, `hisat2-build`,
  `samtools faidx`, GTF→BED12 (`gtfToGenePred | genePredToBed`, RSeQC QC only)

## Environments

Two conda envs (RSeQC is not compatible with the Python in the main env, so it is
called via `conda run -n rseqc ...`):

- **`circrna`**: fastqc, fastp, bwa, hisat2, stringtie, samtools, ciriquant,
  seqtk, ucsc-gtftogenepred, ucsc-genepredtobed
- **`rseqc`**: python=3.11, rseqc>=5

## Before you run

These scripts are **not runnable as-is** — edit them for your environment:

1. **Paths:** every `/scratch/alpine/$USER/...` and `/projects/$USER/...` path is
   a placeholder. SLURM does **not** expand shell variables in `#SBATCH -o/-e`
   directives, so edit those lines to a literal path.
2. **Email:** replace `--mail-user=YOUR_EMAIL@example.com`.
3. **Sample sheet:** the scripts read a `master_samplesheet.csv` with columns
   `patient_id, batch_id, tissue, hpv_status, modality, sample_id, fastq_1,
   fastq_2, survival_time_days, vital_status, notes` and filter to
   `modality == "circrnaseq"`. Adjust the `awk` filters and `--array` range for
   your data.
4. **Modules:** scripts use `module load miniforge`; change to your cluster's
   module system.
5. **Representative sample:** step 0 hard-codes `SAMPLE01_tumor` — set it to one
   of your samples.

## License

[MIT](LICENSE)
