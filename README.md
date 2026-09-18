# custom_circRNAseq_pipeline

A hand-built circRNA-seq analysis pipeline (**CIRI2 + CIRIquant**) for HNSCC
tumor/normal samples: RNase R–treated, ribo-depleted total RNA, paired-end.
It is a set of SLURM (`sbatch`) scripts for a HPC cluster (developed on CU
Boulder's Alpine), plus a written protocol. A parallel
[nf-core/circrna](https://nf-co.re/circrna) run is used to cross-check results.

> **Note:** these scripts are shared as a reference/template. Cluster paths have
> been replaced with `$USER` and emails with `YOUR_EMAIL@example.com` — see
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
| 4c | Consensus BSJ set: merge per-tool BEDs, tag `n_tools`, keep ≥`MIN_TOOLS`; build cohort union BED. **CIRI2-only for now** (`MIN_TOOLS=1`); set `MIN_TOOLS=2` when CIRCexplorer2 / find_circ are added | [build_bsj_consensus.sbatch](build_bsj_consensus.sbatch) |
| 5 | Quantification: CIRIquant (HISAT2 + StringTie linear model, pseudo-circular realignment), prior-guided (`--bed` consensus union), `--library-type 2` | [run_ciriquant_circRNA_py2.sbatch](run_ciriquant_circRNA_py2.sbatch), [ciriquant_config.yaml](ciriquant_config.yaml) |
| 6a | **RNase R effect correction** (total RNA-seq branch): CIRIquant on the matched ribo-depleted total RNA-seq with `--RNaseR` = the step-5 circRNA-seq gtf. Keeps corrected (`ciriquant_py2/`) **and** uncorrected (`circRNAseq/ciriquant_noRNaseRcorrection_py2/`) | [run_ciriquant_RNaseR_corrected_py2.sbatch](run_ciriquant_RNaseR_corrected_py2.sbatch) |
| 6b | **Validation** (total RNA-seq branch): trim → BWA-MEM → CIRI2 on the total RNA-seq (+ FastQC and `samtools flagstat` QC), then intersect with the circRNA-seq consensus → per-circRNA "seen in non-enriched data" confidence flag | [run_fastp_totalRNAseq.sbatch](run_fastp_totalRNAseq.sbatch), [run_bwa_totalRNAseq.sbatch](run_bwa_totalRNAseq.sbatch), [run_ciri2_totalRNAseq.sbatch](run_ciri2_totalRNAseq.sbatch), [run_fastqc_totalRNAseq.sbatch](run_fastqc_totalRNAseq.sbatch), [run_flagstat_totalRNAseq.sbatch](run_flagstat_totalRNAseq.sbatch), [build_bsj_validation.sbatch](build_bsj_validation.sbatch) |
| 6b-chain | Submit 6b pre-processing as a dependency chain: **trim → BWA → CIRI2** (array stages, `aftercorr` — one bad sample drops out instead of blocking all), plus **FastQC** (`afterany` on CIRI2) and **flagstat** (`afterany` on BWA — mapping-rate summary). `bash job_inchain_totalRNAseq.sh` on a login node — submits the five jobs, exits in seconds, then you can log off | [job_inchain_totalRNAseq.sh](job_inchain_totalRNAseq.sh) |
| 7a | Downstream: `prep_CIRIquant` + `prepDE.py` matrices, corrected (primary) + uncorrected (sensitivity) | [build_ciriquant_matrix.sbatch](build_ciriquant_matrix.sbatch) |
| 7b | Downstream: `CIRI_DE_replicate` (edgeR) tumor vs normal, paired by `patient_id`. Annotating circRNAs with the 6b validation flag is still *outlined only* | [run_ciri_de.sbatch](run_ciri_de.sbatch) |

Total RNA-seq per-sample outputs are written under `.../total_RNAseq/` sub-dirs
(`trimmed/total_RNAseq/`, `QC/total_RNAseq/`, `bwa_mapped/total_RNAseq/`,
`ciri2/total_RNAseq/`, `ciriquant/total_RNAseq/`) so they stay separate from the
circRNA-seq outputs.

Full rationale for each step (parameter choices, RNase R caveats, coordinate
conventions) is in
[circRNA_analysis_pipeline_outline.txt](circRNA_analysis_pipeline_outline.txt).

## References / build

- Genome: `GRCh38.primary_assembly.genome.fa` (GENCODE-distributed)
- Annotation: `gencode.v50.annotation.gtf`
- Indexes built once and reused: `bwa index -a bwtsw`, `hisat2-build`,
  `samtools faidx`, GTF→BED12 (`gtfToGenePred | genePredToBed`, RSeQC QC only)

## Environments

Six conda envs total; full package lists, version pins, create commands, and
which script uses which are in
[software_requirements.md](software_requirements.md). Summary:

- **`circrna`**: fastqc, fastp, bwa, hisat2, stringtie, samtools, ciri2, seqtk,
  ucsc-gtftogenepred, ucsc-genepredtobed. General-purpose env for everything
  except invoking CIRIquant's own executables.
- **`CIRIquant_github_env`**: the official CIRIquant v1.1.3 release, installed
  per the developers' own documented method (CIRI-cookbook pinned
  `environment.yml` -- python=2.7.15, pysam==0.15.2, etc.), on its natively
  supported Python 2 interpreter. Replaces an earlier `ciriquant` env
  (bioconda build, patched for Python 3.11) that kept surfacing new upstream
  Python-2-only bugs one at a time across multi-hour runs; py2 has been
  running cleanly and that env/its patches are retired.
- **`ciriDE`**: R + edgeR/limma/statmod/optparse for step 7b's
  `CIRI_DE_replicate` -- kept in its **own** env rather than installed into
  `ciriquant`, because `CIRI_DE_replicate` just shells out to whatever
  `Rscript` is on `PATH` (no rpy2), and CIRIquant's docs-recommended
  `r-base=3.6` pin fails to solve against `ciriquant`'s already-modern
  pysam/openssl stack. Channel order matters here too -- conda-forge must come
  first, or bioconda's r-base build crashes at runtime needing an obsolete
  `libgfortran.so.3`: `conda create -n ciriDE -c conda-forge -c bioconda
  r-base bioconductor-edger bioconductor-limma r-statmod r-optparse`.
- **`rseqc`**: python=3.11, rseqc>=5 -- not compatible with the other envs'
  Python, so it's called via `conda run -n rseqc ...` rather than activated.
- **`sortmerna`**: sortmerna=4.3.7, called via `conda run -n sortmerna ...`.
- **`dupradar`**: bioconductor-dupradar, called via `conda run -n dupradar ...`.

## Before you run

These scripts are **not runnable as-is** — edit them for your environment:

1. **Paths:** every `/scratch/alpine/$USER/...` and `/projects/$USER/...` path is
   a placeholder. SLURM does **not** expand shell variables in `#SBATCH -o/-e`
   directives, so edit those lines to a literal path.
2. **Email:** replace `--mail-user=YOUR_EMAIL@example.com`.
3. **Sample sheet:** the scripts read a `master_samplesheet.csv` with columns
   `patient_id, batch_id, tissue, hpv_status, modality, sample_id, fastq_1,
   fastq_2, survival_time_days, vital_status, notes` and filter to
   `modality == "circrnaseq"` (or `"rnaseq"` for the `*_totalRNAseq.sbatch`
   scripts). Adjust the `awk` filters and `--array` range for your data.
4. **Modules:** scripts use `module load miniforge`; change to your cluster's
   module system.
5. **Representative sample:** step 0 hard-codes `SAMPLE01_tumor` — set it to one
   of your samples.
6. **Total RNA-seq branch (`*_totalRNAseq.sbatch`, `job_inchain_totalRNAseq.sh`,
   `build_bsj_validation.sbatch`):** array ranges are hard-coded to `1-65`
   (rnaseq samples with real paired FASTQ paths); recount with the
   `awk … | wc -l` in each script header if the sample sheet changes.
   `run_ciriquant_RNaseR_corrected_py2.sbatch` defaults `LIBRARY_TYPE=2` for
   the rnaseq libraries (a different prep/batch from the circRNA-seq set) —
   confirmed reverse-stranded for this dataset (both batches checked; see
   step 0).
   `job_inchain_totalRNAseq.sh` submits trim → BWA → CIRI2 (`aftercorr`) plus
   FastQC (`afterany` CIRI2) and flagstat (`afterany` BWA); CIRIquant (6a) and
   the strandedness check are run separately.
   Run it with `bash job_inchain_totalRNAseq.sh` on a login node — it submits the
   five jobs (child submissions use `--export=NONE`) and returns in seconds; the
   queued jobs then outlive your session.

## License

[MIT](LICENSE)
