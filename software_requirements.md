# Software requirements

Every conda env this pipeline touches, what's in each, how to create it, and
which scripts use it. Six envs total. All are activated after
`module load miniforge` (see [Cluster modules](#cluster-modules-not-conda)
below for the one thing loaded outside conda).

## `circrna`

General-purpose env for everything except actually invoking CIRIquant's own
executables (that's `CIRIquant_github_env`, below) or RSeQC/SortMeRNA/dupRadar
(each in their own env for python/R compatibility reasons).

| Package | Notes |
|---|---|
| fastqc | raw-read QC |
| fastp | adapter/quality trim |
| bwa | BWA-MEM mapping (circRNA-seq + total RNA-seq), and CIRI2's required aligner |
| hisat2 | strandedness check, dupRadar's alignment step |
| stringtie | bundled for consistency; not directly invoked outside `ciriquant` env runs |
| samtools | BAM/SAM handling; several scripts instead use the cluster's `module load samtools` (see below) rather than this copy -- check each script's header |
| ciri2 | BSJ detection from BWA SAMs. Not on any default channel -- install with `conda install -n circrna -c bioconda ciri2`; locate it first with [check_ciri2.sh](check_ciri2.sh) |
| seqtk | read subsampling (strandedness check, rRNA screen) |
| ucsc-gtftogenepred, ucsc-genepredtobed | GTF -> BED12, for RSeQC input |
| multiqc | optional -- QC scripts check `command -v multiqc` and skip aggregation if absent |

**Used by:** [run_fastqc_circRNA.sbatch](run_fastqc_circRNA.sbatch),
[run_fastqc_totalRNAseq.sbatch](run_fastqc_totalRNAseq.sbatch),
[run_fastp_circRNA.sbatch](run_fastp_circRNA.sbatch),
[run_fastp_totalRNAseq.sbatch](run_fastp_totalRNAseq.sbatch),
[run_bwa_index.sbatch](run_bwa_index.sbatch),
[run_bwa_circRNA.sbatch](run_bwa_circRNA.sbatch),
[run_bwa_totalRNAseq.sbatch](run_bwa_totalRNAseq.sbatch),
[run_ciri2_circRNA.sbatch](run_ciri2_circRNA.sbatch),
[run_ciri2_totalRNAseq.sbatch](run_ciri2_totalRNAseq.sbatch),
[check_ciri2.sh](check_ciri2.sh),
[run_infer_strandedness_test.sbatch](run_infer_strandedness_test.sbatch) (also
calls into `rseqc`), [check_rrna_totalRNAseq.sbatch](check_rrna_totalRNAseq.sbatch)
(also calls into `sortmerna`),
[run_dupradar_totalRNAseq.sbatch](run_dupradar_totalRNAseq.sbatch) (also calls
into `dupradar`).

## `CIRIquant_github_env`

**Separate from `circrna`.** The official CIRIquant v1.1.3 release, installed
per the developers' own documented method -- CIRI-cookbook's pinned
`environment.yml` (Method 2), which installs `CIRIquant` from PyPI on its
natively supported Python 2 interpreter, alongside the exact tool versions it
was built/tested against:

| Package | Version pin | Notes |
|---|---|---|
| python | 2.7.15 | native -- CIRIquant v1.1.3 explicitly only supports Python 2 |
| bwa | 0.7.17 | |
| hisat2 | 2.2.0 | |
| stringtie | 2.1.1 | also provides `prepDE.py` |
| samtools | >=1.10 | |
| pysam | 0.15.2 | pip-installed, pinned |
| numpy | 1.16.4 | pip-installed, pinned |
| scipy | 1.2.2 | pip-installed, pinned |
| scikit-learn | 0.20.3 | pip-installed, pinned |
| PyYAML | 5.4 | pip-installed, pinned |
| CIRIquant | >=1.1.2 | pip-installed from PyPI; provides `CIRIquant`, `prep_CIRIquant`, `CIRI_DE_replicate` |

Built from [ciriquant_py2_environment.yml](ciriquant_py2_environment.yml):
`mamba env create -f ciriquant_py2_environment.yml`. The resulting env landed
outside the standard `envs_dirs` search path the first time (built via the
CIRI-cookbook's alternate "packed env" tarball route before switching to the
environment.yml route) and had to be registered/moved in; see git history /
conversation notes if setting this up fresh trips the same
`EnvironmentLocationNotFound` surprise.

This env's `activate.d` hooks (from its compiler-toolchain packages, e.g.
`binutils_linux-64`) reference variables like `$ADDR2LINE` without a default
-- harmless normally, but fatal ("unbound variable") under `set -u`. Every
script below wraps its `conda activate CIRIquant_github_env` call in
`set +u` / `set -u` for exactly this reason.

`CIRI_DE_replicate` (step 7b) needs R + edgeR at runtime, but **not inside
this env** -- see the separate [`ciriDE`](#ciride) env below. Reading
CIRIquant's source (`CIRIquant/replicate.py`) shows `CIRI_DE_replicate` just
does `subprocess.call('Rscript .../CIRI_DE.R ...')`: it picks up whatever
`Rscript` is first on `$PATH` at call time, no rpy2/embedding involved, so R
doesn't need to coexist in the same env.

**Used by:** [run_ciriquant_circRNA_py2.sbatch](run_ciriquant_circRNA_py2.sbatch),
[run_ciriquant_RNaseR_corrected_py2.sbatch](run_ciriquant_RNaseR_corrected_py2.sbatch),
[build_ciriquant_matrix.sbatch](build_ciriquant_matrix.sbatch),
[run_ciri_de.sbatch](run_ciri_de.sbatch).

### Retired: `ciriquant` (bioconda, Python 3.11, patched)

An earlier env -- bioconda's CIRIquant build, pinned to Python 3.9 to dodge a
`distutils`-removed-in-3.12 crash, then further patched (see git history:
`patch_ciriquant_all.sh` and its predecessors) for five separate upstream
Python-2-only bugs found one at a time, each after a multi-hour run reached a
new code path (`cmp()`/`LooseVersion` version comparisons, `xrange`/
`.iteritems()`/integer-division/`sorted(cmp=...)` in `circ.py`). After the
fifth fix the patched build was *still* producing bad output on some samples,
while the official Python 2 install (above) ran the full cohort cleanly on
the first attempt with zero patches needed -- so this env, its patch scripts,
and everything under `*_totalRNAseq_independent*` (a third, now-dropped
quantification analysis that only existed for this env) were removed. Kept
here as a note in case anyone wonders why a `ciriquant_config.yaml` template
or old script names turn up in git history referencing a `ciriquant` env that
no longer exists in this repo.

## `rseqc`

| Package | Version pin |
|---|---|
| python | 3.11 |
| rseqc | >=5 |

Not compatible with `circrna`'s Python, so it's called from inside
`run_infer_strandedness_test.sbatch` via `conda run -n rseqc
infer_experiment.py ...` rather than activated directly.

**Used by:** [run_infer_strandedness_test.sbatch](run_infer_strandedness_test.sbatch).

## `sortmerna`

| Package | Version pin | Notes |
|---|---|---|
| sortmerna | 4.3.7 | pinned -- must match the `smr_v4.3_*` SILVA+Rfam database index format the script downloads |

Create with:
```
conda create -n sortmerna -c bioconda -c conda-forge sortmerna=4.3.7
```
Called via `conda run -n sortmerna sortmerna ...`; the calling script only
activates `circrna` (for `seqtk` subsampling) and checks this env exists with
a preflight, failing fast with the create command above if not.

**Used by:** [check_rrna_totalRNAseq.sbatch](check_rrna_totalRNAseq.sbatch).

## `ciriDE`

R + edgeR for `CIRI_DE_replicate` (step 7b) -- kept separate from `ciriquant`
for the solver-conflict reason explained above, not for a Python-version
reason like `rseqc`/`dupradar`.

| Package | Notes |
|---|---|
| r-base | unpinned -- let conda resolve a modern, mutually-compatible version rather than reusing CIRIquant's outdated `r-base=3.6` docs pin, which conflicts when solved against this pipeline's other already-installed packages |
| bioconductor-edger | required by `CIRI_DE_replicate` |
| bioconductor-limma | required by `CIRI_DE_replicate` |
| r-statmod | required by `CIRI_DE_replicate` |
| r-optparse | required by `CIRI_DE_replicate` |

Create with (channel order matters -- conda-forge must come first, or the
solver pulls `r-base` from bioconda's old GCC4-linked build, which crashes at
runtime with `error while loading shared libraries: libgfortran.so.3: cannot
open shared object file` -- that ABI doesn't exist on modern systems.
`bioconductor-edger`/`bioconductor-limma` only exist on bioconda, so they
still resolve from there regardless of order; only `r-base`'s toolchain needs
conda-forge to win):
```
conda create -n ciriDE -c conda-forge -c bioconda \
    r-base bioconductor-edger bioconductor-limma r-statmod r-optparse
```
Unlike `sortmerna`/`dupradar` (called via `conda run -n <env> ...` end to
end), `CIRI_DE_replicate` itself must run as the `CIRIquant_github_env`
Python entry point -- it's the thing that shells out to `Rscript` internally.
So [run_ciri_de.sbatch](run_ciri_de.sbatch) instead prepends
`$(conda info --base)/envs/ciriDE/bin` to `PATH` after activating
`CIRIquant_github_env`, so `CIRI_DE_replicate`'s internal `Rscript` call
resolves into `ciriDE` without switching envs. The preflight check still uses
`conda run -n ciriDE ...`, matching the other envs' pattern.

**Used by:** [run_ciri_de.sbatch](run_ciri_de.sbatch).

## `dupradar`

| Package | Notes |
|---|---|
| bioconductor-dupradar | Bioconductor R package; brings its own R + featureCounts dependency |

Create with:
```
conda create -n dupradar -c bioconda -c conda-forge bioconductor-dupradar
```
Called via `conda run -n dupradar Rscript ...`; the calling script only
activates `circrna` (for `hisat2`/`samtools`) and checks this env exists with
a preflight, failing fast with the create command above if not.

**Used by:** [run_dupradar_totalRNAseq.sbatch](run_dupradar_totalRNAseq.sbatch).

## Cluster modules (not conda)

| Module | Notes |
|---|---|
| miniforge | provides `conda`/`mamba` itself; every script starts with `module load miniforge` before any `conda activate` |
| samtools | loaded directly via `module load samtools` (bypassing any conda env's copy) in [run_bwa_circRNA.sbatch](run_bwa_circRNA.sbatch), [run_bwa_totalRNAseq.sbatch](run_bwa_totalRNAseq.sbatch), [run_flagstat_circRNA.sbatch](run_flagstat_circRNA.sbatch), [run_flagstat_totalRNAseq.sbatch](run_flagstat_totalRNAseq.sbatch) |

Substitute your cluster's own module names if not on CU Boulder Alpine (see
README "Before you run" #4).

## Not software, but required to run anything

Reference genome/annotation (GRCh38 primary assembly + GENCODE v50 GTF) and
their derived indexes (BWA, HISAT2, `.fai`) -- see README
["References / build"](README.md#references--build). None of the above
package installs substitute for having these built first.
