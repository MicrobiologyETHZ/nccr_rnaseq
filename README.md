# NCCR RNA-seq Pipeline

A Snakemake-based pipeline for processing bulk RNA-seq / metatranscriptomic
sequencing data, wrapped in a small command-line tool (`nccrRna`). Each analysis
is exposed as a subcommand that takes a single YAML config file.

The **`metat-bowtie`** workflow (read QC → Bowtie2 alignment → featureCounts) is
the primary, documented workflow and is described in full below. Other workflows
exist and share the same machinery but are still being stabilized — see
[Other workflows](#other-workflows-under-construction) at the bottom.

---

## Installation

1. Clone the repository and `cd` into it.

2. Create and activate the environment that runs the pipeline itself (Snakemake +
   the `nccrRna` CLI). Using `mamba` is recommended.

   ```shell
   mamba env create -f pipeline_environment.yaml
   conda activate pipeline
   ```

3. Install the `nccrRna` command into that environment:

   ```shell
   pip install -e .
   ```

   Confirm it worked:

   ```shell
   nccrRna --help
   ```

> The `pipeline` environment is all you need for `--dry` runs (planning the
> workflow without executing it). To actually run the tools, you also need the
> per-step tool environments — see [Running for real](#3-run).

---

## Test your install (dry run)

The repository ships a tiny example dataset (`test_data/`) and a ready-to-go
config ([`configs/metat-bowtie_config.yaml`](configs/metat-bowtie_config.yaml)).
From the repository root, with the `pipeline` environment active:

```shell
nccrRna metat-bowtie -c configs/metat-bowtie_config.yaml --dry
```

A `--dry` run plans the workflow without running any tools, so this needs only
the `pipeline` environment. You should see a job table (`bowtie_index`,
`bowtie_align`, `bowtie_featurecounts`, `qc`) ending with `This was a dry run`.

**No path editing required.** The example config refers to its paths with
`${NCCR_REPO}`, e.g.:

```yaml
samples:   ${NCCR_REPO}/test_data/samples.csv
refGenome: ${NCCR_REPO}/test_data/ref.fasta
```

`NCCR_REPO` defaults to the repository root, and the pipeline expands
`${ENV_VARS}` and `~` in config values (and in the sample-sheet FASTQ paths)
before running. So the bundled config works as-is from a fresh clone. To point
the same config at a checkout elsewhere, just override the variable:

```shell
NCCR_REPO=/path/to/your/clone nccrRna metat-bowtie -c configs/metat-bowtie_config.yaml --dry
```

> The test FASTQ/reference files are intentionally small. To run them through
> the actual tools (not just `--dry`), replace them with real data — see
> [Running for real](#3-run).

---

## Quick start: `metat-bowtie`

The workflow takes paired-end FASTQ files and a reference genome + annotation,
and produces a per-sample gene count table.

```
raw reads ──► QC/trimming (bbduk) ──► Bowtie2 align ──► featureCounts ──► counts
```

### 1. Prepare a sample sheet

A CSV with one row per FASTQ pair and the columns `sample,unit,fastq_1,fastq_2`:

| sample  | unit | fastq_1                       | fastq_2                       |
|---------|------|-------------------------------|-------------------------------|
| sample1 |      | /abs/path/sample1_1.fq.gz     | /abs/path/sample1_2.fq.gz     |
| sample2 |      | /abs/path/sample2_1.fq.gz     | /abs/path/sample2_2.fq.gz     |

Use absolute paths, or `${ENV_VAR}` / `~` references (these are expanded at run
time — see [Test your install](#test-your-install-dry-run)). You can write this
file by hand, or generate it from a directory of FASTQ files:

```shell
nccrRna samples -i <FASTQ_DIR> -o <OUTPUT.csv> -r1 _1.fq.gz -r2 _2.fq.gz -sn -sd _ -si 1
```

`-r1`/`-r2` are the forward/reverse filename suffixes. `-sn` cleans up the sample
name by splitting the filename on `-sd` (delimiter) and joining the first `-si`
fields — e.g. `sample1_L001_1.fq.gz` with `-sd _ -si 1` gives sample name
`sample1`. (Adapted from the
[nf-core/rnaseq](https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py)
helper.)

### 2. Write a config file

Copy [`configs/metat-bowtie_config.yaml`](configs/metat-bowtie_config.yaml) and
edit the paths to point at your data. Use absolute paths, or `${ENV_VAR}` / `~`
references (expanded at run time) — plain relative paths won't resolve, because
the pipeline runs Snakemake from inside the package directory.

```yaml
# ── Inputs ──
samples:   /abs/path/samples.csv     # CSV: sample,unit,fastq_1,fastq_2
dataDir:   /abs/path/raw             # base dir for reads (used for replicate merge)
outDir:    /abs/path/out             # all outputs go here

# ── Reference ──
refGenome: /abs/path/genome.fasta    # Bowtie2 index is built from this
refAnn:    /abs/path/annotation.gff  # featureCounts annotation (GTF/GFF)

# ── featureCounts ──
feature_type: CDS         # -t : feature type to count (exon, CDS, gene, …)
attribute:    gene_id     # -g : grouping attribute
strand:       0           # -s : 0 unstranded / 1 forward / 2 reverse

# ── Read QC (bbduk) ──
qc:        true           # false = use reads as-is, skip trimming
trimq:     14
minlen:    45
mink:      11
mapq:      20
adapters:  /abs/path/adapters.fa
phix:      /abs/path/phix174_ill.ref.fa.gz

# ── Misc ──
se:               false   # single-end reads?
merge_replicates: false   # merge multiple units per sample before alignment
```

### 3. Run

First, a **dry run** to check the plan (no tools required beyond the `pipeline`
env):

```shell
nccrRna metat-bowtie -c /abs/path/to/metat-bowtie_config.yaml --dry
```

This prints the jobs Snakemake would run. When it looks right, run for real.

To actually execute, Snakemake needs the per-step tool environments. Create them
once from the `envs/` directory (they are referenced by name inside the rules):

```shell
mamba env create -f nccr_rnaseq/envs/preprocessing.yaml   # bbduk (QC)
mamba env create -f nccr_rnaseq/envs/bowtie.yaml          # bowtie2, samtools, subread
```

Then run **locally**:

```shell
nccrRna metat-bowtie -c /abs/path/to/metat-bowtie_config.yaml --local --jobs 4
```

…or submit to a **SLURM cluster** (the default when `--local` is omitted):

```shell
nccrRna metat-bowtie -c /abs/path/to/metat-bowtie_config.yaml --jobs 20 -p <partition>
```

| Flag           | Effect                                                        |
|----------------|---------------------------------------------------------------|
| `--dry`        | Show the jobs/commands without running them                   |
| `--local`      | Run on the current machine instead of submitting to the cluster |
| `--jobs / -j`  | Number of jobs to run/submit in parallel (default 1)          |
| `-p`           | SLURM partition (cluster mode only)                           |

### 4. Outputs

Everything is written under `outDir`:

| Path                                                   | What it is                          |
|--------------------------------------------------------|-------------------------------------|
| `clean_reads/{sample}/{sample}.1.fq.gz`, `.2.fq.gz`    | Adapter/quality-trimmed reads       |
| `bowtie/{sample}/{sample}.bam` (+ `.bai`)              | Sorted, indexed alignments          |
| `bowtie_featurecounts/{sample}/{sample}.count.txt`     | Per-sample featureCounts table      |
| `logs/`                                                | Per-step logs                       |

### Things to know about this workflow

- **QC (bbduk):** With `qc: true`, reads go through adapter trimming and phiX
  removal, plus quality trimming controlled by `trimq`/`minlen`/`mink`/`mapq`.
  Set `qc: false` to skip trimming and align the reads as-is.
- **Alignment filtering:** Bowtie2 output is filtered with `samtools view -q10 -f2`,
  i.e. **only properly-paired reads with mapping quality ≥ 10 are kept**.
  Multi-mapping and low-confidence reads are discarded before counting.
- **`feature_type`/`attribute` must match your annotation.** For many bacterial
  GFFs the right combination is `feature_type: CDS`, `attribute: locus_tag` (or
  `gene_id`). Check your annotation's attribute names if counts come back empty.

---

## Other workflows (under construction)

The CLI exposes several other workflows. They share the same pattern — one
subcommand, one config file (`nccrRna <command> -c <config> [--dry|--local]`) —
but they are **not yet documented or fully validated**, and most are wired to the
ETH/NCCR cluster paths. Treat them as experimental. Minimal, documented config
templates (like the one for `metat-bowtie`) still need to be written for these.

| Command        | Tool / purpose                                | Status                          |
|----------------|-----------------------------------------------|---------------------------------|
| `star`         | STAR + featureCounts (eukaryotic)             | Plans (dry-runs) cleanly        |
| `salmon`       | Salmon transcript quantification (eukaryotic) | Plans cleanly                   |
| `prok`         | BWA + HTSeq-count (prokaryotic)               | Plans cleanly                   |
| `sushi`        | Sushicounter (metatranscriptomics)            | Plans cleanly                   |
| `metasalmon`   | Salmon in metagenomic mode                    | Plans cleanly                   |
| `rnafilter`    | rRNA removal with SortMeRNA                    | Plans cleanly; needs `rnadb`    |
| `mvirs`        | mVIRs prophage detection                      | Plans cleanly                   |
| `coptr`        | CoPTR replication-rate estimation             | Plans cleanly                   |
| `count-te`     | TEtranscripts                                 | Requires `teAnn` in the config  |
| `annotate`     | eggNOG-mapper + dbCAN + cayman                | Requires `.faa`/`.gff3` inputs and a `funanot` env |
| `motus`        | mOTUs profiling                               | **Broken** — target rule disabled |
| `count-reps`   | k-seek satellite-repeat counting              | **Broken** — target rule disabled |

Each workflow only requires the config keys relevant to it (the Snakefile
includes only that workflow's rules), so configs can stay small and
workflow-specific.

---

## Repository layout

```
nccr_rnaseq/
├── main.py                # the nccrRna CLI (one subcommand per workflow)
├── Snakefile              # workflow gating: includes only the selected workflow's rules
├── rules/                 # one .smk file per analysis step
├── envs/                  # conda environment recipes referenced by the rules
└── scripts/               # helper scripts (sample-sheet generation, etc.)
configs/                   # example config files
test_data/                 # small example inputs
```
