# Eukaryotic RNAseq Pipeline

### Installing the pipeline

1. Clone the git repo.
2. Create and activate the conda environment (highly recommend using `mamba`)

```shell
mamba env create -f rnaseq_environment.yaml
conda activate rnaseq
```
3. Install the pipeline

```shell
pip install -e .
```

4. Create the environments (where env_name is the same as the name of the yaml file)
Environment files are stored in `workflow/envs`

```shell
conda config --append envs_dirs <directory_to_store_the_envs>
mamba env create -f env_name.yaml \
--prefix <directory_to_store_theenvs>/env_name
conda config --set env_prompt '({name})'
```


6. Or for internal use:
```shell
conda config --append envs_dirs /nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/conda_envs
```

### Data Input

- Sequencing data.
- Output directory has to exist.
- Config file. Edit config to provide input paths for genomes of interest. Example configs can be found in `configs/` directory.
- Sample sheet: Sample names and paths to each of the files (forward and reverse) in `csv` format. 
Example:

| sample  | unit  | fastq_1                  | fastq_2                 |
|---------|-------|--------------------------|-------------------------|
| Sample1 | L0001 | ./Sample1_L0001_R1.fq.gz | ./Sample1_L001_R2.fq.gz |


- Can be generated using `nccrRna samples -c config.yaml` or

=======
```shell

nccrRna samples -i <FASTQ_DIR> -o <OUTPUT_DIR> -r1 <FWD_SUFFIX> -r2 <RVR_SUFFIX> -sn -sd _ -si 1 

```

Where:
`-r1` and `-r2` are forward and reverse extensions (e.g. _R1.fq.gz and _R2.fq.gz). Add `-sn` if you want to clean up the sample name, i.e. extract sample name from sequencing file name.
If `-sn` is used, need to specify `-sd` (delimiter to split the file name on), and `-si` (index of the last of the elements to be included in sample name)

Example: File name: Sample1_L0001_R1.fq.gz, to extract `Sample1` as sample name specify `-sn -sd _ -si 1`

This code was adopted from the one used by [nf-core rnaseq pipeline](https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py)


### Running the pipeline

0. Preprocessing. All raw sequencing data goes through a preprocessing step. See [here](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html) for more information.

1. First option is to use STAR aligner to align preprocessed reads and featureCounts from subread package to count the number of inserts for each feature (i.e. gene).
You will need the genome sequence in FASTA format, as well as annotation file in gtf format. For many model organisms these can be dowloaded [here]().

```shell
nccrRna star -c <config_file>
```

- By default will try to submit jobs to SGE queue system, add `--local` to run locally. 
- Add `--dry` to see the commands to be run, without running them
- The output will be 

STAR and featureCounts options can be specified in the config file and are discussed in more detail [here](https://methods-in-microbiomics.readthedocs.io/en/latest/transcriptomics/transcriptomics.html)

2. Second option is to use salmon to quantify isoform abundances. For this in addition to genome sequence, you will also need a transcriptome in FASTA format. 


```shell
rnapipe salmon -c <config_file>
```

- By default will try to submit jobs to SGE queue system, add `--local` to run locally. 
- Add `--dry` to see the commands to be run, without running them
- The output will be `salmon` directory, containing a subdirectory for each sample. The quantification file will be in `salmon/Sample1_quant/quant.sf`




