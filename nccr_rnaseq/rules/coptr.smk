

sampleInfo = pd.read_csv(config['samples'])

samples_to_merge = (sampleInfo.loc[sampleInfo.groupby('sample')
                    .unit.filter(lambda x: x.nunique() > 1).index]['sample']
                    .unique())

SAMPLES = pd.read_csv(config['samples'])['sample'].unique()


def getFastq1(wildcards):
    return sampleInfo_merged[sampleInfo_merged['sample'] == wildcards.sample].fastq_1


rule coptr_map:
    input:
        fq1 = getFastq1,
        index_done = f'{config["refGenome"]}.bowtie.index.done',
    output:
        #bam = OUTDIR/'coptr_bams/{sample}.bam',
        marker = touch(OUTDIR/'coptr_bams/{sample}.coptr_map.done')
    params:
        out_dir = OUTDIR/'coptr_bams',
        qerrfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.coptr_map.qerr',
        qoutfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.coptr_map.qout',
        refGenome = config['refGenome'],
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.coptr_map.log'
    conda:
        'coptr'
    threads:
        16
    shell:
        "coptr map {params.refGenome} {input.fq1} {params.out_dir} --threads 16"


# coptr 

rule coptr_extract_estimate:
    input: 
        markers = [OUTDIR/f'coptr_bams/{sample}.coptr_map.done' for sample in SAMPLES]

    output: marker = touch(OUTDIR/'coptr.done')
    params:
        in_dir = OUTDIR/'coptr_bams',
        out_dir = OUTDIR/'coptr_coverage_maps',
        qerrfile = OUTDIR/f'logs/coptr_extract_estimate.qerr',
        qoutfile = OUTDIR/f'logs/coptr_extract_estimate.qout',
        scratch = 6000,
        mem = 7700,
        time = 8400
    log:
        log1 = OUTDIR/'logs/coptr_extract.log',
        log2 = OUTDIR/'logs/coptr_estimate.log'
    conda:
        'coptr'
    threads:
        200
    shell:
        "coptr extract {params.in_dir} {params.out_dir} &> {log.log1}; "
        "coptr estimate {params.out_dir} coptr_estimates &> {log.log2} "


