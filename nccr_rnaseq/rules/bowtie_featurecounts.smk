from pathlib import Path

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])

rule bowtie_index:
    input: config['refGenome']
    output: marker = touch(f'{config["refGenome"]}.bowtie.index.done')
    params:
        prefix = Path(config['refGenome']).stem,
        qerrfile = f'{config["refGenome"]}.bowtie.qerr',
        qoutfile = f'{config["refGenome"]}.bowtie.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'star_salmon'
    threads:
        16
    log:
        OUTDIR/f'logs/{Path(config["refGenome"]).stem}.bowtie.log'
    shell: "bowtie2-build {input} {input} &> {log}"


rule bowtie_align:
    input:
        fq1 = OUTDIR / 'rnasorted/{sample}/{sample}.norna_fwd.fq.gz',
        fq2 = OUTDIR / 'rnasorted/{sample}/{sample}.norna_rev.fq.gz',
        #fq1 = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
        #fq2 = OUTDIR /'clean_reads/{sample}/{sample}.2.fq.gz',
        index_done = f'{config["refGenome"]}.bowtie.index.done',
    output:
        marker = touch(OUTDIR/'bowtie/{sample}/{sample}.bowtie.done'),
        bam = OUTDIR/'bowtie/{sample}/{sample}.bam',
        bai = OUTDIR/'bowtie/{sample}/{sample}.bam.bai'
    params:
        qerrfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.bowtie.qerr',
        qoutfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.bowtie.qout',
        refGenome = config['refGenome'],
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.bowtie.log'
    conda:
        'star_salmon'
    threads:
        16
    shell:
        "bowtie2 -x  {params.refGenome} "
        "-1 {input.fq1} -2 {input.fq2}  | samtools view -q10 -f2 -bh - |samtools sort --reference {params.refGenome} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam};"



rule bowtie_featurecounts:
    input: bam = OUTDIR/'bowtie/{sample}/{sample}.bam',
    output: OUTDIR/'bowtie_featurecounts/{sample}/{sample}.count.txt'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qout',
        annotation = config["refAnn"],
        attribute = config["attribute"],
        strand = config["strand"],
        threads = 32,
        feature_type = config['feature_type'],
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/'logs/{sample}.bowtie_featurecounts.log'
    threads:
        32
    shell:
        "featureCounts -p -T {params.threads} "
        "-a {params.annotation} -o {output} -t {params.feature_type} "
        "-g {params.attribute} {input.bam} -s {params.strand} "
        " &> {log} "
