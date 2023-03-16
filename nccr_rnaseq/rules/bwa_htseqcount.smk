from pathlib import Path

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])

rule bwa_index:
    input: config['refGenome']
    output: f"{config['refGenome']}.bwt",
            marker = touch(f'{config["refGenome"]}.index.done')
    params:
        qerrfile = f'{config["refGenome"]}.bwa.qerr',
        qoutfile = f'{config["refGenome"]}.bwa.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'bwa_htseq'
    threads:
        8
    shell: "bwa index {input}"



rule bwa_align:
    input:
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
        index_done = f'{config["refGenome"]}.index.done',
    output:
        marker = touch(OUTDIR/'bwa/{sample}/{sample}.bwa.done'),
        bam = OUTDIR/'bwa/{sample}/{sample}.bam',
        bai = OUTDIR/'bwa/{sample}/{sample}.bam.bai'
    params:
        qerrfile = OUTDIR/'logs/{sample}.bwa.qerr',
        qoutfile = OUTDIR/'logs/{sample}.bwa.qout',
        refGenome = config["refGenome"],
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.bwa.log'
    conda:
        'bwa_htseq'
    threads:
        8
    shell:
        "bwa mem  -t 4 "
        "-M {params.refGenome} "
        "{input.fq1} {input.fq2} | samtools sort --reference {params.refGenome} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam};"


rule alignment_stats:
    input: OUTDIR/'bwa/{sample}/{sample}.bam'
    output: OUTDIR/'bwa/{sample}/{sample}.bam.stats'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.samtools.stats.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.samtools.stats.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.samtools.stats.log',
    conda:
        'bwa_htseq'
    threads:
        8
    shell:
         "samtools stats {input} > {output}"


rule htseq:
    input: OUTDIR/'bwa/{sample}/{sample}.bam'
    output: OUTDIR /'htseqcount/{sample}/{sample}.htseqcount.txt'
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.htseq.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.htseq.qout',
        scratch=6000,
        mem=7700,
        time=1400,
        gff=config["refGff"],
        strand=config["strand"],
        attribute=config["attribute"],
        feature=config['feature_type']
    conda: 'bwa_htseq'
    log:
        log = OUTDIR/'logs/{sample}.htseq-count.log'
    threads: 8
    shell:
        "htseq-count -f bam -r pos -s {params.strand} -a 10 -t {params.feature} "
        "-i {params.attribute} -m union {input} {params.gff} > {output}"