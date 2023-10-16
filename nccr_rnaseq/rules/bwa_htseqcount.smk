from pathlib import Path

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])


rule bwa_index:
    input: config['refGenome']
    output: f"{config['refGenome']}.bwt",
            marker = touch(f'{config["refGenome"]}.bwa_index.done')
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

rule bwa_index_transcriptome:
    input: config['transcriptome']
    output: f"{config['transcriptome']}.bwt",
            marker = touch(f'{config["transcriptome"]}.bwa_index.done')
    params:
        qerrfile = f'{config["transcriptome"]}.bwa.qerr',
        qoutfile = f'{config["transcriptome"]}.bwa.qout',
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
        index_done = f'{config["refGenome"]}.bwa_index.done',
    output:
        marker = touch(OUTDIR/'bwa/{sample}/{sample}.bwa.done'),
        bam = OUTDIR/'bwa/{sample}/{sample}.bam',
        bai = OUTDIR/'bwa/{sample}/{sample}.bam.bai'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.bwa.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.bwa.qout',
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



rule bwa_for_sushi_fwd:
    input:
        fq1 = OUTDIR/'rnasorted/{sample}/{sample}.norna_fwd.fq.gz',
        #fq2 = OUTDIR/'rnasorted/{sample}/{sample}.norna_rev.fq.gz',
        index_done = f'{config["transcriptome"]}.bwa_index.done',
    output:
        r1 = OUTDIR/'sushibam/{sample}/{sample}.r1.bam',
       # r2 = OUTDIR/'sushibam/{sample}/{sample}.r2.bam
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa_fwd.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa_fwd.qout',
        refGenome=config["transcriptome"],
        scratch=6000,
        mem=7700,
        time=1400
    log:
        log = OUTDIR / 'logs/{sample}.sushi_bwa_fwd.log'
    conda:
        'sushi'
    threads:
        16
    shell:
        "bwa mem -a -t {threads} {params.refGenome} {input.fq1} 2>> {log.log}| samtools view -F 4 -bh - > {output.r1} "
        #"bwa mem -a -t {threads} {params.refGenome} {input.fq2} 2>> {log.log}| samtools view -F 4 -bh - > {output.r2} "

rule bwa_for_sushi_rev:
    input:
        fq2 = OUTDIR/'rnasorted/{sample}/{sample}.norna_rev.fq.gz',
        index_done = f'{config["transcriptome"]}.bwa_index.done',
    output:
       r2 = OUTDIR/'sushibam/{sample}/{sample}.r2.bam'
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa_rev.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa_rev.qout',
        refGenome=config["transcriptome"],
        scratch=6000,
        mem=7700,
        time=1400
    log:
        log = OUTDIR / 'logs/{sample}.sushi_bwa_rev.log'
    conda:
        'sushi'
    threads:
        16
    shell:
        "bwa mem -a -t {threads} {params.refGenome} {input.fq2} 2>> {log.log}| samtools view -F 4 -bh - > {output.r2} "


rule sushi_merge:
    input:
        r1 = OUTDIR / 'sushibam/{sample}/{sample}.r1.bam',
        r2 = OUTDIR / 'sushibam/{sample}/{sample}.r2.bam'
    output:
        bam = OUTDIR / 'sushibam/{sample}/{sample}.bam'
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_bwa.qout',
        scratch=6000,
        mem=7700,
        time=1400
    conda:
        'sushi'
    threads:
        16
    log:
        log=OUTDIR / 'logs/{sample}.sushi_bwa.log'
    shell:
        "sushicounter  mergebam -i1 {input.r1} -i2 {input.r2}  -t {threads} -o {output.bam} -i 0.95 -a 45 -c 0.8 -u &> {log.log} "

rule sushi_count:
    input:
        bam = OUTDIR / 'sushibam/{sample}/{sample}.bam'
    output:
        countsfile = OUTDIR/'sushibam/{sample}/{sample}.shushicounts'
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_count.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}.sushi_count.qout',
        inserts = 1,
        bases = 1,
        scratch=6000,
        mem=7700,
        time=1400
    conda:
        'sushi'
    threads:
        16
    log:
        log=OUTDIR / 'logs/{sample}.sushi_count.log'
    shell:
        "sushicounter counter -pe {input.bam} -o {output.countsfile} " 
        "-i 0.95 -c 0.8 -a 45 -m mem  &> {log.log} "

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