rule te_counts:
    input: bam = OUTDIR/'bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
    output: touch(OUTDIR/'tetranscripts/{sample}/{sample}.te.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.te.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.te.qout',
        annotation = config["refAnn"],
        out_dir = lambda wildcards: OUTDIR/f'tetranscripts/{wildcards.sample}',
        te_annotation = config["teAnn"],
        threads = 32,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'te'
    log: OUTDIR/'logs/{sample}.te.log'
    threads:
        32
    shell: 
        "TEcount  -b {input.bam} --GTF {params.annotation} --TE {params.te_annotation} "
        "--outdir {params.out_dir} --sortByPos "


rule kseek:
    input: fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
    output: OUTDIR/'kseek/{sample}.rep.total'
    params:
        kdir = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/soft/k-seek",
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kseek.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kseek.qout',
        prefix = lambda wildcards: OUTDIR/f'kseek/{wildcards.sample}',
        threads = 16,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'call_variants'
    log: OUTDIR/'logs/{sample}.kseek.log'
    threads:
        16
    shell: 
        "zcat {input} |perl {params.kdir}/k_seek.pl - {params.prefix}"