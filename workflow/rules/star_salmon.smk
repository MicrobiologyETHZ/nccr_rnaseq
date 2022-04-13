from pathlib import Path

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])


rule STAR_index:
    input: config['refGenome']
    output: touch(Path(config['refGenome']).parent/f'{config["projectName"]}.star.index.done')
    params:
        qerrfile = OUTDIR/'logs/STAR.index.qerr',
        qoutfile = OUTDIR/'logs/STAR.index.qout',
        genomeDir = config["genomeDir"],
        genome_no_suffix = Path(config['refGenome']).parent/Path(config['refGenome']).stem,
        annotation = config["refAnn"],
        overhang = config["overhang"],
        genomeSAindexNbases = config["genomeSAindexNbases"],
        threads = 32,
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/'logs/STAR.index.log'
    threads:
        32
    shell:
        '''
        #!/bin/bash
        command="
        if [[ {input} = *.gz ]]; then
          gunzip {input}
          STAR --runThreadN {params.threads} --runMode genomeGenerate \
        --genomeFastaFiles {params.genome_no_suffix} --genomeDir {params.genomeDir} \
        --genomeSAindexNbases {params.genomeSAindexNbases} \
        --sjdbGTFfile {params.annotation} --sjdbOverhang {params.overhang}
          gzip {params.genome_no_suffix};
        else
          STAR --runThreadN {params.threads} --runMode genomeGenerate \
        --genomeFastaFiles {input} --genomeDir {params.genomeDir} \
        --genomeSAindexNbases {params.genomeSAindexNbases} \
        --sjdbGTFfile {params.annotation} --sjdbOverhang {params.overhang}
        fi
        ";
        eval "$command"
        '''


           #"--sjdbGTFtagExonParentTranscript Parent " # For GFF3 annotations

def star_input(wildcards):
    if not config['se']:
        return [f'{OUTDIR}/clean_reads/{wildcards.sample}/{wildcards.sample}.1.fq.gz',
                f'{OUTDIR}/clean_reads/{wildcards.sample}/{wildcards.sample}.2.fq.gz']
    else:
        return f'{OUTDIR}/clean_reads/{wildcards.sample}/{wildcards.sample}.1.fq.gz'

rule STAR_align:
    input: index_done = Path(config['refGenome']).parent/f'{config["projectName"]}.star.index.done',
        fq = star_input
    output: marker = touch(OUTDIR/'bam/{sample}/{sample}.done'),
         bam = OUTDIR/'bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.star_align.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.star_align.qout',
        files = lambda wildcards, input: " ".join(input[1]),
        genomeDir = config["genomeDir"],
        annotation = config["refAnn"],
        overhang = config["overhang"],
        maxIntron=config['maxIntron'],
        prefix = lambda wildcards: OUTDIR/f'bam/{wildcards.sample}/{wildcards.sample}_',
        threads = 8,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/'logs/{sample}.star_align.log'
    threads:
        32
    shell: "STAR --runThreadN {params.threads} "
           "--readFilesIn {input.fq} "
           "--readFilesCommand gunzip -c " #format for paired end
           "--genomeDir {params.genomeDir} "
           "--sjdbGTFfile {params.annotation} "
           "--sjdbOverhang {params.overhang} "
           "--outFileNamePrefix {params.prefix} "
           "--outSAMtype BAM SortedByCoordinate "
           "--outSAMunmapped Within "
           "--outSAMattributes Standard "
           "--alignIntronMax {params.maxIntron} "
           "--quantMode GeneCounts &> {log} "



rule featureCounts:
    input: bam = OUTDIR/'bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
    output: OUTDIR/'counts/{sample}/{sample}.count.txt'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qout',
        annotation = config["refAnn"],
        attribute = config["attribute"],
        strand = config["strand"],
        threads = 32,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/'logs/{sample}.featurecounts.log'
    threads:
        32
    shell:
        "featureCounts -p -T {params.threads} " 
        "-a {params.annotation} -o {output} " 
        "-g {params.attribute} {input.bam} -s {params.strand} &> {log}"



rule kallisto_index:
    input: tna = config['transcriptome']
    output: marker = touch(OUTDIR/f"{config['projectName']}.kallisto.index.done"),
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/kallisto_index.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/kallisto_index.qout',
        kallistoIdx = config["kallistoIdx"],
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        '../envs/kallisto.yaml'
    log: OUTDIR/'logs/kallisto_index.log'
    threads:
        32
    shell:
        "kallisto index -i {params.kallistoIdx} {input.tna} &> {log} "


rule kallisto_quant:
    input: index_done = OUTDIR/f"{config['projectName']}.kallisto.index.done",
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2= OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz'
    output: marker = touch(OUTDIR/'kallisto/{sample}.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kallisto_quant.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kallisto_quant.qout',
        kallistoIdx = config["kallistoIdx"],
        prefix = lambda wildcards: OUTDIR/f'kallisto/{wildcards.sample}',
        scratch = 6000,
        threads = 8,
        mem = 8000,
        time = 1400
    conda:
        '../envs/kallisto.yaml'
    log: OUTDIR/'logs/{sample}.kallisto_quant.log'
    threads:
        32
    shell:
         "kallisto quant -i {params.kallistoIdx} -o {params.prefix} "
         "-t {params.threads} -b 100 <(zcat {input.fq1}) <(zcat {input.fq2}) &> {log}"


rule salmon_index:
    input: tna = config['transcriptome'],
        genome = config['refGenome']
    output: marker = touch(Path(config['transcriptome']).parent/f"{config['projectName']}.salmon.index.done"),
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/salmon_index.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/salmon_index.qout',
        salmonIdx = f"{Path(config['transcriptome']).parent}/{Path(config['transcriptome']).stem}_salmon_index",
        gentrome = f"{Path(config['transcriptome']).parent}/{Path(config['transcriptome']).stem}.gentrome.gz",
        decoys = f"{Path(config['transcriptome']).parent}/decoys.txt",
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/'logs/salmon_index.log'
    threads:
        32
    shell:
        'grep "^>" <(gunzip -c {input.genome}) | cut -d " " -f 1 > {params.decoys};'
        'sed -i.bak -e \'s/>//g\' {params.decoys}; '
        'cat {input.tna} {input.genome} > {params.gentrome}; '
        'salmon index -t {params.gentrome} -d {params.decoys} -i {params.salmonIdx} -p 32 &> {log}'


if not config['se']:
    rule salmon_run:
        input: fq1 = OUTDIR/"clean_reads/{sample}/{sample}.1.fq.gz",
            fq2= OUTDIR / "clean_reads/{sample}/{sample}.2.fq.gz",
            indx_marker = Path(config['transcriptome']).parent/f"{config['projectName']}.salmon.index.done"
        output:  OUTDIR/"salmon/{sample}_quant/quant.sf"
        params:
            qerrfile = lambda wildcards: OUTDIR/f'logs/salmon_index.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/salmon_index.qout',
            salmonIdx = f"{Path(config['transcriptome']).parent}/{Path(config['transcriptome']).stem}_salmon_index",
            out_dir = lambda wildcards: OUTDIR/f'salmon/{wildcards.sample}_quant',
            scratch = 6000,
            mem = 8000,
            time = 1400
        conda:
            'star_salmon'
        log: OUTDIR/'logs/{sample}_quant.log'
        threads:
            32
        shell:
            "salmon quant -i {params.salmonIdx} -l A -1 {input.fq1} -2 {input.fq2} "
            "-p 8 --validateMappings --gcBias  -o {params.out_dir} &> {log}"

else:
    rule salmon_run:
        input: fq1=OUTDIR / "clean_reads/{sample}/{sample}.1.fq.gz",
            indx_marker=Path(config['transcriptome']).parent/f"{config['projectName']}.salmon.index.done"
        output: OUTDIR / "salmon/{sample}_quant/quant.sf"
        params:
            qerrfile=lambda wildcards: OUTDIR / f'logs/salmon_index.qerr',
            qoutfile=lambda wildcards: OUTDIR / f'logs/salmon_index.qout',
            salmonIdx=f"{Path(config['transcriptome']).parent}/{Path(config['transcriptome']).stem}_salmon_index",
            out_dir=lambda wildcards: OUTDIR / f'salmon/{wildcards.sample}_quant',
            scratch=6000,
            mem=8000,
            time=1400
        conda:
            'star_salmon'
        log: OUTDIR /'logs/{sample}_quant.log'
        threads:
            32
        shell:
            "salmon quant -i {params.salmonIdx} -l A -r  {input.fq1} "
            "-p 32 --validateMappings --gcBias  -o {params.out_dir} &> {log}"

