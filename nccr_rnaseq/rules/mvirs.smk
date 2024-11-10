sampleInfo = pd.read_csv(config['samples'])



SAMPLES = pd.read_csv(config['samples'])['sample'].unique()


def getFastq1(wildcards):
    return sampleInfo[sampleInfo['sample'] == wildcards.sample].fastq_1.values

def getFastq2(wildcards):
    return sampleInfo[sampleInfo['sample'] == wildcards.sample].fastq_2.values


print(sampleInfo[sampleInfo['sample']== 'STAU21-2_21MS06-EcoHS-Glu-7-0-ae-37-1_ISOG_subsample'].fastq_2.values)

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])

rule mvirs_index:
    input: config['refGenome']
    output: marker = touch(f'{config["refGenome"]}.mvirs.index.done')
    params:
        qerrfile = f'{config["refGenome"]}.mvirs.qerr',
        qoutfile = f'{config["refGenome"]}.mvirs.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'mvirs'
    threads:
        16
    log:
        OUTDIR/f'logs/{Path(config["refGenome"]).stem}.mvirs.log'
    shell: "mvirs index -f {input} &> {log}"


rule mvirs_opr:
    input:
        #fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        #fq2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
        fq1= getFastq1,
        fq2 = getFastq2,
        index_done = f'{config["refGenome"]}.mvirs.index.done',
    output:
        marker = touch(OUTDIR/'mvirs_out/{sample}.mvirs_oprs.done')
    params:
        refGenome = config['refGenome'],
        out_prefix = lambda wildcards: OUTDIR/f'mvirs_out/{wildcards.sample}',
        qerrfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.mvirs_oprs.qerr',
        qoutfile =  lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.mvirs_oprs.qout',
        db = config['refGenome'],
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.coptr_map.log'
    conda:
        'mvirs'
    threads:
        16
    shell:
        "mvirs oprs -f {input.fq1} -r {input.fq2} -db {params.refGenome} -t 16 -ml 300 -o {params.out_prefix}"
