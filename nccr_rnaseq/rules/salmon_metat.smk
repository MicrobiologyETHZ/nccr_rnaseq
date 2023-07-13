rule salmon_index_metat:
    input:
        genome = config['refGenome']
    output: marker = touch(f'{config["refGenome"]}.salmon_index.done'),
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{Path(config["refGenome"]).stem}_salmon_index.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{Path(config["refGenome"]).stem}_salmon_index.qout',
        salmonIdx = f"{Path(config['refGenome']).parent}/{Path(config['refGenome']).stem}_salmon_index",
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'star_salmon'
    log: OUTDIR/f'logs/{Path(config["refGenome"]).stem}_salmon_index.log'
    threads:
        32
    shell:
        'salmon index -t {input.genome} -i {params.salmonIdx} -p 32 &> {log}'


rule salmon_run_metat:
        input:
            fwd = OUTDIR/'rnasorted/{sample}/{sample}.norna_fwd.fq.gz',
            rev = OUTDIR/'rnasorted/{sample}/{sample}.norna_rev.fq.gz',
            indx_marker = f'{config["refGenome"]}.salmon_index.done'
        output:  OUTDIR/"salmon_metat/{sample}_quant/quant.sf"
        params:
            qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}_salmon_metat.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}_salmon_metat.qout',
            salmonIdx = f"{Path(config['refGenome']).parent}/{Path(config['refGenome']).stem}_salmon_index",
            out_dir = lambda wildcards: OUTDIR/f'salmon_metat/{wildcards.sample}_quant',
            scratch = 6000,
            mem = 8000,
            time = 1400
        conda:
            'star_salmon'
        log: OUTDIR/'logs/{sample}_salmon_quant.log'
        threads:
            32
        shell:
            "salmon quant -i {params.salmonIdx} -l A -1 {input.fwd} -2 {input.rev} "
            "-p 8 --validateMappings --gcBias  -o {params.out_dir} &> {log}"