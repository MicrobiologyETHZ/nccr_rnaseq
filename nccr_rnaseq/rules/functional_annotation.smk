rule eggnog_mapper:
    input:
        faa = OUTDIR/"{sample}.faa"
    output:
        annotations = OUTDIR / "eggnog/{sample}/{sample}.emapper.annotations",
        seed_orthologs = OUTDIR / "eggnog/{sample}/{sample}.emapper.seed_orthologs",
        done = touch(OUTDIR / "eggnog/{sample}/{sample}.eggnog.done")
    params:
        qerrfile = lambda wildcards: OUTDIR / f"logs/{wildcards.sample}.eggnog.qerr",
        qoutfile = lambda wildcards: OUTDIR / f"logs/{wildcards.sample}.eggnog.qout",
        output_dir = lambda wildcards: OUTDIR / f"eggnog/{wildcards.sample}",
        temp_dir = lambda wildcards: OUTDIR / f"eggnog/{wildcards.sample}/temp",
        data_dir = config.get('eggnog_data_dir', '/science/ansintsova/eggnog-data/'),
        sample_name = lambda wildcards: f"{wildcards.sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.eggnog.log"
    threads: 16
    conda: "funanot"
    shell:
        """
        mkdir -p {params.output_dir} {params.temp_dir}
        
        emapper.py -o {params.sample_name} \
            -m diamond \
            --data_dir {params.data_dir} \
            -i {input.faa} \
            --output_dir {params.output_dir} \
            --temp_dir {params.temp_dir} \
            --cpu 16 \
            &> {log}
        
        """


rule dbcan_annotation:
    input:
        faa = OUTDIR/"{sample}.faa",
        gff = OUTDIR/"{sample}.gff3"
    output:
        overview = OUTDIR / "dbcan/output_{sample}/overview.txt",
        done = touch(OUTDIR / "dbcan/output_{sample}/{sample}.dbcan.done")
    params:
        qerrfile = lambda wildcards: OUTDIR / f"logs/{wildcards.sample}.dbcan.qerr",
        qoutfile = lambda wildcards: OUTDIR / "logs/{wildcards.sample}.dbcan.qout",
        output_dir = lambda wildcards: OUTDIR / f"dbcan/output_{wildcards.sample}",
        db_dir = config.get('dbcan_db_dir', '/science/ansintsova/dbCAN/db'),
        sample_name = lambda wildcards: f"{wildcards.sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.dbcan.log"
    threads: 16
    conda: "funanot"
    shell:
        """
        mkdir -p {params.output_dir}
        
        run_dbcan {input.faa} protein \
            -c {input.gff} \
            --dbcan_thread 16 \
            --tf_cpu 16 \
            --stp_cpu 16 \
            --dia_cpu 16 \
            --hmm_cpu 16 \
            --cgc_substrate \
            --out_dir {params.output_dir} \
            --db_dir {params.db_dir} \
            &> {log}
    
        """


# This has not been tried or tested
rule cayman_annotation:
    input:
        faa = OUTDIR/"{sample}.faa"
    output:
        results = OUTDIR / "cayman/{sample}/{sample}.cayman",
        done = touch(OUTDIR / "cayman/{sample}/{sample}.cayman.done")
    params:
        qerrfile = lambda wildcards: OUTDIR / f"logs/{wildcards.sample}.cayman.qerr",
        qoutfile = lambda wildcards: OUTDIR / f"logs/{wildcards.sample}.cayman.qout",
        output_dir = lambda wildcards:OUTDIR / f"cayman/{wildcards.sample}",
        db_dir = config.get('cayman_db_dir', '/nfs/cds-1/IMB-databases/CAYMAN/v3'),
        cutoffs = config.get('cayman_db_dir', '/nfs/cds-1/IMB-databases/CAYMAN/v3') + "/cutoffs.csv",
        sample_name = lambda wildcards: f"{wildcards.sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.cayman.log"
    threads: 16
    conda: "funanot"
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        
        cayman annotate_proteome \
            --cutoffs {params.cutoffs} \
            {params.db_dir} \
            {input.faa} \
            --threads 16 \
            -o {params.sample_name}.cayman \
            &> {log}
        
        """

