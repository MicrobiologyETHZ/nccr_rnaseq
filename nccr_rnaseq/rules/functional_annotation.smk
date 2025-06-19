rule eggnog_mapper:
    input:
        faa = "{sample}.faa"
    output:
        annotations = OUTDIR / "eggnog/{sample}/{sample}.emapper.annotations",
        seed_orthologs = OUTDIR / "eggnog/{sample}/{sample}.emapper.seed_orthologs",
        done = touch(OUTDIR / "eggnog/{sample}/{sample}.eggnog.done")
    params:
        qerrfile = OUTDIR / "logs/{sample}.eggnog.qerr",
        qoutfile = OUTDIR / "logs/{sample}.eggnog.qout",
        output_dir = OUTDIR / "eggnog/{sample}",
        temp_dir = OUTDIR / "eggnog/{sample}/temp",
        data_dir = config.get('eggnog_data_dir', '/science/ansintsova/eggnog-data/'),
        sample_name = "{sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.eggnog.log"
    threads: 16
    conda: funanot
    shell:
        """
        mkdir -p {params.output_dir} {params.temp_dir}
        
        emapper.py -o {params.sample_name} \
            -m diamond \
            --data_dir {params.data_dir} \
            -i {input.faa} \
            --output_dir {params.output_dir} \
            --temp_dir {params.temp_dir} \
            --cpu {threads} \
            &> {log}
        
        """


rule dbcan_annotation:
    input:
        faa = "{sample}.faa",
        gff = "{sample}.gff3"
    output:
        overview = OUTDIR / "dbcan/output_{sample}/overview.txt",
        done = touch(OUTDIR / "dbcan/output_{sample}/{sample}.dbcan.done")
    params:
        qerrfile = OUTDIR / "logs/{sample}.dbcan.qerr",
        qoutfile = OUTDIR / "logs/{sample}.dbcan.qout",
        output_dir = OUTDIR / "dbcan/output_{sample}",
        db_dir = config.get('dbcan_db_dir', '/science/ansintsova/dbCAN/db'),
        sample_name = "{sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.dbcan.log"
    threads: 16
    conda: funanot
    shell:
        """
        mkdir -p {params.output_dir}
        
        run_dbcan {input.faa} protein \
            -c {input.gff} \
            --dbcan_thread {threads} \
            --tf_cpu {threads} \
            --stp_cpu {threads} \
            --dia_cpu {threads} \
            --hmm_cpu {threads} \
            --cgc_substrate \
            --out_dir {params.output_dir} \
            --db_dir {params.db_dir} \
            &> {log}
    
        """


# This has not been tried or tested
rule cayman_annotation:
    input:
        faa = "{sample}.faa"
    output:
        results = OUTDIR / "cayman/{sample}/{sample}.cayman",
        done = touch(OUTDIR / "cayman/{sample}/{sample}.cayman.done")
    params:
        qerrfile = OUTDIR / "logs/{sample}.cayman.qerr",
        qoutfile = OUTDIR / "logs/{sample}.cayman.qout",
        output_dir = OUTDIR / "cayman/{sample}",
        db_dir = config.get('cayman_db_dir', '/nfs/cds-1/IMB-databases/CAYMAN/v3'),
        cutoffs = config.get('cayman_db_dir', '/nfs/cds-1/IMB-databases/CAYMAN/v3') + "/cutoffs.csv",
        sample_name = "{sample}",
        threads = 16,
        scratch = 8000,
        mem = 16000,
        time = 1400
    log:
        OUTDIR / "logs/{sample}.cayman.log"
    threads: 16
    conda: funanot
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        
        cayman annotate_proteome \
            --cutoffs {params.cutoffs} \
            {params.db_dir} \
            {input.faa} \
            --threads {threads} \
            -o {params.sample_name}.cayman \
            &> {log}
        
        """

