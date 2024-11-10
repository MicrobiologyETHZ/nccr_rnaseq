from pathlib import Path
import pandas as pd

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])
import sys


SAMPLES = pd.read_csv(config['samples'])['sample'].unique()

def get_fwd_files():
    fastq1 = [f"{OUTDIR}/'clean_reads/{sample}/{sample}.1.fq.gz" for sample in SAMPLES]
    fwd_libs = [f"--pe-1 {i} {fastq}" for i, fastq in enumerate(fastq1, 1)]
    return " ".join(fwd_libs)

def get_rvr_files():
    fastq2 = [f"{OUTDIR}/'clean_reads/{sample}/{sample}.2.fq.gz" for sample in SAMPLES]
    rev_libs = [f"--pe-2 {i} {fastq}" for i, fastq in enumerate(fastq2, 1)]
    return " ".join(rev_libs)


rule rnaspades:
    input:
        fwd_libs = get_fwd_files,
        rvr_libs = get_rvr_files,
    output: 
        marker = touch(OUTDIR/f'rnaspades/{config["projectName"]}.rnaspade.done')

    params:
        outdir = OUTDIR/'rnaspades',
        qoutfile = OUTDIR /f'logs/{config["projectName"]}.rnaspades.qout',
        qerrfile = OUTDIR /f'logs/{config["projectName"]}.rnaspades.qerr',
        scratch = 500,
        mem = 8000,
        time = 235
        conda:
            "assembly"
        log:
            log = OUTDIR /'logs/rnaspades.qc.log'
        threads:
            16
        shell:
           "rnaspades.py -t 16  "
            " {input.fwd_libs} {input.rvr_libs} "
            "-o {params.outdir} &> {log.log} "