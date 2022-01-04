import sys
from pathlib import Path


def parse_gtf(gtf_file, tx2file, gene_name="gene_id", transcript_name='transcript_id'):
    tx2gene = {}
    with open(gtf_file, 'r') as fh:
        for line in fh.readlines():
            if line.startswith("#"):
                continue
            gene_id = None
            tx_id = None
            attribute = line.split("\t")[8]
            try:
                if gene_name in line:
                    gene_id = attribute.split(gene_name)[1].split(';')[0].strip()
                else:
                    gene_id = attribute.split("gene_id")[1].split(';')[0].strip()
                tx_id = attribute.split(transcript_name)[1].split(';')[0].strip()
            except IndexError:
                continue
            finally:
                if gene_id and tx_id:
                    tx2gene[tx_id] = gene_id

    with open(tx2file, 'w') as fo:
        fo.write('TXNAME,GENEID\n')
        for k, v in tx2gene.items():
            fo.write(f'{k},{v}\n')

if __name__ == "__main__":
    gtf_file = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        gene_name = sys.argv[2]
    else:
        gene_name = 'gene_id'
    if len(sys.argv) == 4:
        transcript_name = sys.argv[3]
    else:
        transcript_name = 'transcript_id'
    parse_gtf(gtf_file, gtf_file.with_suffix((".tx2gene.csv")), gene_name, transcript_name)
