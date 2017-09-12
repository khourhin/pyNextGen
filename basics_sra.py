#! /usr/bin/env python3

import subprocess
from multiprocessing import Pool
from functools import partial

import click

import basics_ensembl as be
from mylog import get_logger


logger = get_logger(__file__, __name__)

# USAGE EX:
#for sra in $SRAS; do python ~/code/batsche/methyleu/get_gene_specific_sra.py $sra gene_list_file -f 10000; done


def get_bam_by_sra_and_gene(sra, gene, flank):
    """From a SRA accession and an couple gene_id/gene/entry from the
    basics_ensembl.parse_id_name_list function, fetch the
    corresponding bam interval.

    """
    logger.debug('Fetching SRA: {0}, for gene: {1} with flanking of {2}'.format(sra, gene, flank))
    
    # sratool-kit use ncbi "chr" notation for chromosome names
    # So add it if not present
    if not gene.contig.startswith('chr'):
        chro = 'chr' + gene.contig
        
    bam = "{sra}_{gene_name}.bam".format(sra=sra, gene_name=gene.name)
    cmd = 'sam-dump --aligned-region {chr}:{start}-{end} {sra} | samtools view -bh > {bam}'.format(chr=chro,
                                                                                                   start=gene.start - flank,
                                                                                                   end=gene.end + flank,
                                                                                                   bam=bam,
                                                                                                   sra=sra)
    logger.debug(cmd)
    subprocess.check_output(cmd, shell=True)


def get_all_bams(sras, gene_file, flank, threads):
    """
    """
    p = Pool(threads)
    
    genes = [g.strip() for g in open(gene_file, 'r')]

    gene_entries = be.parse_id_name_list(genes)

    for gene in gene_entries:

        p.map(partial(get_bam_by_sra_and_gene, gene=gene, flank=flank), sras)


@click.command()
@click.argument('sras', nargs=-1)
@click.argument('gene_file')
@click.option('--flank', type=click.INT, help='The number of bases to add on each side of the gene', default=0)
@click.option('--threads', type=click.INT, help='Number of threads', default=1)
def main(sras, gene_file, flank, threads):
    """Download from the sra database the data corresponding to the
    accession(s) SRA for the genes present in the GENE_LIST
    file. GENE_LIST file is a txt file with a list of Ensembl or Gene
    Symbol IDs (one per line).
    """
    logger.info("Fetching data for {0} SRAs librairies:{1}. Using {2} thread(s). Flanking region: {3}b.".format(len(sras), sras, threads, flank))
    get_all_bams(sras, gene_file, flank, threads)

if __name__ == "__main__":
    main()
