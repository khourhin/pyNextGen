#! /usr/bin/env python3

import mygene
from mylog import get_logger
import click
import sys

logger = get_logger(__file__, __name__)

# USAGE EXAMPLE (from stdin)
# cut -f2 gene_list_common_aslt_event_EB.tab | python ~/Programming/pyNextGen/identify_my_gene.py -

# WARNING SPECIFICI Working for hg19 now !!!

# For now working only with entrez or ensembl IDs, but can be much
# more with mg.querymany. But easily modified.


def post_output(res, kwargs):
    """Return the output to file or to stdout depending on the CLI options
    """
    if kwargs['outfile']:
        res.to_csv(kwargs['outfile'])
    else:
        res.to_csv(sys.stdout)

        
def annotate_genes(genes, kwargs=None):
    """Use entrez or ensembl ids list and print annotations for each ids
    """

    mg = mygene.MyGeneInfo()
    res = mg.getgenes(
        genes, fields='symbol,name,taxid,genomic_pos_hg19', as_dataframe=True)

    return res
    post_output(res, kwargs)


def search_genes(genes, kwargs):
    """Query genes based on common gene name (eg 'CD44')"""
    
    mg = mygene.MyGeneInfo()
    res = mg.querymany(
        genes, scopes='symbol', fields='ensembl.gene,symbol,name,taxid,genomic_pos_hg19', as_dataframe=True, species=kwargs['species'])
    post_output(res, kwargs)


@click.command()
@click.argument('genes', type=click.File('r'))
@click.option('--query', '-q', is_flag=True, help='For querying by common genes names (eg CD44)')
@click.option('--species', '-s', help='Specify the species to use (eg "human", "10090")')
@click.option('--outfile', '-o', type=click.File('w'), help='Output file')
def main(genes, query, **kwargs):

    if query:
        logger.info("Querying common gene names in {} to mygene.info".format(genes.name))
        search_genes(genes, kwargs)
    else:
        logger.info("Annotating entrez or ensembl ids in {} with mygene.info".format(genes.name))
        annotate_genes(genes, kwargs)
        
if __name__ == '__main__':
    main()
