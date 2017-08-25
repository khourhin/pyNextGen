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


def annotate_genes(genes):
    """Use entrez or ensembl ids list and print annotations for each
    """

    mg = mygene.MyGeneInfo()
    res = mg.getgenes(genes, fields='symbol,name,taxid,genomic_pos_hg19', as_dataframe=True)

    return(res.to_csv(sys.stdout))

    
@click.command()
@click.argument('genes', type=click.File('r'))
def main(genes):
    
    annotate_genes(genes)

if __name__ == '__main__':
    main()
