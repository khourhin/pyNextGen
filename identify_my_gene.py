import mygene
import mylog
import click

logger = mylog.get_logger(__file__, __name__)

# WARNING SPECIFICI Working for hg19 now !!!
# But easily modified


def annotate_genes(genes):
    """Use entrez or ensembl ids list and print annotations for each
    """

    mg = mygene.MyGeneInfo()
    res = mg.getgenes(genes, fields='symbol,name,genomic_pos_hg19', as_dataframe=True)

    print(res)

    
@click.command()
@click.argument('genes', type=click.File('r'))
def main(genes):
    
    annotate_genes(genes)

if __name__ == '__main__':
    main()
