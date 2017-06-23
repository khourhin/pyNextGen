from pyensembl import EnsemblRelease
import re
import argparse
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


data = EnsemblRelease(release=75, species="homo_sapiens")


def parse_ID_name_list(mixed_list):
    """
    From a mixed list of ensembl IDs and associated gene names,
    retrieve information from the corresponding gene
    """
    log.info('Total number of genes in the list: {}'.format(len(mixed_list)))

    bad_names_count = []
    multiple_IDs_per_name = []
    a = re.compile('ENSG\d{11}')
    
    for gene_id in mixed_list:
        # the item is an ensembl gene id

        # Check if the name is an Ensembl gene ID (ENSG)
        if a.match(gene_id):
            g_entry = data.gene_by_id(gene_id)

        # Else if the name is an associated gene name
        else:
            
            try:
                g_entry = data.genes_by_name(gene_id)

                # Check if unique match associated gene name / ensembl ID
                if len(g_entry) == 1:
                    g_entry = g_entry[0]
                else:
                    multiple_IDs_per_name.append(gene_id)
                    g_entry = None
                
            except ValueError:
                g_entry = None
                bad_names_count.append(gene_id)

        yield g_entry
                
    log.warn('Number of genes failing entry lookup (probably alias gene problem): {0} {1} [DISCARDED]'.format(len(bad_names_count),
                                                                                                              bad_names_count))
    log.warn('Number of gene names associated with more than 1 ensembl ID: {0} {1} [DISCARDED]'.format(len(multiple_IDs_per_name),
                                                                                                       multiple_IDs_per_name))

    
def print_gene_list(ginfos):
    """
    Print the results obtained by parse_ID_name_list
    """
    
    for i in ginfos:
        if i:
            print(i.contig, i.start, i.end, i.strand, i.gene_id, i.name, sep='\t')

            
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gene_names', help='A file with one string per line with either Ensembl IDs or associated gene names')
    args = parser.parse_args()

    ginfos = parse_ID_name_list([x.strip() for x in open(args.gene_names)])
    print_gene_list(ginfos)

    
if __name__ == '__main__':
    main()
