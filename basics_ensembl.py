#! /usr/bin/env python3

import argparse
import re
from pyensembl import EnsemblRelease
from mylog import get_logger

logger = get_logger(__file__, __name__)

### THIS IS ONLY FOR HUMAN release 75 !!!
ensembl_release = EnsemblRelease(release=75, species="homo_sapiens")


def parse_id_name_list(mixed_list):
    """
    From a mixed list of ensembl IDs and associated gene names,
    retrieve information from the corresponding gene
    """
    logger.info('Total number of genes in the list: {}'.format(len(mixed_list)))

    bad_names_count = []
    multiple_ids_per_name = []
    a = re.compile('ENSG\d{11}')

    for gene_id in mixed_list:
        # the item is an ensembl gene id

        # Check if the name is an Ensembl gene ID (ENSG)
        if a.match(gene_id):
            g_entry = ensembl_release.gene_by_id(gene_id)

        # Else if the name is an associated gene name
        else:

            try:
                g_entry = ensembl_release.genes_by_name(gene_id)

                # Check if unique match associated gene name / ensembl ID
                if len(g_entry) == 1:
                    g_entry = g_entry[0]
                else:
                    multiple_ids_per_name.append(gene_id)
                    g_entry = None

            except ValueError:
                g_entry = None
                bad_names_count.append(gene_id)

        yield g_entry

    logger.warn('Number of genes failing entry lookup (probably alias gene problem): {0} {1} [DISCARDED]'.format(len(bad_names_count),
                                                                                                              bad_names_count))
    logger.warn('Number of gene names associated with more than 1 ensembl ID: {0} {1} [DISCARDED]'.format(len(multiple_ids_per_name),
                                                                                                       multiple_ids_per_name))


def print_gene_list(gene_entries):
    """
    Print the results obtained by parse_ID_name_list
    """

    for gene in gene_entries:
        if gene:
            print(gene.contig, gene.start, gene.end,
                  gene.strand, gene.gene_id, gene.name, sep='\t')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('gene_names',
                        help='A file with one string per line with either \
                        Ensembl IDs or associated gene names')
    args = parser.parse_args()

    logger.warning('This script is using h.sapiens release 75')
    
    gene_entries = parse_id_name_list(
        [x.strip() for x in open(args.gene_names)])
    print_gene_list(gene_entries)
    

if __name__ == '__main__':
    main()
