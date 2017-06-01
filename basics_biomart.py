import argparse
import logging as log
from pybiomart import Dataset

log.basicConfig(filename='example.log', level=log.INFO)
log.getLogger().addHandler(log.StreamHandler())

# SO FAR WORKING ONLY FOR HUMANS 

def get_ensembl_table():

    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    
    table = dataset.query(attributes=['ensembl_gene_id',
                                      'external_gene_name',
                                      'unigene'])
    return table

def filter_on(table, column='Unigene ID'):
    pass
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('', help='')
    args = parser.parse_args()




