import argparse
from pybiomart import Dataset
from mylog import get_logger

logger = get_logger(__file__, __name__)

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
