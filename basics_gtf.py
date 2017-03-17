import argparse
import logging
import gffutils
import os

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# STUB

def create_db(gtf, dbfn):
    """
    From a 'gtf' file, create a 'dbfn' sqlite database.
    """

    log.info('Creating db')

    gffutils.create_db(gtf, dbfn=dbfn, force=True, # Delete db if already exist
                       merge_strategy='merge',
                       id_spec={'exon': 'exon_id',
                                'gene': 'gene_id',
                                'transcript': 'transcript_id'},
                       disable_infer_transcripts=True,
                       disable_infer_genes=True)
    log.info('db done')



def filter_by_attribute(attr_in, dbfn):
    db = gffutils.FeatureDB(dbfn)

    for feat in db.all_features():
        if feat.attributes['gene_biotype'] != [attr_in]:
            continue
        
        print(feat)
    

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='')
    args = parser.parse_args()

    # Create db if not in current folder
    dbfn = os.path.basename(args.gtf) + '.db'
    if not os.path.isfile(dbfn):
        create_db(args.gtf, dbfn)
    
    filter_by_attribute('protein_coding', dbfn)

if __name__ == '__main__':
    
    main()
    
