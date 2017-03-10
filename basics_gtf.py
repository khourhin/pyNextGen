import argparse
import logging
import gffutils

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# STUB

def create_db(gtf, dbfn='test.db'):
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



def filter_by_attribute(attr_in, dbfn='test.db'):
    db = gffutils.FeatureDB(dbfn)

    for feat in db.all_features():
        if feat.attributes['gene_biotype'] != [attr_in]:
            continue
        
        print(feat)
    

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='')
    args = parser.parse_args()

    # create_db(args.gtf)
    filter_by_attribute('protein_coding')

if __name__ == '__main__':
    
    main()
    
