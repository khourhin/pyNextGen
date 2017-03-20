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


def get_introns(dbfn, out_gtf):
    
    db = gffutils.FeatureDB(dbfn)
    log.info('Calculating introns')
    introns = db.create_introns()

    # Write gff will not work here, so using a simple loop
    with open(out_gtf, 'w') as f:
        for i in introns:
            f.write(str(i) + '\n')


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--filter', '-f', action='store_true', help='')
    group.add_argument('--introns', '-i', action='store_true', help='')
    args = parser.parse_args()

    # Create db if not in current folder
    dbfn = os.path.basename(args.gtf) + '.db'
    if not os.path.isfile(dbfn):
        create_db(args.gtf, dbfn)

    if args.filter:
        filter_by_attribute('protein_coding', dbfn)
    if args.introns:
        get_introns(dbfn, 'introns.gtf')


if __name__ == '__main__':
    
    main()
    
