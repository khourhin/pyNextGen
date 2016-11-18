import logging
import gffutils

# NOT FULLY WORKING (CANNOT GET THE GTF FILE WITH INTRONS)

# ASSUMPTIONS:
# This is set for ENSEMBL gtf files

# Source: http://pythonhosted.org/gffutils/

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def create_db(gtf, dbfn):
    """
    From a 'gtf' file, create a 'dbfn' sqlite database.
    """

    log.info('Creating db')

    gffutils.create_db(gtf, dbfn=dbfn, force=True,
                       merge_strategy='merge',
                       id_spec={'exon': 'exon_id',
                                'gene': 'gene_id',
                                'transcript': 'transcript_id'},
                       disable_infer_transcripts=True,
                       disable_infer_genes=True)
    log.info('db done')

    # The parameters are coming from:
    # http://pythonhosted.org/gffutils/ AND
    # https://www.biostars.org/p/152517/


def gffutils_demo(gene_id, db):
    """
    A set of commands to illustrate gffutils functionalities
    """
    # Extract a gene by "gene_id"
    gene = db[gene_id]

    # Show its attributes
    gene.attributes.items()
    
    # Get one
    gene.attributes['gene_name']

    # Get all exons for a gene:
    for i in db.children(gene, order_by='start'):
        print(i)


def write_gff(db, fn):
    """
    Write the gff database 'db' to a gff file 'fn'
    The gff files will be ordered by start coordinates.

    TO CHECK: gff or gtf produced ?
    """

    log.info('writing db to file {}'.format(fn))
    
    with open(fn, 'w') as f:
        for i in db.all_features(order_by="start"):
            f.write(str(i) + '\n')


if __name__ == '__main__':

    log.warning('This module is designed with ENSEMBL GTFs in mind. Do not know how this would work with UCSC ones.')

    # From a gtf file from ensembl create another gtf with introns
    # added

    gtf = 'demo_data/C_elegans_small.gtf'
    dbfn = 'c_elegans.db'
    out_gtf = 'c_elegans_introns.gtf'
    
    create_db(gtf, dbfn)
    db = gffutils.FeatureDB(dbfn)
    
    # Generate introns from the db
    log.info('Calculating introns')
    introns = db.create_introns()

    # Write gff will not work here, so using a simple loop
    with open(out_gtf, 'w') as f:
        for i in introns:
            f.write(str(i) + '\n')
    
