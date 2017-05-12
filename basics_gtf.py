import argparse
import logging
import gffutils
import os

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


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

def get_exons(dbfn, genes):
    """
    Print out all the exons for a list of genes
    """
    db = gffutils.FeatureDB(dbfn)
    
    # Header
    print("exon_name\tgene_name\tensID\texonID\tchrom\tstart\tstop\tlength\tstrand")
 
    for gene in genes:
        log.debug("Printing exon for gene: {}".format(gene))
        ex_start = None
        ex_stop = None
        count = 0
        
        for exon in db.children(gene, featuretype='exon', order_by='start'):
            # skip if exon got the same start and stop than the previous entry
            if exon.start == ex_start and exon.stop == ex_stop:
                continue

            count += 1
            ex_start = exon.start
            ex_stop = exon.stop
            
            print("{exon_name}\t{gene_name}\t{ensID}\t{exonID}\t{chrom}\t{start}\t{stop}\t{length}\t{strand}".format(
                exon_name=exon.attributes['gene_name'][0] + "_exon_" + str(count),
                gene_name=exon.attributes['gene_name'][0],
                ensID=exon.attributes['gene_id'][0],
                exonID=exon.attributes['exon_id'][0],
                chrom=exon.chrom,
                start=exon.start,
                stop=exon.stop,
                length=exon.stop - exon.start + 1,
                strand=exon.strand
            ))
        # Or complete report
        # print(dir(exon))

            
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='')
    group = parser.add_mutually_exclusive_group()

    group.add_argument('--filter', '-f',
                       help='Specify a gene_biotype to select. Ex: protein_coding.')

    group.add_argument('--introns', '-i', action='store_true',
                       help='Print out all introns corresponding to the exons in the gtf file.')

    group.add_argument('--exons', '-e',
                       help='Print out all the exons for the list of genes specified.', nargs='+')
    args = parser.parse_args()

    # Create db if not in current folder
    dbfn = args.gtf + '.db'
    if not os.path.isfile(dbfn):
        create_db(args.gtf, dbfn)
    else:
        log.info("Using already available database: {}".format(dbfn))
    
    if args.filter:
        filter_by_attribute(args.filter, dbfn)
    if args.introns:
        get_introns(dbfn, 'introns.gtf')
    if args.exons:
        get_exons(dbfn, args.exons)


if __name__ == '__main__':
    
    main()
    
