import argparse
import logging as log
import pysam
from collections import Counter
import pandas as pd

log.basicConfig(filename='example.log', level=log.INFO)
log.getLogger().addHandler(log.StreamHandler())


# IN DVPT

# get_read_by_id is HIGHLY SUBOPTIMAL (found some index methods in
# pysam but did not manage to make it work to accelerate the acces to
# read by names)


def openBam(bam):
    with pysam.AlignmentFile(bam, 'rb') as f:
        for read in f:
            yield read

            
def get_read_by_id(id_list, bam):
    """
    Get the alignement from a bam file given an read id_list
    """

    for read in openBam(bam):
        if read.query_name in id_list:
            yield read

            
def count_bases(bam):
    """
    Describe the proportion of A in primary alignments (we are
    trying here to segregate between polyAs and genomes As so need to
    remove the clipped part of the alignments before computing
    nucleotide proportions in reads)
    """
    for read in openBam(bam):
        if not read.is_secondary:
            count = Counter(read.query_alignment_sequence)
            yield(count)

            
def overall_bases_composition(bam):
    counts = pd.DataFrame(count_bases(bam))
    log.info('Total bases composition:\n{}'.format(counts.sum() / sum(counts.sum())))
    return counts

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='The path to a bam file')
    args = parser.parse_args()
    
    overall_bases_composition(args.bam)
