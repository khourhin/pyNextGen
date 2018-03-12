import argparse
from collections import Counter
from mylog import get_logger
import pysam
#import pysamstats
from utils import apply_threads, simplify_path
import os

logger = get_logger(__file__, __name__)

# IN DVPT

# get_read_by_id is HIGHLY SUBOPTIMAL (found some index methods in
# pysam but did not manage to make it work to accelerate the acces to
# read by names)


def open_bam(bam):
    with pysam.AlignmentFile(bam, 'rb') as f:
        for i, read in enumerate(f):
            if i % 1000000 == 0:
                logger.info("Reached {} reads.".format(i))
            yield read

            
def get_polyNs_read(bam, poly_len, poly_char):
    """ Get reads with a polyN tail
    """
    
    for read in open_bam(bam):
        if read.seq.endswith(poly_char * poly_len):
            yield read

            
def get_read_by_id(id_list, bam):
    """
    Get the alignement from a bam file given an read id_list
    """

    for read in open_bam(bam):
        if read.query_name in id_list:
            yield read


def count_mapped_bases(bam):
    """
    Describe the proportion of each bases in a single read in primary
    alignments (we are trying here to segregate between polyAs and
    genomes As so need to remove the clipped part of the alignments
    before computing nucleotide proportions in reads)
    """

    for read in open_bam(bam):
        if not read.is_secondary:
            count = Counter(read.query_alignment_sequence)
            yield(count)


def overall_mapped_bases_composition(bam):

    all_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}

    for count in count_mapped_bases(bam):
        for k in all_counts:
            all_counts[k] += count[k]

    all_counts = {k: all_counts[k] /
                  sum(all_counts.values()) for k in all_counts}

    for k, v in all_counts.items():
        logger.info('{0}: {1}'.format(k, v))


### BROCKEN FOR NOW:
# Last checked, pysamstats was giving an error at import
        
def write_region_coverage(interval, bam, flank):
    """
    From a bam file, get the coverage for each bases in an interval
    (as in pyBedtools interval ie with interval.chrom, interval.start,
    interval.end)
    """

    sam = pysam.AlignmentFile(bam, 'rb')
    logger.debug("Interval: {}:{}-{}".format(interval.chrom,
                                             interval.start, interval.end))

    with open('coverage_' + simplify_path(bam) + "_" + interval.name, 'w') as f:
        for base_cov in pysamstats.stat_coverage(sam,
                                                 chrom=interval.chrom,
                                                 start=interval.start - flank,
                                                 end=interval.end + flank,
                                                 truncate=True, pad=True):
            f.write('{chrom}\t{pos}\t{reads_all}\t{reads_pp}\n'.format(**base_cov))

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='The path to a bam file')
    parser.add_argument('-b', '--base_comp',
                        help='Produce composition of the mapped bases',
                        action='store_true')
    
    args = parser.parse_args()

    if args.base_comp:
        overall_mapped_bases_composition(args.bam)

