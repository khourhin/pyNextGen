#! /usr/bin/env python3

# Now using pybedtools
from pybedtools import BedTool
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import click
import logging
import os
from itertools import combinations
from multiprocessing import Pool
import pandas as pd
from mylog import get_logger

logger = get_logger(__file__, __name__)


def get_genome_coverage(bedObj, genome_size):
    """
    Get an estimation of the proportion of the genome covered by all
    the intervals present in the bed object.
    """

    bedObj = bedObj.sort().merge()
    cov = sum([len(i) for i in bedObj]) / genome_size

    return cov


def bed_stats(bedObj, genome_size):
    """
    Produce various stats on a bed file
    """

    lengths = [len(i) for i in bedObj]
    pd.DataFrame({"lengths":lengths}).plot()
    plt.show()
     
    stats = {
        'length': len(bedObj),
        'smallest': min(lengths),
        'largest': max(lengths),
        'mean': np.mean(lengths),
        'genome_coverage': get_genome_coverage(bedObj, genome_size)
    }

    
    [logger.info('{0}:{1}'.format(k, v)) for k, v in stats.items()]

    
def count_intersections(bed_tuple):
    """Return a normalized (by the total number of feature of both beds)
    intersetion counts
    """

    bed1 = bed_tuple[0]
    bed2 = bed_tuple[1]
    
    logger.info('Comparing: {0} VS {1}'.format(bed1.fn, bed2.fn))

    bed_inter = bed1.intersect(bed2)

    # Normnalized number of intersections
    # (divided by the total number of features in both files)
    
    return(bed1.fn,
           bed2.fn,
           len(bed_inter) / (len(bed1) + len(bed2)))


def count_all_intersections(beds):

    p = Pool(20)
    inter_df = p.map(count_intersections, combinations(beds, 2))
    inter_df = pd.DataFrame(inter_df, columns=['bed1', 'bed2', 'normCounts'])

    return(inter_df)

    
@click.command()
@click.argument('beds', required=True, nargs=-1)
@click.option('--genome_size', type=int)
def main(beds, genome_size):

    beds = (BedTool(bed) for bed in beds)
    inter_df = count_all_intersections(beds)
    inter_df.to_csv('out2.csv')

    
if __name__ == "__main__":
    main()

    
    # FOR LEGACY, TO RECONNECT PAST FUNCTION (NOT WORKING NOW)
    
    # parser = argparse.ArgumentParser()
    # parser.add_argument('bed', help='A bed file')
    # parser.add_argument('genome_size', type=int, help='The total number of bases of the considered genome')
    # args = parser.parse_args()

    # bedObj = BedTool(args.bed)
    # bed_stats(bedObj, args.genome_size)
