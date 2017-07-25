# Now using pybedtools
from pybedtools import BedTool
import argparse
import logging as log
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

log.basicConfig(filename='example.log', level=log.INFO)
log.getLogger().addHandler(log.StreamHandler())


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

    
    [log.info('{0}:{1}'.format(k, v)) for k, v in stats.items()]

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', help='A bed file')
    parser.add_argument('genome_size', type=int, help='The total number of bases of the considered genome')
    args = parser.parse_args()

    bedObj = BedTool(args.bed)
    bed_stats(bedObj, args.genome_size)
