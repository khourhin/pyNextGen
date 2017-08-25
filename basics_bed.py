# remember that pybedtools exists !!
import random
import sys
import numpy as np
import matplotlib.pyplot as plt
from mylog import get_logger

logger = get_logger(__file__, __name__)


class BedInterval():
    def __init__(self, bedline_list):
        self.chrom = bedline_list[0]
        self.start = int(bedline_list[1])
        self.end = int(bedline_list[2])
        self.length = self.end - self.start

    def rand_cut(self, fragl):
        """
        Cut a random fragment fo length "fragl" from the interval
        """
        try:
            c_start = random.randint(self.start, self.end - fragl)
        except ValueError as e:
            raise(e)

        return BedInterval([self.chrom, c_start, c_start + fragl])

    def __repr__(self):
        return "{0}\t{1}\t{2}".format(self.chrom, self.start, self.end)


class BedFile():
    def __init__(self, bedfile):
        self.intervals = [BedInterval(line.split())
                          for line in open(bedfile, 'r')]
        self.name = bedfile

    def rand_interval(self, min_size=None):
        """
        Return an randomly chose interval from a bed list.  If min_size is
        specified, intervals will be sampled until they can form an
        interval of minimal size min_size
        """

        sel_int = random.choice(self.intervals)
        if min_size:
            while sel_int.length < min_size:
                sel_int = random.choice(self.intervals)

        return sel_int

    def get_stats(self):
        """
        Get stats on the bed file
        """
        inters_len = [inter.length for inter in self.intervals]
        stats = {
            'nInt': len(self.intervals),
            'minLen': min(inters_len),
            'maxLen': max(inters_len),
            'meanLen': np.mean(inters_len),
            'medianLen': np.median(inters_len),
        }
        return stats

    def __repr__(self):
        return "Bedfile Object for file {}".format(self.name)


if __name__ == "__main__":

    # From 2 beds:
    # Get the lengths of the intervals in ref_bed
    # Produce random intervals from sampled_bed with
    # same lengths as ref_bed intervals

    ref_bed = sys.argv[1]
    sampled_bed = sys.argv[2]

    ref_bed = BedFile(ref_bed)
    sampled_bed = BedFile(sampled_bed)

    for inter in ref_bed.intervals:
        sample = sampled_bed.rand_interval(inter.length)
        print(sample.rand_cut(inter.length))
