# This is to remember that pybedtools exists !!
from pybedtools import BedTool
import random
import sys


# Dirty function ! TO cleanup !
def get_random_intervals(in_bed, n, min_size, max_size):
    """
    Get a bed with n random intervals extracted from the intervals
    present in in_bed. Random intervals generated will have a size
    between min_size and max_size
    """

    # Randomly select lines from the in_bed
    num_lines = sum(1 for line in open(in_bed, 'r'))
    line_sel = random.sample(range(num_lines), n)

    with open(in_bed, 'r') as f:
        lcount = 0
        for line in f:
            if lcount not in line_sel:
                lcount += 1
                continue

            chrom, start, end = line.split()[:3]

            # Define the random interval
            r_size = random.randint(min_size, max_size)
            try:
                r_start = random.randint(int(start), int(end) - r_size)
                r_end = r_start + r_size

                # If base interval smaller than desired length, return complete interval
            except ValueError as e:
                r_start = start
                r_end = end
            
            print("{0}\t{1}\t{2}".format(chrom, r_start, r_end))
            lcount += 1

if __name__ == "__main__":

    in_bed = sys.argv[1]
    
    get_random_intervals(in_bed, 2000, 520, 4469)
