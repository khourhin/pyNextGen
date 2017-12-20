import logging
import os

from itertools import combinations
from multiprocessing.dummy import Pool
import subprocess
import pandas as pd

from functools import partial

from mylog import get_logger

logger = get_logger('', __name__, logging.DEBUG)


class Bed(object):
    """Bed object 
    """
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)

    def __len__(self):
        return sum(1 for _ in open(self.path))

    def __repr__(self):
        return '<Bed object: {}>'.format(self.name)

    def head(self, nlines=10):
        with open(self.path, 'r') as f:
            for i in range(nlines):
                print(f.readline().strip())

    def to_dataframe(self, index_col=0):
        return pd.DataFrame.from_csv(self.path, sep='\t', header=None, index_col=index_col)
    
    def total_bases(self):
        nbases = 0
        with open(self.path) as bed:
            for line in bed:
                line = line.split()
                nbases += int(line[2]) - int(line[1])

        return nbases

    def stats(self):
        """
            Produce all stats for a bed file
        """
    
        return {'Nbases': self.total_bases(),
                'Nintervals': len(self),
                'path': self.path,
                'name': self.name,
                'bedobj': self}


    def intersect(self, bed_obj, outfolder='bed_outfolder', supp_args=''):
        """ Make default bedtools intersect of two beds"""

        os.makedirs(outfolder, exist_ok=True)
        inter_path = os.path.join(outfolder,
                                  self.name + '-inter-' + bed_obj.name)
        cmd = 'bedtools intersect -a {0} -b {1} {2} | bedtools sort > {3}'.format(self.path,
                                                                  bed_obj.path,
                                                                  supp_args,
                                                                  inter_path)

        subprocess.call(cmd, shell=True)

        return Bed(inter_path)

    def subtract(self, bed_obj, outfolder='bed_outfolder', supp_args=''):
        """ Make default bedtools subtract of bed1 - bed2"""

        os.makedirs(outfolder, exist_ok=True)
        subtract_path = os.path.join(outfolder,
                                     self.name + '-minus-' + bed_obj.name)
        cmd = 'bedtools subtract -a {0} -b {1} {2} > {3}'.format(self.path,
                                                                 bed_obj.path,
                                                                 supp_args,
                                                                 subtract_path)

        subprocess.call(cmd, shell=True)

        return Bed(subtract_path)

    def fisher(self, bed_obj, genome_chro_size, supp_args=''):
        """ Make bedtools fisher test
        """

        cmd = 'bedtools fisher -a {0} -b {1} -g {2} {3}'.format(self.path,
                                                                bed_obj.path,
                                                                genome_chro_size,
                                                                supp_args)
        
        output = subprocess.check_output(cmd, shell=True).decode('utf-8')

        # Extract test values at the end of fisher output
        header = output.split('\n')[-3].split('\t')
        values = output.split('\n')[-2].split('\t')

        res = {k:float(v) for k,v in zip(header, values)}

        res['bed1'] = self.name
        res['bed2'] = bed_obj.name

        if len(res) != 6:
            raise IOError('Not correct fisher output for comparison:\n{}'
                          .format(cmd))
                
        return res

    def multi_shuffle(self, nshuffle, genome_chro_size, supp_args='', outfolder='bed_outfolder'):
        """ Bedtools shuffle multiple times
        """
        shuffle_beds = []
        
        for i in range(nshuffle):
            shuffle_path = os.path.join(outfolder, self.name + '-shuffle-' + str(i))
            cmd = 'bedtools shuffle -i {0} -g {1} {2} > {3}'.format(self.path,
                                                                    genome_chro_size,
                                                                    supp_args,
                                                                    shuffle_path)

            logger.debug('Running: {}'.format(cmd))
            subprocess.call(cmd, shell=True)
            shuffle_beds.append(Bed(shuffle_path))

        return shuffle_beds

    def coverage(self, bed_obj, outfolder='bed_outfolder', supp_args=''):
        """Bedtools coverage
        """
        coverage_path =  os.path.join(outfolder,
                                      self.name + '-coveredby-' + bed_obj.name)
        cmd = 'bedtools coverage -a {0} -b {1} > {2}'.format(self.path,
                                                             bed_obj.path,
                                                             coverage_path)

        subprocess.call(cmd, shell=True)

        return Bed(coverage_path)
        
    
def merge_bed_list(bed_tuple, outfolder='bed_outfolder'):
    """Concatenate, sort and merge a list of bed files, removing one bed
    at the time.  Name the resulting bed according to the only which
    was not in the concatenation.

    """

    del_bed_name = bed_tuple[0]
    bed_list = bed_tuple[1]
    
    merge_path = os.path.join(outfolder,
                              'all-merged-except-' + del_bed_name)
    
    cmd = 'cat {0} | bedtools sort | bedtools merge > {1}'\
          .format(' '.join([bed.path for bed in bed_list]), merge_path)
    
    subprocess.call(cmd, shell=True)
    
    return Bed(merge_path)


def get_all_merged_beds(bed_list, nthreads=1, outfolder='bed_outfolder'):
    """Threaded merge_bed_list
    """

    os.makedirs(outfolder, exist_ok=True)

    pool = Pool(nthreads)

    subsetted_bed_lists = {}
    
    for i, bed in enumerate(bed_list):

        # Remove one bed at the time and pass the list of beds to concatenate
        # bed_list[:] is to generate a copy to not pop out from actual bed_list
        l = bed_list[:]
        del_bed = l.pop(i)

        subsetted_bed_lists[del_bed.name] = l

    func = partial(merge_bed_list, outfolder=outfolder)
    
    res = pool.map(func, [(k, i) for k, i in subsetted_bed_lists.items()])

    return res


def jaccard_index(inter_stats_df, raw_beds_stats_df, count_column='Nbases'):
    """Function to get the jaccard index for each intersection, based on
    the values specified in column specified in count_column"""

    res = []

    for index, row in inter_stats_df.iterrows():
        bed1 = row.bed1
        bed2 = row.bed2
    
        bed1_count = raw_beds_stats_df.at[bed1, count_column]
        bed2_count = raw_beds_stats_df.at[bed2, count_column]
    
        jaccard_index = row[count_column] / (bed1_count + bed2_count - row[count_column])
    
        res.append({'jaccard': jaccard_index, 'bed1': bed1, 'bed2': bed2})

    return pd.DataFrame(res)


################################################################################
# FOR LEGACY, parallel execution is now done in the notebook
################################################################################

# def parallel_bedop(beds_list, method, nthreads=1):
#     """ Apply with threads the method for a list of Bed object
#     """
    
#     p = Pool(nthreads)
#     res = []

#     if method == 'intersect':

#         #This is a trick to be able to use the object method in threads
#         def threaded_intersect(bed1, bed2): return bed1.intersect(bed2)
#         res = p.starmap(threaded_intersect, combinations(beds_list, 2))

#     if method == 'subtract':
#         for i, bed in enumerate(beds_list):
#             # Copy the list
#             l = beds_list[:]

#             # Remove the bed considered to apply method against all other beds
#             l.pop(i)
#             res += p.map(getattr(bed, method), [vs_bed for vs_bed in l])

#     return res

################################################################################
