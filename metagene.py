import click
from basics_bam import write_region_coverage
from mylog import get_logger
import glob
import os
from pybedtools import BedTool
from multiprocessing import Pool
import pandas as pd
from utils import simplify_path
import numpy as np

logger = get_logger(__file__, __name__)


def get_all_coverages(bams, bed, flank, threads):
    """
    Get coverages for all regions in the 'bed' file for each bam
    present in 'bam_dir'
    """

    p = Pool(threads)
    
    logger.debug('Checking bams: {}'.format(bams))

    for bam in bams:
    
        p.starmap(write_region_coverage, ((interval, bam, flank) for interval in bed))
        logger.info("DONE Coverage for {}".format(bam))


def coverage_plot(cov_file, flank):

    df = pd.read_csv(cov_file, sep='\t', usecols=['reads_all'])
    df['bins'] = [x * 2 for x in df.reads_all if x.index]
    print(df)
    
    # bins = np.linspace(df.reads_all,)
    # pos = np.digitize(df.reads_all, bins)
    
    # df['pos'] = pos
    # df.columns = [simplify_path(cov_file), 'pos']
        
    
@click.command()
@click.argument('bam_dir', type=click.Path(exists=True))
@click.argument('bed', type=click.Path(exists=True))
@click.option('--threads', type=int, default=1, help='Number of threads to use')
@click.option('--flank', type=int, default=0, help='Number of bases to add for the flanking area')
def main(**kwargs):
    
    logger.debug(kwargs)
    bed = BedTool(kwargs['bed'])
    bams = glob.glob(os.path.join(kwargs['bam_dir'], '*.bam'))
    
#    get_all_coverages(bams, bed, kwargs['flank'], kwargs['threads'])

    cov_file = '/home/ekornobis/Programming/pyNextGen/testing/j1'
    coverage_plot(cov_file, kwargs['flank'])

    
if __name__ == '__main__':
    main()
