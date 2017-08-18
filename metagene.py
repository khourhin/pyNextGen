import click
from basics_bam import write_region_coverage
from mylog import get_logger
import glob
import os
from pybedtools import BedTool
from multiprocessing import Pool


logger = get_logger(__file__, __name__)


def get_all_coverages(bams, bed, threads):
    """
    Get coverages for all regions in the 'bed' file for each bam
    present in 'bam_dir'
    """

    p = Pool(threads)
    
    logger.debug('Checking bams: {}'.format(bams))

    for bam in bams:
    
        p.starmap(write_region_coverage, ((interval, bam) for interval in bed))
            
    
@click.command()
@click.argument('bam_dir', type=click.Path(exists=True))
@click.argument('bed', type=click.Path(exists=True))
@click.option('--threads', type=int, default=1)
def main(**kwargs):
    
    logger.debug(kwargs)
    bed = BedTool(kwargs['bed'])
    bams = glob.glob(os.path.join(kwargs['bam_dir'], '*.bam'))
    
    get_all_coverages(bams, bed, kwargs['threads'])

    
if __name__ == '__main__':
    main()
