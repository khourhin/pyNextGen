import click
from basics_bam import get_region_coverage
from mylog import get_logger
import glob
import os
from pybedtools import BedTool

logger = get_logger(__file__, __name__)


def get_regions_coverages(bams, bed):
    """
    Get coverages for all regions in the 'bed' file for each bam
    present in 'bam_dir'
    """
    logger.debug('Checking bams: {}'.format(bams))

    for bam in bams:
        for interval in bed:

            with open('test' + interval.name, 'w') as f:
            
                for i in get_region_coverage(interval, bam):
                    f.write(i)




    
@click.command()
@click.argument('bam_dir', type=click.Path(exists=True))
@click.argument('bed', type=click.Path(exists=True))
def main(**kwargs):
    
    logger.debug(kwargs)
    bed = BedTool(kwargs['bed'])
    bams = glob.glob(os.path.join(kwargs['bam_dir'], '*.bam'))

    get_region_coverage(bed[0], bams[0])
        
    
#    get_regions_coverages(kwargs['bam_dir'], bed)

    
if __name__ == '__main__':
    main()
