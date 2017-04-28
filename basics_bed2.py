# Now using pybedtools
from pybedtools import BedTool
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

bam1 = '/home/ekornobis/analysis/demo/c_elegans/star_out/bams/SRR019721.bam'
bam2 = '/home/ekornobis/analysis/demo/c_elegans/star_out/bams/SRR019722.bam'


a = BedTool(bam1).bam_to_bed(bed12=True)
a = a.merge(s=True)

b= BedTool(bam2).bam_to_bed()
b = b.merge(s=True)

log.info('Bams import DONE')

c = a.intersect(b, s=True)

log.info('Intersection DONE')

log.info('DONE')
