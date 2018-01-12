# A library to detect divergent transcripttion

import logging
from mylog import get_logger
import inspect
import os
import subprocess

logger = get_logger('', __name__, logging.DEBUG)


def bam_filtering(bam):
    """ Filters bam file based on mapping, mapping of mate, splicing
    """

    filt_bam = os.path.splitext(bam)[0] + "_filt.bam"
    
    # Filtering bam, removing unmapped, mate unmapped, spliced, secondary mapping
    cmd = "samtools view -h -F 0x4 -F 0x8 -F 0x100 {0} | awk '{{if ($1 ~ /^@/) {{print}} else if ($6 !~ /N/) {{print}}}}' | samtools view -bh > {1}".format(bam, filt_bam)

    subprocess.check_output(cmd, shell=True)
    logger.debug('DONE: {}'.format(cmd))

    return filt_bam


def bam_to_fragments(bam):
    """From mapped reads in a bam, get bed intervals of the pairs joined
    into fragments
    """
    bedpe = os.path.splitext(bam)[0] + "_frag.bed"
    # Converting to fragments and bed format
    cmd = "bedtools bamtobed -bedpe -mate1 -i {0} > {1} 2> bedpe.log""".format(bam, bedpe)
    
    subprocess.check_output(cmd, shell=True)
    logger.debug('DONE: {}'.format(cmd))

    return bedpe


def bedpe_to_fragbed(bedpe, max_insert_size=500):

    diff_chroms_count = 0
    frag_lengths = []

    fragbed = os.path.splitext(bedpe)[0] + "_clean.bed"

    frag_filt_n = 0
    
    with open(fragbed, 'w') as f:
        with open(bedpe, 'r') as bedin:
            for l in bedin:
                # Check if different chromosomes 
                if l.split()[0] != l.split()[3]:
                    diff_chroms_count += 1
                    continue
            
                coords = [int(l.split()[i]) for i in [1,2,4,5]]
                #print(coords)
                start = min(coords) 
                stop = max(coords)

                frag_l = stop - start
                frag_lengths.append(frag_l)

                # This is setting the correct strand (strand of the second of pair)
                if frag_l  <= max_insert_size:
                    f.write('{0}\t{1}\t{2}\t{3}\t.\t{4}\n'.format(l.split()[0], start, stop, l.split()[6], l.split()[9]))
                else:
                    frag_filt_n += 1
        
    logger.info("Number of pairs mapping on different chromosomes: {}".format(diff_chroms_count))
    logger.info("Min frag_length: {0}; Max frag_length: {1};".format(min(frag_lengths), max(frag_lengths)))
    logger.info("Number of fragments filtered out because insert size > {0}b: {1}".format(max_insert_size, frag_filt_n))

    return fragbed


def identify_divtrans(bed):
    """From a bed file of fragments infer regions of divergent
    transcription
    """
    pass

