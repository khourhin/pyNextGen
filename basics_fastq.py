import logging
import numpy
import pysam
import collections
import os.path
import basics_nuc_seq as bns
from multiprocessing import Pool

logger = logging.getLogger(__name__)


def fastq_stats(fastq):
    """
    Get stats from a fastq (or gz compressed fastq) file

    Return a dictionnary of stats
    """

    rds_len = []
    GC_content = []
    N_count = 0

    logger.info("Producing stats for: %s" % fastq)

    with pysam.FastxFile(fastq) as fq:
        for read in fq:
            rds_len.append(len(read.sequence))
            GC_content.append(bns.get_seq_GC(read.sequence))
            N_count += read.sequence.count("N")

    res_dict = collections.OrderedDict([
        ("FileName", os.path.basename(fastq)),
        ("Nreads(M)", len(rds_len) / 1.0e6),
        ("Nbases(G)", sum(rds_len) / 1.0e9),
        ("Ns", N_count),
        ("MinLen", min(rds_len)),
        ("MaxLen", max(rds_len)),
        ("MeanLen", numpy.mean(rds_len)),
        ("StdevLen", numpy.std(rds_len)),
        ("MeanGC", numpy.mean(GC_content)),
    ])
    return res_dict


def print_all_stats(fastqs_list, threads=10):
    """
    Iterate fastq_stats over a list of fastq files.

    Print a csv
    """

    res_list = []

    p = Pool(threads)
    res_list = p.map(fastq_stats, [fq for fq in fastqs_list])
        
    # header
    print(",".join(res_list[0].keys()))

    # Rows
    for res in res_list:
        print(",".join([str(x) for x in res.values()]))


    # header
#    print(",".join(res_list[0].keys()))

    # Rows
    # for res in res_list:
    #     print(",".join([str(x) for x in res.values()]))
