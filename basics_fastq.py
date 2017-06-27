import logging
import numpy
import pysam
import collections
import os.path
import basics_nuc_seq as bns
from multiprocessing import Pool
import argparse
import os

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

handler = logging.FileHandler('basics_fastq.log')
handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(handler)


def iterate_fastq(fastq):
    """ Iterate through a fastq file """

    with pysam.FastxFile(fastq) as fq:
        for read in fq:
            yield read


def filter_fastq(fastq, filter_fun):
    """ Filter a fastq file according to the specified filter """

    filtered_reads = set()

    for read in iterate_fastq(fastq):
        if filter_fun(read):
            filtered_reads.add(read.name)
            log.debug(read)

    log.info('Done filtering {}'.format(fastq))
    return filtered_reads


def has_polyA(read, n=6):
    """
    SPECIFIC !!! For now only check if the read ends with n "A"
    Intended for stranded specific data
    """
    if read.sequence.endswith("A" * n):
        return True


def has_polyT(read, n=6):
    """
    SPECIFIC !!! For now only check if the read starts with n "T"
    Intended for stranded specific data
    """
    if read.sequence.startswith("T" * n):
        return True


def print_filter_fastq(fastq, read_set):
    """
    From a fastq file and a list of read name, produce a filtered
    fastq file with only the selected reads in.
    Output fastq will be with the suffixe '_filtered.fq'
    """

    fastq_out = fastq.split(os.extsep)[0] + "_filtered.fq"

    with open(fastq_out, 'w') as fout:
        for read in iterate_fastq(fastq):
            if read.name in read_set:
                fout.write(str(read) + '\n')

    log.info("Written filtered fastq to {}".format(fastq_out))


def fastq_stats(fastq):
    """
    Get stats from a fastq (or gz compressed fastq) file

    Return a dictionnary of stats
    """

    rds_len = []
    GC_content = []
    N_count = 0

    logger.info("Producing stats for: %s" % fastq)

    for read in iterate_fastq(fastq):
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


def filter_polyA_reverse_stranded_reads(fq1, fq2):
    """
    Produce fastqs file with only reads containing (in either of the
    pair) a polyA tail. As it is IT WILL ONLY WORK WITH REVERSE
    STRANDED PAIRED DATA.
    """

    log.info("Starting polyA selection with pair: {0}; {1}".format(fq1, fq2))
    read_set1 = filter_fastq(fq1, has_polyT)
    log.info('Found {0} polyT hits in {1}'.format(len(read_set1), fq1))
    read_set2 = filter_fastq(fq2, has_polyA)
    log.info('Found {0} polyA hits in {1}'.format(len(read_set2), fq2))
    read_set = read_set1 | read_set2

    print_filter_fastq(fq1, read_set)
    print_filter_fastq(fq2, read_set)


def apply_threads(func, arg_list, nthreads=10):
    p = Pool(nthreads)
    p.starmap(func, arg_list)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--stats', '-s', nargs='+',
                        help='At least 1 fastq file to produce stats from')
    parser.add_argument('--fastq1', '-f1', nargs='+',
                        help='FOR POLYA_ONLY_CHECK_CODE. A fastq file. First of pair. In the path, Dots "." should only separate extensions !')
    parser.add_argument('--fastq2', '-f2', nargs='+',
                        help='FOR POLYA_ONLY_CHECK_CODE. A fastq file. Second of pair. In the path, Dots "." should only separate extensions !')
    parser.add_argument('--threads', '-t', help='Number of threads to use.')
    args = parser.parse_args()

    if args.stats:
        print_all_stats(args.stats)

    elif len(args.fastq1) != len(args.fastq2):
        raise IOError(
            "The number of files specified with --fastq1 and --fastq2 should be equal.")

    else:
        apply_threads(filter_polyA_reverse_stranded_reads,
                      zip(args.fastq1, args.fastq2))

    
