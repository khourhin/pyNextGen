import sys
import os
import pysam
import logging
import numpy as np
import subprocess
from collections import Counter
from itertools import chain
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def dedup_fasta(fasta, fasta_out=None):
    """
    From a fasta file return a dictionary with the duplicate sequences
    merged under a single entry (the name of the sequence will be the
    concatenation of all the previous sequence names).
    WARNING: Output dictionnary is as:
    keys=sequences
    values=sequence names

    Produce a fasta outfile if specified
    """

    dedup_dict = {}
    for seq in pysam.FastxFile(fasta):
        if seq.sequence in dedup_dict:
            dedup_dict[seq.sequence].append(seq.name)
        else:
            dedup_dict[seq.sequence] = [seq.name]

    reads_per_seq = Counter(map(lambda x: len(x), dedup_dict.values()))
    log.info('Number of copies per reads:')
    log.info(reads_per_seq)

    # Plot log(x) + 1 to see count =1 and large values
    plt.bar(reads_per_seq.keys(), np.log(reads_per_seq.values()) + 1)
    plt.xlabel("Number of copies")
    plt.ylabel("Log(Frequencies) + 1")
    plt.savefig("Number_of_copies_per_unique_read.png")
    plt.close()

    # Create a fasta file with only unique sequences
    if fasta_out:
        with open(fasta_out, 'w') as f:
            for seq in dedup_dict:
                f.write('>' + '-'.join(dedup_dict[seq]) + '\n')
                f.write(seq + '\n')

    return dedup_dict


def get_putative_exons(bam):
    """
    Get a bed file with putative exons from a bam file
    """

    # Index the bam file
    pysam.index(bam)
    
    cmd = 'bedtools bamtobed -split -i {0} | bedtools sort -i - | bedtools merge -i - > {1}'.format(bam, config_potential_exon)
    subprocess.call(cmd, shell=True)


def get_overall_exon_counts(bam):
    """
    Get the total number of reads which are supporting each exons
    """

    cmd = 'bedtools multicov -split -bams {0} -bed {1}'.format(bam, config_potential_exon)
    out = subprocess.check_output(cmd, shell=True)
    out = out.strip().split('\n')
    out = [x.split('\t') for x in out]
    log.info("Number of exons: {}".format(len(out)))
    log.info("Exon coverage distribution: {}".format(out))

    # xs = map(int, chain.from_iterable([i[1:3] for i in out]))
    # ys = map(int, chain.from_iterable([[i[3], i[3]] for i in out]))

    # plt.plot(xs, ys)
    # plt.show()


def get_exon_count_per_read(bam):
    """
    From a bam file with all reads aligned to the genome, generate one
    bam file per read and counts per exons This might be unecessary
    (this is done to produce the count matrix for clustering)
    """

    sam = pysam.AlignmentFile(bam, 'rb')
    fcount = 0
    count_dict = {}

    for rd in sam.fetch(until_eof=True):
        bam_path = os.path.join(config_bam_out, str(fcount) + '.bam')

        with pysam.AlignmentFile(bam_path, 'wb', template=sam) as bam_out:
            bam_out.write(rd)

        fcount += 1

        log.debug("Indexing {}".format(bam))
        pysam.index(bam_path)

        cmd = 'bedtools multicov -split -bams {0} -bed {1}'.format(bam_path, config_potential_exon)
        out = subprocess.check_output(cmd, shell=True).strip('\n')

        if rd in count_dict:
            raise IOError('Duplicate read names detected')

        count_dict[rd.query_name] = ''.join([x.split('\t')[3]
                                             for x in out.strip().split('\n')])

        with open(config_count_out, 'a') as f:
            f.write(rd.query_name + '\t')
            f.write('\t'.join(list(count_dict[rd.query_name])) + '\n')

    log.info('Total sequences put in single indexed bams: {}'.format(fcount))
    return count_dict


def cluster_by_exon_composition(bam, count_dict):
    """
    From a bam file with all reads aligned to the genome and a
    dictionnary (k: read name, v: exon composition), get bams grouped
    by exon composition
    """

    exon_compo = set(count_dict.values())

    for ex in exon_compo:

        seq_in_clust = [k for k in count_dict if count_dict[k] == ex]
        log.info('For cluster {0}: {1} sequence(s)'.format(ex, len(seq_in_clust)))

        bam_path = os.path.join(config_bam_clusters, 'cluster' + ex + '.bam')
        fas_path = os.path.join(config_fas_clusters, 'cluster' + ex + '.fas')        
        sam = pysam.AlignmentFile(bam, 'rb')

        with pysam.AlignmentFile(bam_path, 'wb', template=sam) as bam_out:

            for rd in sam.fetch(until_eof=True):
                if rd.query_name in seq_in_clust:
                    bam_out.write(rd)

                    with open(fas_path, 'a') as fas:
                        fas.write('>' + rd.query_name + '\n')
                        fas.write(rd.seq + '\n')

        pysam.index(bam_path)


def prepare_outfolder():
    log.info('Creating output folders in {}'.format(config_outfolder))
    os.makedirs(config_outfolder)
    os.makedirs(config_bam_clusters)
    os.makedirs(config_fas_clusters)
    os.makedirs(config_bam_out)


if __name__ == '__main__':
    # Deduplicate a fasta file
    # dedup_fasta(sys.argv[1], sys.argv[2])    

    config_bam = sys.argv[1]
    config_outfolder = sys.argv[2]

    config_potential_exon = os.path.join(config_outfolder, 'potential_exons.bed')
    config_bam_out = os.path.join(config_outfolder, 'bam_tmp')
    config_bam_clusters = os.path.join(config_outfolder, 'bam_clusters')
    config_fas_clusters = os.path.join(config_outfolder, 'fas_clusters')
    config_count_out = os.path.join(config_outfolder, 'counts')
    
    prepare_outfolder()
    get_putative_exons(config_bam)
    get_overall_exon_counts(config_bam)
    count_dict = get_exon_count_per_read(config_bam)
    cluster_by_exon_composition(config_bam, count_dict)
