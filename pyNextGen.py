#! /usr/bin/python

import logging
import basics_fastq as bfq
import sys
import clusters_from_blast as cfb


# Checked 'click' library for parameters handling. Sounds super cool
# but visibly cannot deal with multiple arguments (max number not
# predefined) without having to specify a - flag each time. ie cannot
# do: "pyNextGen data/*"

# TODO Use threading

def fastq_stats(fastqs):
    bfq.print_all_stats(fastqs)


def blast_clustering(blastout, fasta, bam, nExons):

    clus_dict = cfb.get_clusters_from_blast(blastout)
    cfb.get_cluster_stat(clus_dict)
    cfb.get_fastas_from_clusters(fasta, clus_dict, nExons)
#    cfb.get_bams_from_cluster(bam, clus_dict, nExons)

if __name__ == "__main__":

    # logging.basicConfig(level=logging.ERROR)
    logging.basicConfig(level=logging.INFO)

    # Fastq Stats
    # fastq_stats(sys.argv[1:])

    # The blast output in outfmt 6
    blastout = sys.argv[1]
    # The fasta file to cluster
    fasta = sys.argv[2]
    bam = sys.argv[3]
    nExons = int(sys.argv[4])
    blast_clustering(blastout, fasta, bam, nExons)
