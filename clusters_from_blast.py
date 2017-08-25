import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import pysam
import os
import sys
import basics_bam as bb
from mylog import get_logger

logger = get_logger(__file__, __name__)



######## BADDDDDDDDDD
#### CHANGE THE RANGE IN GET_ISOFORM CODE for it to be dynamic

def get_isoform_code(exon_list, nExons):

    exon_l = [int(x.strip("exon")) - 1 for x in exon_list]
    isoform = [1 if x in exon_l else 0 for x in range(nExons)]
    return isoform


def get_clusters_from_blast(blast_out):
    """
    From blastn results against a db of exons returns a dictionary:
    key: exon composition (in format '01011101...')
    value: list of seqids with this exon composition
    """
    
    res_dict = {}
    
    with open(blast_out) as f:
        for line in f:
            qname = line.split()[0]
            if qname in res_dict:

                res_dict[qname].add(line.split()[1])
            else:
                res_dict[qname] = set([line.split()[1]])

    #print res_dict
    res = set(frozenset(i)
              for i in res_dict.values())

    clus_dict = {}

    for i in res:
        clus_dict[i] = [x for x in res_dict if res_dict[x] == i]

    freqs = Counter([len(x) for x in clus_dict.values()])
    logger.info(freqs)

    plt.bar(freqs.keys(), np.log(freqs.values()) + 1, width=1)
    plt.xlabel("Number of read in a cluster")
    plt.ylabel("Log(Frequencies) + 1")
    plt.savefig(config_fig_clust_hist)
    plt.close()
    
    return clus_dict


def get_clusters_from_transcript_blast(blastout):
    """
    From blastn results against a db of assembled transcripts returns
    a dictionary:
    key: ref transcript id
    value: list of seqids which have a best hit against this transcript

    This assumes (it was checked first) that the first hit is always the 
    best according to bitscore
    """

    seq_dict={}

    for line in open(blastout, 'r'):
        seq_id = line.split()[0]
        hit_id = line.split()[1]

        # Takes only the first hits
        if seq_id not in seq_dict:
            seq_dict[seq_id] = hit_id

    # Inverting the dict to have clus_name -> [seqs]
    print("INverting dict")
    clus_dict = {}
    for k, v in seq_dict.iteritems():
        clus_dict[v] = clus_dict.get(v, [])
        clus_dict[v].append(k)

    logger.info('# Clusters: {}'.format(len(clus_dict)))

    return clus_dict

def get_cluster_stat(clus_dict):
    """
    Produce statistics from a cluster dict obtained by get_clusters_from_blast
    """
    
    clusters_sup1 = [clus for clus in clus_dict if len(clus_dict[clus]) > 1]
    logger.info("Number of blast clusters generated: {}".format(len(clusters_sup1)))

    
def get_fastas_from_clusters(in_fasta, clus_dict, nExons=None):

    os.makedirs(config_fasout_dir)
    
    for clus_name in clus_dict:
        clus_size = len(clus_dict[clus_name])
        
        # If clusters have been defined by exon composition (get_clusters_from_blast)
        if nExons:
            iso_codes = get_isoform_code(clus_name, nExons)
            fas_out = os.path.join(config_fasout_dir,
                                   "clus_{0}_{1}seqs.fas".format("".join(str(x)
                                                                         for x in iso_codes), clus_size))
            
        # If clusters have been defined by transcripts hit (get_clusters_from_transcript_blast)
        else:
            fas_out = os.path.join(
                config_fasout_dir,
                "clus_{0}_{1}seqs.fas".format(clus_name, clus_size))

        with open(fas_out, 'w') as fout:
            
            with pysam.FastxFile(in_fasta) as fa:
                for seq in fa:
                    if seq.name in clus_dict[clus_name]:
                        fout.write(">" + seq.name + "\n")
                        fout.write(seq.sequence + "\n")

                        
def get_bams_from_cluster(bam_in, clus_dict, nExons):

    os.makedirs(config_bamout_dir)

    for clus_name in clus_dict:
        clus_size = len(clus_dict[clus_name])
        bam_out = os.path.join(config_bamout_dir,
                               "clus_{0}_{1}seqs.bam".format("".join(str(x)
                                                                     for x in get_isoform_code(clus_name, nExons)),
                                                             clus_size))

        read_sel = bb.get_read_by_id(clus_dict[clus_name], bam_in)
        bb.print_bam(read_sel, bam_in, bam_out)


if __name__ == "__main__":
    
    # The blast output in outfmt 6
    blastout = sys.argv[1]
    # The fasta file to cluster
    fasta = sys.argv[2]
    bam = sys.argv[3]
    nExons = int(sys.argv[4])

    config_fig_clust_hist='test.png'
    config_fasout_dir= blastout + '_clusters_fas'
    config_bamout_dir= blastout + '_clusters_bam'

    # Cluster from exon blast
    clus_dict = get_clusters_from_blast(blastout)
    get_cluster_stat(clus_dict)
    get_fastas_from_clusters(fasta, clus_dict, nExons)

    # Cluster from transcript blast
    # clus_dict = get_clusters_from_transcript_blast(blastout)
    # get_fastas_from_clusters(fasta, clus_dict)
