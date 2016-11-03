import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import pysam
import os
import logging
import sys
import basics_bam as bb

######## BADDDDDDDDDD
#### CHANGE THE RANGE IN GET_ISOFORM CODE for it to be dynamic

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

config_fig_clust_hist='test.png'
config_fasout_dir='clusters_fas'
config_bamout_dir='clusters_bam'

def get_isoform_code(exon_list, nExons):

    exon_l = [int(x.strip("exon")) - 1 for x in exon_list]
    isoform = [1 if x in exon_l else 0 for x in range(nExons)]
    return isoform


def get_clusters_from_blast(blast_out):
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
    log.info(freqs)

    plt.bar(freqs.keys(), np.log(freqs.values()) + 1, width=1)
    plt.xlabel("Number of read in a cluster")
    plt.ylabel("Log(Frequencies) + 1")
    plt.savefig(config_fig_clust_hist)
    plt.close()
    
    return clus_dict


def get_cluster_stat(clus_dict):
    clusters_sup1 = [clus for clus in clus_dict if len(clus_dict[clus]) > 1]
    print len(clusters_sup1)

    
def get_fastas_from_clusters(in_fasta, clus_dict, nExons):

    os.makedirs(config_fasout_dir)
    
    for clus_name in clus_dict:
        clus_size = len(clus_dict[clus_name])
        fas_out = os.path.join(config_fasout_dir,
                               "clus_{0}_{1}seqs.fas".format("".join(str(x)
                                                                     for x in get_isoform_code(clus_name)),
                                                             clus_size))

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
                                                                     for x in get_isoform_code(clus_name)),
                                                             clus_size))

        read_sel = bb.get_read_by_id(clus_dict[clus_name], bam_in)
        bb.print_bam(read_sel, bam_in, bam_out)

