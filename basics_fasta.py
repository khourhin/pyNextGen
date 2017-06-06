import argparse
import pysam
import logging
import basics_nuc_seq as bns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# TODO: check if dedup_fasta can be improved using yield as in
# /home/ekornobis/code/allemand/gphn/7.2_generate_ref_before_rsem.py

def dedup_fasta(fasta, fasta_out=None, graph=True):
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

    if graph:
        labels, values = zip(*reads_per_seq.items())
        index = np.arange(len(labels))
        plt.bar(index, [np.log(x) for x in values])
        plt.xlabel("Number of copies")
        plt.ylabel("Log(Frequencies) + 1")
        plt.savefig("Number_of_copies_per_unique_read.png")
        plt.show()
    
    # Create a fasta file with only unique sequences
    if fasta_out:
        with open(fasta_out, 'w') as f:
            for seq in dedup_dict:
                f.write('>' + '-'.join(dedup_dict[seq]) + '\n')
                f.write(seq + '\n')

    return dedup_dict


def fasta_stats(fasta):
    log.info("Producing stats for: %s" % fasta)

    seq_len = []
    GC_content = []
    N_count = 0

    with pysam.FastxFile(fasta) as fa:
        for seq in fa:
            seq_len.append(len(seq.sequence))
            GC_content.append(bns.get_seq_GC(seq.sequence))
            N_count += seq.sequence.count("N")

    plt.hist(seq_len)
    plt.show()

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('fastaIn', help='')
    parser.add_argument('fastaOut', help='')
    parser.add_argument('--noGraph', help='', action='store_false')

    args = parser.parse_args()

#    fasta_stats(args.fas)
    dedup_fasta(args.fastaIn, args.fastaOut, args.noGraph)
    
