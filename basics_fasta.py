import sys
import pysam
import logging
import basics_nuc_seq as bns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

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

    plt.hist(seq_len, bins=range(0,300,10))
    plt.show()


if __name__ == "__main__":
    fas = sys.argv[1]

    #fasta_stats(fas)
    dedup_fasta(fas, 'dedup.fas')
