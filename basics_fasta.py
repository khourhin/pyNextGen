#! /usr/bin/env python

import argparse
import pysam
import logging
import basics_nuc_seq as bns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import click
import os
import pandas as pd

logger = logging.getLogger(os.path.basename(__file__) + " - " +  __name__)
logger.setLevel(logging.DEBUG)

handler = logging.FileHandler('/home/ekornobis/logs/common.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(handler)


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
    logger.info('Number of copies per reads:')
    logger.info(reads_per_seq)

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

@click.command()
@click.argument('fasta')
def fasta_stats(fasta):
    logger.info("Producing stats for: %s" % fasta)

    fasta_df = pd.DataFrame({
        'seq_names': [seq.name for seq in pysam.FastxFile(fasta)],
        'seq_len': [len(seq.sequence) for seq in pysam.FastxFile(fasta)],
        'GC_content': [bns.get_seq_GC(seq.sequence) for seq in pysam.FastxFile(fasta)],
        'Ns': [seq.sequence.count('N') for seq in pysam.FastxFile(fasta)],
    })

    print(fasta_df)
    
    plt.hist(fasta_df['seq_len'])
    plt.show()

    
if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.add_argument('fastaIn', help='')
    # parser.add_argument('fastaOut', help='')
    # parser.add_argument('--noGraph', help='', action='store_false')

    # args = parser.parse_args()

    # fasta_stats(args.fas)
#    dedup_fasta(args.fastaIn, args.fastaOut, args.noGraph)
    fasta_stats()
    
