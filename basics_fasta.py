#! /usr/bin/env python

import pysam
import basics_nuc_seq as bns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import click
import os
import pandas as pd
from mylog import get_logger

logger = get_logger(__file__, __name__)

# TODO: check if dedup_fasta can be improved using yield as in
# /home/ekornobis/code/allemand/gphn/7.2_generate_ref_before_rsem.py


class Fasta(object):
    """A Fasta sequence file object
    """
    def __init__(self, path):
        self.path = path

    def __repr__(self):
        return '<Fasta Object for: {}'.format(self.path)

    def get_stats(self):
        """
        Produce statistics on the sequences in the fasta file
        """
        
        logger.info("Producing stats for: %s" % self.path)

        fasta_df = pd.DataFrame({
            'seq_names': [seq.name for seq in pysam.FastxFile(self.path)],
            'seq_len': [len(seq.sequence) for seq in pysam.FastxFile(self.path)],
            'GC_content': [bns.get_seq_GC(seq.sequence) for seq in pysam.FastxFile(self.path)],
            'Ns': [seq.sequence.count('N') for seq in pysam.FastxFile(self.path)],
        })

        return fasta_df

    def summary(self):
        """
        Produce a statistic summary for the Fasta file
        """
        
        fasta_df = self.get_stats()
        plt.hist(fasta_df['seq_len'], bins=200)
        plt.title('Sequence length frequencies')
        plt.show()


class TrinityFasta(Fasta):
    """
    A Fasta file out of Trinity de-novo assembler
    """
    
    def get_stats(self):
        fasta_df = super().get_stats()


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
    logger.info(reads_per_seq)
    logger.info('Number of copies per reads:')

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
def main(fasta):
    
    f = Fasta(fasta)

    f.summary()

    
if __name__ == "__main__":

    main()
    
