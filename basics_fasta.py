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

# TO IMPROVE:
# Fasta method to_dict


class Fasta(object):
    """A Fasta sequence file object
    """
    def __init__(self, path):

        if os.path.isfile(path):
            self.path = path
        else:
            raise IOError('No such file: {}'.format(path))

    def __repr__(self):
        return '<Fasta Object for: {}'.format(self.path)

    def to_dict(self, simple_ID=True):
        """
        Get a dictionnary from fasta file with:
        key = seq_id
        value = seq

        If simple_ID is True, then the seq_id in the dictionnary will be
        the original one truncated after the first white space.
        """
        with open(self.path, "r") as f:
            seq_d = {}
            
            # split each time a ">" in encountered
            fasta = f.read().split('>')[1:]
            
            # split the strings in 3 items tuple: seq id, the sep, and the seq
            tmp = [seq.partition("\n") for seq in fasta]
            
            # build a dictionnary with key=seq id, value=seq
            for i in range(len(tmp)):
            
                if simple_ID:
                    # for ex: for trinity output, split()[0] removes the len, path infos
                    seq_d[tmp[i][0].split()[0]] = tmp[i][2].replace("\n", "").upper()                
                else:
                    seq_d[tmp[i][0]] = tmp[i][2].replace("\n", "").upper()

        # Check if no empty sequences
        empty_seqs = [k for k in seq_d if len(seq_d[k]) == 0]
        if empty_seqs:
            raise IOError("Something seems wrong with this sequence:\n" 
                          + "\n".join(empty_seqs))
        
        return seq_d
    
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
        fasta_df = fasta_df.set_index('seq_names')

        return fasta_df

    def summary(self):
        """
        Produce a statistic summary for the Fasta file
        """
        
        fasta_df = self.get_stats()
        plt.subplot(1,2,1)
        plt.boxplot(fasta_df['seq_len'])
        plt.title('Sequence length boxplot')
        plt.subplot(1,2,2)
        plt.hist(fasta_df['seq_len'], bins=100, log=True)
        plt.title('Sequence length Histogram(log)')
        plt.show()

        print(fasta_df.describe())

        return fasta_df


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
    
