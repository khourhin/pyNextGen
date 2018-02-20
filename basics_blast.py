#! /usr/bin/env python3

import os
import subprocess

import click
import pandas as pd

from mylog import get_logger

logger = get_logger(__file__, __name__)


class Blasting(object):
    """Automatize the blasting workflow
    """

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return '<Blasting object query:{query}; db:{db}>'.format(query=self.query,
                                                                 db=self.db)

    def detect_db_type(self):
        """ Automatically detect the type of db ie either nucl or prot"""
        # TO DO
        pass
                
    def format_results(self):
        """Properly format the blast results

        WARNING: For now, only the outfmt=6 case is treated"""

        blast_df = pd.read_csv(self.outfile, sep='\t',
                               header=None, names=['qseqid', 'sseqid', 'pident', 'length',
                                                   'mismatch', 'gapopen', 'qstart', 'qend',
                                                   'sstart', 'send', 'evalue', 'bitscore'])
        
        return blast_df

    def format_db(self):
        """
        Format a fasta file for having them as db for blast+
        """

        # Check if formatting db files already existing
        db_exts = ['.phr', '.pin', '.psq'] if self.dbtype == 'prot' else ['.nhr', '.nin', '.nsq']
        db_exist = all([os.path.isfile(self.db + s) for s in db_exts])
        
        if db_exist:
            logger.info('Database already formatted for file: {}'.format(self.db))
            
        else:
            logger.info('Formating database: {}'.format(self.db))
            cmd = 'makeblastdb -dbtype {dbtype} -in {fasta}'.format(dbtype=self.dbtype, fasta=self.db)
            subprocess.check_output(cmd, shell=True)


    def blast(self, only_best=False, sup_args=''):
        """
        Run the blast analysis:

        TO IMPROVE

        ASSUMPTIONS:
        - The query is always nucleotides
        - So this method performs blastx if dbtype=prot and blastn if dbtype=nucl
        """

        if self.dbtype == 'prot':
        
            cmd = ('blastx -query {query} -db {db} -out {outfile} -outfmt {outfmt} -evalue {evalue} -num_threads {threads}'
                   .format(query=self.query, db=self.db, outfile=self.outfile,
                           outfmt=self.outfmt, evalue=self.evalue, threads=self.cpus))

        elif self.dbtype == 'nucl':
            cmd = ('blastn -query {query} -db {db} -out {outfile} -outfmt {outfmt} -evalue {evalue} -num_threads {threads}'
                   .format(query=self.query, db=self.db, outfile=self.outfile,
                           outfmt=self.outfmt, evalue=self.evalue, threads=self.cpus))

        if self.only_best:
            cmd = cmd + ' -max_target_seqs 1'

        logger.info('Executing: {cmd}'.format(cmd=cmd))
        
        try:
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e.output.decode())

    def run(self):

        if os.path.isfile(self.outfile):
            logger.info('Blast results already computed apparently, so using available file: {}'.format(self.outfile))
        else:
            self.format_db()
            self.blast()

        res_df = self.format_results()
        return res_df


# FOR LEGACY, BROCKEN BY ADDING BLasting self attributes for outfmt, evalue, etc ...

# @click.command()
# @click.option('--query', '-q', required=True, type=click.Path(exists=True, resolve_path=True))
# @click.option('--db', '-d', required=True, type=click.Path(exists=True, resolve_path=True))
# @click.option('--dbtype', '-t', default='nucl', type=click.Choice(['nucl', 'prot']))
# @click.option('--outfile', '-o', default='blast_out.txt', type=click.Path(resolve_path=True))
# @click.option('--cpus', '-c', default=1, type=click.INT)
# @click.option('--onlybest', default=False, flag_value=True)
# def main(**kwargs):
#     """A blast interface"""
    
#     b = Blasting(**kwargs)
#     b.run()
    
    
# if __name__ == '__main__':
#     main()
