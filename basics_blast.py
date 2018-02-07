#! /usr/bin/env python3

import sys
import os
import subprocess

import click
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
    
    def format_db(self):
        """
        Format a fasta file for having them as db for blast+
        """

        # Check if formatting db files already existing
        db_exts = ['.phr', '.pin', '.psq'] if self.dbtype=='prot' else ['.nhr', '.nin', '.nsq']
        db_exist = [os.path.isfile(self.db + s) for s in db_exts]
        
        if db_exist:
            logger.info('Database already formatted for file: {}'.format(self.db))
            
        else:
            cmd = 'makeblastdb -dbtype {dbtype} -in {fasta}'.format(dbtype=self.dbtype, fasta=self.db)    
            subprocess.check_output(cmd, shell=True)

    def blastx(self, onlyBest=False):
        """
        Run the blastx analysis
        """

        cmd = ('blastx -query {query} -db {db} -out {outfile} -outfmt {outfmt} -evalue {evalue} -num_threads {threads}'
               .format(query=self.query, db=self.db, outfile=self.outfile,
                       outfmt='6', evalue='1e-6', threads=self.cpus))

        if self.onlybest:
            cmd = cmd + ' -max_target_seqs 1'

        logger.info('Executing: {cmd}'.format(cmd=cmd))
        subprocess.check_output(cmd, shell=True)
        
            
    def run(self):

        self.format_db()
        self.blastx()

@click.command()
@click.option('--query', '-q', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('--db', '-d', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('--dbtype', '-t', default='nucl', type=click.Choice(['nucl', 'prot']))
@click.option('--outfile', '-o', default='blast_out.txt', type=click.Path(resolve_path=True))
@click.option('--cpus', '-c', default=1, type=click.INT)
@click.option('--onlybest', default=False, flag_value=True)
def main(**kwargs):
    """A blast interface"""
    
    b = Blasting(**kwargs)
    b.run()
    
    
if __name__ == '__main__':
    main()

################################################################################
# LEGACY

# ASSUMPTIONS:
# Blast results are in -outfmt 6 format

# blast_out = "/home/ekornobis/analysis/allemand/v4/blast5_sup27.tab"

# # Plotting blast results
# blast_res = np.genfromtxt(blast_out, dtype=None,
#                   names="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore")


# res_l = blast_res['length']
# plt.hist(res_l, range(0, max(res_l), 10))
# plt.show()


# def do_blastP(query, db, outfile, onlyBest=False):
#     """
#     Launched a blastp analysis for the [query] against the [db]
#     """
#     cmd = [PATH_BLAST + 'blastp', '-query', query, '-db', db, '-out', outfile,
#            '-outfmt', '6', '-evalue', '1e-6']
#     if onlyBest:
#         cmd = cmd + [ '-max_target_seqs', '1' ]
#     proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
#     out, err = proc.communicate()    
