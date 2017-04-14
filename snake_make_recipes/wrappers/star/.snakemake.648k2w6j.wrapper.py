
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/ekornobis/programs/anaconda3/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x01X\x06\x00\x00\x00configq\x04}q\x05X\x04\x00\x00\x00ruleq\x06X\x0b\x00\x00\x00star_singleq\x07X\x06\x00\x00\x00outputq\x08csnakemake.io\nOutputFiles\nq\t)\x81q\nX\x1e\x00\x00\x00star/SRR019721/Aligned.out.bamq\x0ba}q\x0cX\x06\x00\x00\x00_namesq\r}q\x0esbX\x05\x00\x00\x00inputq\x0fcsnakemake.io\nInputFiles\nq\x10)\x81q\x11X\x15\x00\x00\x00reads/SRR019721.fq.gzq\x12a}q\x13(h\r}q\x14X\x06\x00\x00\x00sampleq\x15K\x00K\x01\x86q\x16sh\x15csnakemake.io\nNamedlist\nq\x17)\x81q\x18h\x12a}q\x19h\r}q\x1asbubX\x03\x00\x00\x00logq\x1bcsnakemake.io\nLog\nq\x1c)\x81q\x1dX\x17\x00\x00\x00logs/star/SRR019721.logq\x1ea}q\x1fh\r}q sbX\t\x00\x00\x00resourcesq!csnakemake.io\nResources\nq")\x81q#(K\x01K\x01e}q$(X\x06\x00\x00\x00_coresq%K\x01h\r}q&(X\x06\x00\x00\x00_nodesq\'K\x00N\x86q(h%K\x01N\x86q)uh\'K\x01ubX\x06\x00\x00\x00paramsq*csnakemake.io\nParams\nq+)\x81q,(X\x00\x00\x00\x00q-X\x05\x00\x00\x00indexq.e}q/(h.h.h\r}q0(X\x05\x00\x00\x00extraq1K\x00N\x86q2h.K\x01N\x86q3uh1h-ubX\t\x00\x00\x00wildcardsq4csnakemake.io\nWildcards\nq5)\x81q6X\t\x00\x00\x00SRR019721q7a}q8(h\r}q9X\x06\x00\x00\x00sampleq:K\x00N\x86q;sh\x15h7ubub.')
######## Original script #########
__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.sample)
assert n == 1 or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

if snakemake.input.sample[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""


outprefix = os.path.dirname(snakemake.output[0]) + "/"


shell(
    "STAR "
    "{snakemake.params.extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.index} "
    "--readFilesIn {snakemake.input.sample} "
    "{readcmd} "
    "--outSAMtype BAM Unsorted "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "{log}")

