
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/ekornobis/programs/anaconda3/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X\x1e\x00\x00\x00star/SRR019721/Aligned.out.bamq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x07\x00\x00\x00threadsq\nK\x01X\t\x00\x00\x00wildcardsq\x0bcsnakemake.io\nWildcards\nq\x0c)\x81q\rX\t\x00\x00\x00SRR019721q\x0ea}q\x0f(h\x08}q\x10X\x06\x00\x00\x00sampleq\x11K\x00N\x86q\x12sX\x06\x00\x00\x00sampleq\x13h\x0eubX\x05\x00\x00\x00inputq\x14csnakemake.io\nInputFiles\nq\x15)\x81q\x16X\x15\x00\x00\x00reads/SRR019721.fq.gzq\x17a}q\x18(h\x08}q\x19h\x13K\x00N\x86q\x1ash\x13h\x17ubX\x06\x00\x00\x00paramsq\x1bcsnakemake.io\nParams\nq\x1c)\x81q\x1d(X\x05\x00\x00\x00indexq\x1eX\x00\x00\x00\x00q\x1fe}q (h\x1eh\x1eh\x08}q!(h\x1eK\x00N\x86q"X\x05\x00\x00\x00extraq#K\x01N\x86q$uh#h\x1fubX\x06\x00\x00\x00configq%}q&X\x03\x00\x00\x00logq\'csnakemake.io\nLog\nq()\x81q)X\x17\x00\x00\x00logs/star/SRR019721.logq*a}q+h\x08}q,sbX\t\x00\x00\x00resourcesq-csnakemake.io\nResources\nq.)\x81q/(K\x01K\x01e}q0(h\x08}q1(X\x06\x00\x00\x00_coresq2K\x00N\x86q3X\x06\x00\x00\x00_nodesq4K\x01N\x86q5uh4K\x01h2K\x01ubX\x04\x00\x00\x00ruleq6X\x0b\x00\x00\x00star_singleq7ub.')
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

