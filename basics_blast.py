import matplotlib.pyplot as plt
import numpy as np

# This is a dirty stub !

# ASSUMPTIONS:
# Blast results are in -outfmt 6 format

blast_out = "/home/ekornobis/analysis/allemand/v4/blast5_sup27.tab"

# Plotting blast results
blast_res = np.genfromtxt(blast_out, dtype=None,
                  names="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore")


res_l = blast_res['length']
plt.hist(res_l, range(0, max(res_l), 10))
plt.show()
