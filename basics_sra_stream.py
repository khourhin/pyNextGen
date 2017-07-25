# NOT WORKING, VISIBLY MISSING A NGS LIBRARY

from srastream import SraReader


sra='SRR019721'


with SraReader(sra) as reader:
    for frags in reader:
        print('\n'.join(str(read) for read in reads))
