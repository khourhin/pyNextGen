complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def rev_complement(seq):

    return ''.join([complements[x] for x in reversed(seq)])

