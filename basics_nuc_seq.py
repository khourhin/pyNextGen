#-------------------------------------------------------------------------------
def get_seq_GC(seq):
    """
    From a nucleotide sequence return the GC content
    """
    seq = seq.lower()
    n_C = seq.count("c")
    n_G = seq.count("g")
    
    GC_content = ( n_C + n_G ) / float(len(seq))
    return GC_content


def get_bases_composition(seq):
    """
    Proportion of A,T,G,C bases in a sequence
    """
    seq = seq.lower()
    n_A = seq.count("a")
    n_T = seq.count("t")
    n_G = seq.count("g")
    n_C = seq.count("c")

    return {'A': n_A, 'T': n_T, 'G': n_G, 'C': n_C}
        
#-------------------------------------------------------------------------------
def get_all_GCs(seq_d):
    """
    From a seq_dict (fas2dict)
    Return all GCs as dictionnary with seqid as keys
    """
    GCs_d = {k: get_seq_GC( seq_d[k] ) for k in seq_d }
    
    return GCs_d


def get_all_bases_compositions(seq_d):
    """
    From a seq_dict, return all the bases compositions
    """

    bases_comp_d = {k: get_bases_composition( seq_d[k] ) for k in seq_d }

    return bases_comp_d

#-------------------------------------------------------------------------------
def get_all_lens(seq_d):
    """
    From a seq_dict (fas2dict):
    Return all lengths as dictionnary with seqid as keys
    """
    lens_d = {k: len( seq_d[k] ) for k in seq_d }

    return lens_d

#-------------------------------------------------------------------------------
def get_N50(seq_d):
    """
    From a seq_dict (fas2dict):
    Return the N50
    """
    import numpy 

    lens_d = get_all_lens(seq_d)
    N50_l = []

    for i in lens_d.values():
        for j in range(i):
            N50_l.append(i)

    N50 = numpy.median(N50_l)
    return N50

#-------------------------------------------------------------------------------

    
