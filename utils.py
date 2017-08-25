import os
from mylog import get_logger

logger = get_logger(__file__, __name__)

def apply_threads(func, arg_list, nthreads=10):
    p = Pool(nthreads)
    p.starmap(func, arg_list)


def simplify_path(filename, prefix='', suffix=''):
    """
    Produce a name for an output file
    """
    
    output_name = os.path.basename(filename)
    output_name = os.path.splitext(output_name)[0]
    output_name = prefix + output_name + suffix

    return output_name
