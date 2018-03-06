import os
from mylog import get_logger

logger = get_logger(__file__, __name__)

def apply_threads(func, arg_list, nthreads=10):
    p = Pool(nthreads)
    p.starmap(func, arg_list)


def simplify_path(filename, prefix='', suffix='', keep_path=False, check_exist=True):
    """
    Produce a name for an output file
    """

    if not keep_path:
        output_name = os.path.basename(filename)
    else:
        output_name = filename

    # Remove all extensions
    ext = None
    while ext != '':
        output_name, ext = os.path.splitext(output_name)

    output_name = prefix + output_name + suffix

    if check_exist and os.path.exists(output_name):
        logger.warning('File already exists: {}'.format(output_name))
    
    return output_name

