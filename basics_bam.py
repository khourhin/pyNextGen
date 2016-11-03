import pysam

# IN DVPT

# get_read_by_id is HIGHLY SUBOPTIMAL (found some index methods in
# pysam but did not manage to make it work to accelerate the acces to
# read by names)


def get_read_by_id(id_list, bam):
    """
    Get the alignement from a bam file given an read id_list
    """

    read_sel = []
    with pysam.AlignmentFile(bam, 'r') as bam_in:
        for rd in bam_in.fetch(until_eof=True):
            if rd.query_name in id_list:
                read_sel.append(rd)

    return read_sel


def print_bam(read_sel, bam_in, bam_out):
    """
    From a list of pysam reads, write out a bam file. bam_in is used
    only for its header.
    """

    with pysam.AlignmentFile(bam_out, 'wb',
                             template=pysam.AlignmentFile(bam_in)) as bam_out:
        for rd in read_sel:
            bam_out.write(rd)
