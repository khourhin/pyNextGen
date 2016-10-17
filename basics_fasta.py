import pysam
import logging
import basics_nuc_seq as bns

logger = logging.getLogger(__name__)


def fasta_stats(fasta):
    logger.info("Producing stats for: %s" % fasta)

    seq_len = []
    GC_content = []
    N_count = 0

    with pysam.FastxFile(fasta) as fa:
        for seq in fa:
            seq_len.append(len(seq.sequence))
            GC_content.append(bns.get_seq_GC(seq.sequence))
            N_count += seq.sequence.count("N")

