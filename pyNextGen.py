#! /usr/bin/python

import basics_fastq as bfq
import sys

# Checked 'click' library for parameters handling. Sounds super cool
# but visibly cannot deal with multiple arguments (max number not
# predefined) without having to specify a - flag each time. ie cannot
# do: "pyNextGen data/*"

# TODO Use threading


def fastq_stats(fastqs):
    bfq.print_all_stats(fastqs)

if __name__ == "__main__":
    fastq_stats(sys.argv[1:])
