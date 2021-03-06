#! /usr/bin/env python3

# table input example
# | CR_Inp_1.bam     | CR_Inp     |
# | CR_Inp_2.bam     | CR_Inp     |
# | CR_Inp_3.bam     | CR_Inp     |
# | CR_Inp_nas_1.bam | CR_Inp_nas |
# | CR_Inp_nas_2.bam | CR_Inp_nas |
# | CR_Inp_nas_3.bam | CR_Inp_nas |

# TO ADD
# MAKE BAM INDEX FIRST (ALWAYS FORGETTING ABOUT THOSE

import os
import argparse
import subprocess
import itertools
from mylog import get_logger
import logging

logger = get_logger(__file__, __name__, log_level=logging.DEBUG)


def parse_input(input_tab):
    """Parse input for majiq. Input should have:
    $1: file paths
    $2: group name
    """
    inp = {}
    groups = set()
    with open(input_tab, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
        
            file_path = line.split()[0]
            group = line.split()[1]
            name = os.path.splitext(os.path.basename(file_path))[0]

            if name in inp:
                raise IOError('Duplicates file names (first column in {})'.format(input_tab))

            inp[name]=(file_path, group)
            groups.add(group)
            
    return (inp, groups)

def create_conf(conf_file, inp):

    with open(conf_file, 'w') as f:
        f.write("""[info]
readlen={0}
samdir={1}
genome={2}
genome_path={3}
{4}
[experiments]\n""".format(args.read_len, args.bam_dir, os.path.basename(args.gfas), args.gfas, args.stranded))

        for grp in groups:
            f.write('{0}={1}\n'.format(grp, ','.join([x
                                                      for x in inp
                                                      if inp[x][1] == grp])))

def main(args, inp, groups, control, debug=False):

    if not debug:
        os.makedirs(args.outdir)
        os.chdir(args.outdir)

    logger.debug("Summary Bams/groups: {}".format(inp))
    logger.debug("Groups used: {}".format(groups))
    logger.debug("Strandness switch used: {}".format(args.stranded == "type=strand-specific"))
    
    if not debug:
        create_conf(args.majiq_conf, inp)
    
    cmd = "majiq build {0} -conf {1} --nthreads {2} --output .".format(args.gff, args.majiq_conf, args.threads)
    logger.info(cmd)

    if not debug:
        subprocess.check_output(cmd, shell=True)

    for grp in groups:

        
        dot_majiq = ' '.join([x + '.majiq.hdf5' for x in inp if inp[x][1] == grp])
        
        cmd =  'majiq psi {0} --nthreads {1} --output psi_{2} --name {2}'.format(dot_majiq, args.threads, grp)
        logger.info(cmd)
        
        if not debug:
            subprocess.check_output(cmd, shell=True)

        cmd = 'voila psi psi_{0}/{0}.deltapsi.voila -splice-graph splicegraph.hdf5 -o voila_{0}'.format(grp)
        logger.info(cmd)

        if not debug:
            # Apparently this step might not be necessary
            #subprocess.check_output(cmd, shell=True)
            pass


    groups.remove(control)
    for i in groups:
        j = control
        logger.debug('Comparing "{0}" VS "{1}" (Control will be "{1}" then ;)'.format(i, j))
        
        dot_majiq_i = ' '.join([x + '.majiq.hdf5' for x in inp if inp[x][1] == i])
        dot_majiq_j = ' '.join([x + '.majiq.hdf5' for x in inp if inp[x][1] == j])

        cmd = 'majiq deltapsi -grp1 {majiq_i} -grp2 {majiq_j} --names {i} {j} --nthreads {t} --output dpsi_{i}_VS_{j}'.format(majiq_i=dot_majiq_i, majiq_j=dot_majiq_j, i=i, j=j, t=args.threads)
        logger.info(cmd)

        if not debug:
            subprocess.check_output(cmd, shell=True)
            
        cmd = 'voila deltapsi dpsi_{i}_VS_{j}/{i}_{j}.deltapsi.voila --splice-graph splicegraph.hdf5 -o voila_{i}_VS_{j}'.format(i=i, j=j)

        logger.info(cmd)
        if not debug:
            subprocess.check_output(cmd, shell=True)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta','-m', required=True, help='A table with meta information on the bams ie $1: File names; $2: groups')
    parser.add_argument('--bam_dir','-b', required=True, help='The absolute path to the bam directory')
    parser.add_argument('--gfas','-g', required=True, help='The path to the genome fasta file')
    parser.add_argument('--gff', '-a', required=True, help='Annotation GFF file')
    parser.add_argument('--read_len', '-l', required=True, help='Length of the reads (set to maximum read length in your dataset)')
    parser.add_argument('--stranded', '-s', action='store_const', const='type=strand-specific', default='',
                        help='Specify if reads are stranded. Default is unstranded.')
    parser.add_argument('--threads', '-t', help='Number of threads to use', default=5, type=int)    
    parser.add_argument('--outdir', '-o', help='Path to output directory (NO SLASH AT THE END FOR NOW!)', default='majiq_out', type=str)    
    parser.add_argument('--majiq_conf', '-c', help='Path to the create majiq.conf file (not important at all...)', default='majiq.conf', type=str)
    parser.add_argument('--control', help='Name of the group which will be used as control', required=True)
    parser.add_argument('--debug', '-d', help='', action='store_true')
    args = parser.parse_args()

    inp, groups = (parse_input(args.meta))

    if args.control not in groups:
        raise IOError('The control group chosen "{0}" is not present in the groups specified: {1}'.format(args.control, groups))
    
    main(args, inp, groups, args.control, args.debug)
