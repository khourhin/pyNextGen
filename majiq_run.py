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

def main(args, inp, groups, debug=False):

    if not debug:
        os.makedirs(args.outdir)
        os.chdir(args.outdir)

    print(inp)
    print(groups)
    print(args.stranded)
    
    if not debug:
        create_conf(args.majiq_conf, inp)
    
    cmd = "majiq build {0} -conf {1} --nthreads {2} --output .".format(args.gff, args.majiq_conf, args.threads)
    print(cmd)

    if not debug:
        subprocess.check_output(cmd, shell=True)

    for grp in groups:

        
        dot_majiq = ' '.join([x + '.majiq' for x in inp if inp[x][1] == grp])
        dot_splicegraph = ' '.join([x + '.splicegraph' for x in inp if inp[x][1] == grp])
        
        cmd =  'majiq psi {0} --nthreads {1} --output psi_{2} --name {2}'.format(dot_majiq, args.threads, grp)
        print(cmd)
        
        if not debug:
            subprocess.check_output(cmd, shell=True)

        cmd = 'voila psi psi_{0}/{0}_psigroup.pickle -splice-graphs1 {1} -o voila_{0}'.format(grp, dot_splicegraph)
        print(cmd)

        if not debug:
            subprocess.check_output(cmd, shell=True)

    for i,j in itertools.combinations(groups, 2):
        print(i,j)
        dot_majiq_i = ' '.join([x + '.majiq' for x in inp if inp[x][1] == i])
        dot_majiq_j = ' '.join([x + '.majiq' for x in inp if inp[x][1] == j])

        cmd = 'majiq deltapsi -grp1 {majiq_i} -grp2 {majiq_j} --names {i} {j} --nthreads {t} --output dpsi_{i}_{j}'.format(majiq_i=dot_majiq_i, majiq_j=dot_majiq_j, i=i, j=j, t=args.threads)
        print(cmd)

        if not debug:
            subprocess.check_output(cmd, shell=True)

        dot_splicegraph_i = ' '.join([x + '.splicegraph' for x in inp if inp[x][1] == i])
        dot_splicegraph_j = ' '.join([x + '.splicegraph' for x in inp if inp[x][1] == j])
            
        cmd = 'voila deltapsi dpsi_{i}_{j}/{i}_{j}.deltapsi_quantify.pickle -splice-graphs1 {splice_i} -splice-graphs2 {splice_j} -o voila_{i}_{j}'.format(i=i, j=j, splice_i=dot_splicegraph_i, splice_j=dot_splicegraph_j)

        print(cmd)
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
    parser.add_argument('--debug', '-d', help='', action='store_true')
    args = parser.parse_args()

    inp, groups = (parse_input(args.meta))
    main(args, inp, groups, args.debug)
