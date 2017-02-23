#### TO EDIT !!!!!!!!!!!

# table input example
# | CR_Inp_1.bam     | CR_Inp     |
# | CR_Inp_2.bam     | CR_Inp     |
# | CR_Inp_3.bam     | CR_Inp     |
# | CR_Inp_nas_1.bam | CR_Inp_nas |
# | CR_Inp_nas_2.bam | CR_Inp_nas |
# | CR_Inp_nas_3.bam | CR_Inp_nas |

READ_LEN = 50 # Set to maximum read length in your dataset
SAM_DIR = '/home/ekornobis/analysis/batsche/badeg/star_out_blattler/bams'
GENOME_NAME = 'hg19'
GENOME_PATH = '/home/ekornobis/data/genomes/ensembl/h_sapiens/hg37/Homo_sapiens.GRCh37.75.dna_sm.toplevel.fa'
GFF = '/home/ekornobis/data/genomes/ensembl/h_sapiens/hg37/majiq/ensembl_ENST_hg19.gff3'

NTHREADS=30
# NO SLASH AT THE END !
OUTDIR='majiq_blattler_1.0'
CONF_FILE='majiq.conf'

#STR_SPE='type=strand-specific'
# or
STR_SPE=''

# TO ADD
# MAKE BAM INDEX FIRST (ALWAYS FORGETTING ABOUT THOSE

import os
import argparse
import subprocess

def parse_input(input_tab):
    """Parse input for majiq. Input should have:
    $1: file paths
    $2: group name
    """
    inp = {}
    groups = set()
    with open(input_tab, 'r') as f:
        for line in f:
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
[experiments]\n""".format(READ_LEN, SAM_DIR, GENOME_NAME, GENOME_PATH, STR_SPE))

        for grp in groups:
            f.write('{0}={1}\n'.format(grp, ','.join([x
                                                      for x in inp
                                                      if inp[x][1] == grp])))

def main(inp, groups):
            
#    os.makedirs(OUTDIR)
    os.chdir(OUTDIR)

    print(inp)
    print(groups)
    
    create_conf(CONF_FILE, inp)
    
    cmd = "majiq build {0} -conf {1} --nthreads {2} --output .".format(GFF, CONF_FILE, NTHREADS)
    print(cmd)

#    subprocess.check_output(cmd, shell=True)

    for grp in groups:

        
        dot_majiq = ' '.join([x + '.majiq' for x in inp if inp[x][1] == grp])
        dot_splicegraph = ' '.join([x + '.splicegraph' for x in inp if inp[x][1] == grp])
        
        cmd =  'majiq psi {0} --nthreads {1} --output psi_{2} --name {2}'.format(dot_majiq, NTHREADS, grp)
        print(cmd)
 #       subprocess.check_output(cmd, shell=True)

        cmd = 'voila psi psi_{0}/{0}_psigroup.pickle -splice-graphs1 {1} -o voila_{0}'.format(grp, dot_splicegraph)
        print(cmd)
  #      subprocess.check_output(cmd, shell=True)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='')
    args = parser.parse_args()

    inp, groups = (parse_input(args.input))
    
    main(inp, groups)
