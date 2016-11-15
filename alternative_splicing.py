import logging
import subprocess
import os

# For starter only work with single reads


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def get_gdir(gfas):
    """
    Return the path to the gdir folder for star
    """
    gdir = os.path.join(os.path.split(gfas)[0],
                        'star_index_' + os.path.basename(gfas))

    return gdir


def star_index(gfas, nthreads):
    """
    Produce the star index if not present already in a folder in the same
    folder as gfas.
    """

    gdir = get_gdir(gfas)
    os.makedirs(gdir)

    # This is a syntax to concatenate strings (can add '+' for verbose)
    cmd = ('STAR '
           '--runMode genomeGenerate '
           '--genomeDir {0} '
           '--genomeFastaFiles {1} '
           '--runThreadN {2}').format(gdir, gfas, nthreads)
    log.info('Executing cmd:' + cmd)

    p = subprocess.check_output(cmd, shell=True)
    print(p)

def star_align_1pass(gfas, gtf, read_list1, nthreads):
    """
    Produce the read alignment with STAR
    """

    gdir = get_gdir(gfas)

    if not os.path.exists(gdir):
        star_index(gfas, nthreads)
    else:
        log.info("Using index already present in {}".format(gdir))

    for rd in read_list1:
        outdir = 'star_out'
        if os.path.exists(outdir):
            raise IOError("The directory {} already exists. To not overwrite it, Star alignment is canceled.".format(outdir))

        os.makedirs(outdir)
        prefix = os.path.join(outdir, os.path.splitext(os.path.basename(rd))[0])
        cmd = ('STAR '
               '--genomeDir {0} '
               '--sjdbGTFfile {1} '               
               '--readFilesIn {2} '
               '--outFileNamePrefix {3} '
               '--outSAMtype BAM SortedByCoordinate '
               '--runThreadN {4}').format(gdir, gtf, rd, prefix, nthreads)
        log.info('Executing cmd:' + cmd)

        p = subprocess.check_output(cmd, shell=True)
        print(p)

if __name__ == '__main__':

    GFAS = '/home/ekornobis/data/genomes/c_elegans/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa'
    NTHREADS = 20
    READS=['/home/ekornobis/data/demo_data/illumina/c_elegans/SRR019721.fastq']
    GTF = '/home/ekornobis/data/genomes/c_elegans/Caenorhabditis_elegans.WBcel235.86.gtf'

    star_align_1pass(GFAS, GTF, READS, NTHREADS)
