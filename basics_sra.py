import argparse
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# USAGE EX:
#cd /home/ekornobis/analysis/batsche/methyleu/2.0_pileuping
#SRAS=$(ls pups/ | cut -f1,1 -d'_' | sort -u)
#grep CD44 ~/data/genomes/ensembl/h_sapiens/hg37/Homo_sapiens.GRCh37.75.gtf | head -n1 > j1
#for sra in $SRAS; do python ~/code/batsche/methyleu/get_gene_specific_sra.py $sra j1 -f 10000; done

class Gene():
    """
    Gene object parsed from a gtf entry
    """

    def __init__(self, gtf_entry):
        gtf_entry = gtf_entry.split('\t')
        self.chrom = gtf_entry[0]
        self.start = int(gtf_entry[3])
        self.end = int(gtf_entry[4])

        # Parsing gtf 
        attrs = [x.strip().split(' ') for x in gtf_entry[-1].split(';')[:-1]]
        attrs = {x[0]:x[1].strip('"') for x in attrs}
        
        self.ensemblID = attrs['gene_id']
        self.gene_name = attrs['gene_name']
        
    def __str__(self):

        string = "ID: {ensemblID}\tNAME: {gene_name}\t{chr}:{start}-{end}".format(ensemblID=self.ensemblID,
                                                                                 gene_name=self.gene_name,
                                                                                 chr=self.chrom,
                                                                                 start=self.start,
                                                                                 end=self.end)
        return string

def get_bam_for_gene(sra, gene, flank):
    """
    """
    log.info(gene)
    bam = "{sra}_{gene_name}.bam".format(sra=sra, gene_name=gene.gene_name)
    cmd = 'sam-dump --aligned-region {chr}:{start}-{end} {sra} | samtools view -bh > {bam}'.format(chr=gene.chrom,
                                                                                                   start=gene.start - flank,
                                                                                                   end=gene.end + flank,
                                                                                                   bam=bam,
                                                                                                   sra=sra)
    log.info(cmd)
    subprocess.check_output(cmd, shell=True)
    
def get_all_bams(sra, gene_tab, flank):
    """
    """
    
    with open(gene_tab, 'r') as f:
        for line in f:
            g = Gene(line)
            get_bam_for_gene(sra, g, flank)
        
def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument('sra', help='The SRR identifier to get the data from.')
    parser.add_argument('gene_tab', help='A gtf file extract with only gene rows for the selected genes')
    parser.add_argument('--flank', '-f', help='Number of bases to add on each flank of the gene interval', default=0, type=int)

    args = parser.parse_args()

    get_all_bams(args.sra, args.gene_tab, args.flank)

if __name__ == "__main__":
    main()
