import logging
import click
import json

# Usage
# cd /home/ekornobis/Programming/pyNextGen/snake_make_recipes/test
# python ../config_generator.py -s stranded -p single -b ../../demo_data/ALL/bams -f ../../demo_data/ALL/fqs -g ~/data/genomes/ensembl/h_sapiens/hg37/Homo_sapiens.GRCh37.75.gtf
# snakemake --cores 20 -s ../post_bam.sm --configfile snakemake.json
# Limitations
# Check the CMD_OPTIONS dictionnary for adding other programs


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.FileHandler('/home/ekornobis/logs/config_generator.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

CMD_OPTIONS = {
    'qorts': {'single': '--singleEnded',
              'paired': '',
              'unstranded': '',
              'stranded': '--stranded --stranded_fr_secondstrand',
              'reverse': '--stranded',
    },
    'qualimap': {'single': '',
                 'paired': '-pe',
                 'unstranded': '-p non-strand-specific',
                 'stranded': '-p strand-specific-forward',
                 'reverse': '-p strand-specific-reverse',
    },
    'featureCount': {'single': '',
                     'paired': '-p',
                     'unstranded': '-s 0',
                     'stranded': '-s 1',
                     'reverse': '-s 2',
    },
}


class Job():
    """Object representing a job submission to be used in a snakemake
    pipeline run. Basically this object gather all parameters passed
    to the command line and generate a json file to be used as a
    --configfile in snakemake
    """
    
    def __init__(self, **kwargs):
        self.params = kwargs

        # Building argument string  
        for prog in CMD_OPTIONS:
            self.params[prog] = ' '.join([CMD_OPTIONS[prog][self.params[p]] for p in ['pair', 'strand']])

            
    def __repr__(self):
        return "Job object: {0}".format(self.params)

    
    def generate_config(self, json_out):

        with open(json_out, 'w') as out:
            out.write(json.dumps(self.params, indent=4))


@click.command()
@click.option('--bamdir','-b', default=None, type=click.Path(exists=True, resolve_path=True),
              help='Path to the bam directory.')
@click.option('--fqdir', '-f', default=None, type=click.Path(exists=True, resolve_path=True), required=True,
              help='Path to the fastqs directory.')
@click.option('--strand', '-s', type=click.Choice(['unstranded', 'stranded', 'reverse']), required=True,
              help='Strandedness: 0:not stranded; 1:stranded; 2: reverse stranded.')
@click.option('--pair', '-p', default=False, type=click.Choice(['single', 'paired']), required=True,
              help='Paired end: either True or False.')
@click.option('--json_out', '-o', default='snakemake.json', type=str,
              help='Path to the json output.')
@click.option('--gtf', '-g', default=None, type=click.Path(exists=True, resolve_path=True), required=True,
              help='Path to the gtf annotation file.')
@click.option('--star_index','-i', default=None, type=click.Path(exists=True, resolve_path=True),
              help='Path to the star index directory.')

def main(**kwargs):
    """A simple parser to get NGS library information
    """
    
    job = Job(**kwargs)
    logger.info("Initializing: {0}".format(job))
    job.generate_config(kwargs['json_out'])

    
if __name__ == '__main__':
    main()
