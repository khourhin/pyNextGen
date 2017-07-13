import logging
import click
import json

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.FileHandler('config_generator.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Job():
    def __init__(self, **kwargs):
        self.params = kwargs
        
    def __repr__(self):
        return "Job object: {0}".format(self.params)

    def generate_config(self, json_out):

        with open(json_out, 'w') as out:
            out.write(json.dumps(self.params, indent=4))


def init_job(**kwargs):
    """A simple parser to get NGS library information
    """
    
    job = Job(**kwargs)
    logger.info("Initializing: {0}".format(job))
    return job


@click.command()
@click.option('--bamdir', default=None, type=click.Path(exists=True), required=True, help='Path to the bam directory.')
@click.option('--fqdir', default=None, type=click.Path(exists=True), required=True, help='Path to the fastqs directory.')
@click.option('--strand', default=0, type=int, required=True, help='Strandedness: 0:not stranded; 1:stranded; 2: reverse stranded.')
@click.option('--paired', default=False, type=bool, required=True, help='Paired end: either True or False.')
@click.option('--json_out', default='snakemake.json', type=str, help='Path to the json output.')
@click.option('--gtf', default=None, type=click.Path(exists=True), required=True, help='Path to the gtf annotation file.')
def main(**kwargs):
    """A simple parser to get NGS library information
    """
    
    job = init_job(**kwargs)
    job.generate_config(kwargs['json_out'])

    
if __name__ == '__main__':
    main()
