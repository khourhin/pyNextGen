import urllib
import argparse
from bs4 import BeautifulSoup
import json
import logging

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def get_geo_id(gse):
    """Retrieve UID id from the Gene Expression Omnibus
    corresponding to the dataset id 'gse'
    """

    params = urllib.parse.urlencode({'db': 'gds',
                                     'term': gse + '[ACCN]',
                                     'usehistory': 'y',
                                     'retmode': 'json',
                                     'email': 'ekornobis@gmail.com',
                                     'tool': 'pythonScript'})
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{}'.format(params)
    log.info('URL Query: {}'.format(url))
    
    with urllib.request.urlopen(url) as response:
        res = json.loads(response.read().decode('utf-8'))
        
        log.debug(json.dumps(res, indent=4, sort_keys=True))
        webenv = res['esearchresult']['webenv']
        query_key = res['esearchresult']['querykey']
        
        return webenv, query_key


def get_sras_ids(webenv, query_key):
    """From a previous entrez query using get_geo_id, extract from gds
    identifiers the corresponding sras identifier and return them
    """
    
    params = urllib.parse.urlencode({'dbfrom': 'gds',
                                     'db': 'sra',
                                     'WebEnv': webenv,
                                     'query_key': query_key,
                                     'cmd': 'neighbor_history',
                                     'retmode': 'json',
                                     'email': 'ekornobis@gmail.com',
                                     'tool': 'pythonScript'})
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?{}'.format(params)

    log.info('URL Query: {}'.format(url))
    
    with urllib.request.urlopen(url) as response:
        
        res = json.loads(response.read().decode('utf-8'))
        log.debug(json.dumps(res, indent=4, sort_keys=True))
        webenv = res['linksets'][0]['webenv']
        query_key = res['linksets'][0]['linksetdbhistories'][0]['querykey']

        
    params = urllib.parse.urlencode({'db': 'sra',
                                     'WebEnv': webenv,
                                     'query_key': query_key,
                                     'usehistory': 'y',
                                     'retmode': 'json',
                                     'email': 'ekornobis@gmail.com',
                                     'tool': 'pythonScript'})
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{}'.format(params)

    with urllib.request.urlopen(url) as response:
        
        res = json.loads(response.read().decode('utf-8'))

        # for entry in res['result']:
        #     print(res['result'][entry]['runs'])
        res = res['result']
        log.debug(json.dumps(res, indent=4, sort_keys=True))


    # params = urllib.parse.urlencode({'db': 'sra',
    #                                 'save': 'efetch',
    #                                 'rettype': 'runinfo',
    #                                 'WebEnv': webenv,
    #                                 'query_key': query_key,})
    # url = 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi'.format(params)
        
    # with urllib.request.urlopen(url) as response:
        
    #     res = response.read()
    #     print(res)

    # urllib.request.urlretrieve(url, filename='test.csv')
        
        
        
def download_sras(sra_ids):
        """From SRA ids obtained with get_sras_ids, download the
        corresponding sras
        """
        
        params = urllib.parse.urlencode({'db':'sra', 'save':'efetch', 'rettype':'runinfo', 'WebEnv':webenv, 'query_key':query_key,})
    
        url = 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?{}'.format(params)
        
        

def parse_arguments():
        parser = argparse.ArgumentParser(description='')

        return parser.parse_args()
        

if __name__ == '__main__':

    web_env, query_key = get_geo_id('GSE66763')
    get_sras_ids(web_env, query_key)
