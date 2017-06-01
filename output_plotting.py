import argparse
import logging as log
import pandas as pd
import filehandling
import matplotlib.pyplot as plt

log.basicConfig(filename='example.log', level=log.INFO)
log.getLogger().addHandler(log.StreamHandler())

def plot_featureCount_output(prefix):
    """
    Visualization for the porcentage of successfull annotation/counts
    obtained from a featureCount run.
    prefix: the output prefix given in FeatureCount with the option -o 
    """

    feat_summary = prefix + '.summary'
    log.info('Checking: {}'.format(feat_summary))

    df = pd.read_csv(feat_summary, sep='\t', index_col=0)
    df.columns = [filehandling.simplify_path(x) for x in df.columns]

    df.plot()
    plt.show()

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--featureCount', help='Specify the prefix used for the FeatureCount output')
    args = parser.parse_args()

    plot_featureCount_output(args.featureCount)
