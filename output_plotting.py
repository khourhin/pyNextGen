import argparse
import logging as log
import pandas as pd
import filehandling
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

log.basicConfig(filename='example.log', level=log.DEBUG)
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

def compare_sets(set_file):

    set_list = []
    with open(set_file) as f:
        for line in f:
            set_list.append(set(line.split()))
        
    log.debug('Set List:{}'.format(set_list))
    
    if len(set_list) == 2:
        venn2(set_list)
    if len(set_list) == 3:
        venn3(set_list)
    plt.show()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--featureCount', help='Specify the prefix used for the FeatureCount output')
    parser.add_argument('--setCompare', help='Specify a file with the sets to compare, one on each line')
    args = parser.parse_args()

    if args.featureCount:
        plot_featureCount_output(args.featureCount)

    if args.setCompare:
        compare_sets(args.setCompare)
