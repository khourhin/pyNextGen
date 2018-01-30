import pandas as pd
import numpy as np
import os
import math
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# REQUIREMENTS
# The DR tables should be created with the new version of DR notebook:
# http://localhost:8888/notebooks/Programming/jupyter/good_practices/DEG_DESeq2_1.0.ipynb


class Deseq2Results(object):
    """Documentation for Deseq2Results

    """
    def __init__(self, csv_file):
        self.path = csv_file
        self.name = os.path.splitext(os.path.basename(csv_file))[0]
        self.df = self.read_de_res(csv_file)
        self.add_fold_change()

    def __repr__(self):
        return 'Deseq2_results: {}'.format(self.name)

    def read_de_res(self, csv_file):

        df = pd.read_csv(csv_file, index_col=0)
        
        if list(df.columns) != ['external_gene_name', 'chromosome_name', 'start_position',
                                'end_position', 'strand', 'baseMean', 'log2FoldChange', 'lfcMLE',
                                'lfcSE', 'stat', 'pvalue', 'padj']:
            raise ValueError('The DR table does not have the correct headers.')
    
        return df

    def add_fold_change(self):
        """ Transform logFC in FC and add a column top df
        """

        self.df['fold_change'] = [math.copysign(2**x, x) for x in self.df['log2FoldChange']]
    
    def filter_fdr_fc(self, fdr_thres=0.05, fc_thres=0):
        """filter a DESeq2 result with FDR and logFC thresholds
        """

        fc_filt = self.df['fold_change'].abs() >= fc_thres
        fdr_filt = self.df['padj'] < fdr_thres

        self.df = self.df[fc_filt.values & fdr_filt.values]

    def extract_gene_names_for_GO(self, df):
        """ Return a file with one Gene IDs per line for further GO analysis
        """
               
        outfile = os.path.splitext(self.name) + "geneIDs_list_for_GO.txt"
        with open(outfile, 'w') as f:
            for gene_ID in df.index:
                f.write(gene_ID + '\n')

    def summary(self):
        """ Describe the DR results
        STUB
        """

        print('Comparison: {}'.format(self.name))
        print('Up regulated: {}'.format((self.df.fold_change > 0).sum()))
        print('Down regulated: {}'.format((self.df.fold_change < 0).sum()))
        print('\n')
              
    def compare(self, d2res, verbose=False):
        """Compare two Deseq2Results objects.  For now this will plot venn
        diagrams for the comparisons considering UP and DOWN regulated
        genes separately

        If verbose, will return a dictionnary of dataframes with
        subset of genes for each comparisons
        """

        # Compare up regulated genes lists
        up1 = set(self.df.loc[self.df.fold_change > 0 ].index)
        up2 = set(d2res.df.loc[d2res.df.fold_change > 0 ].index)

        plt.figure()
        plot_venn(up1, up2, self.name, d2res.name, 'UP regulated')

        # Compare up regulated genes lists
        down1 = set(self.df.loc[self.df.fold_change < 0 ].index)
        down2 = set(d2res.df.loc[d2res.df.fold_change < 0 ].index)

        plt.figure()
        plot_venn(down1, down2, self.name, d2res.name, 'DOWN regulated')

        if verbose:
            return {self.name + '_specific_up': self.df.loc[up1 - up2],
                    d2res.name + '_specific_up': d2res.df.loc[up2 - up1],
                    self.name + '_specific_down': self.df.loc[down1 - down2],
                    d2res.name + '_specific_down': d2res.df.loc[down2 - down1],
                    self.name + '_' + d2res.name + '_common_up': self.df.loc[up1.intersection(up2)],
                    self.name + '_' + d2res.name + '_common_down': self.df.loc[down1.intersection(down2)]
            }

        
def plot_venn(set1, set2, label1, label2, title):
    v = venn2([set1, set2], set_labels=[label1, label2])
    plt.title(title)

# BROKEN         
def get_all_results(deseq2_de_files, FDR_thres=0.05, logFC_thres=0):

    
    res_dict = {}
    
    for de_file in deseq2_de_files:
        
        df = read_de_res(de_file)
        df_filt = filter_FDR_logFC(df, FDR_thres=FDR_thres, logFC_thres=logFC_thres)
        res_dict[os.path.basename(df.index.name)] = describe_DEG(df_filt)
        extract_gene_names_for_GO(df_filt)

    return pd.DataFrame(res_dict).T
