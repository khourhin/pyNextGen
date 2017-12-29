import pandas as pd
import os
import math

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

    def read_de_res(self, csv_file):

        df = pd.read_csv(csv_file, index_col=0)
        
        if list(df.columns) != ['external_gene_name', 'chromosome_name', 'start_position',
                                'end_position', 'strand', 'baseMean', 'log2FoldChange', 'lfcMLE',
                                'lfcSE', 'stat', 'pvalue', 'padj']:
            raise ValueError('The DR table does not have the correct headers.')
    
        return df

    def filter_fdr_fc(self, fdr_thres=0.05, fc_thres=0):
        """filter a DESeq2 result with FDR and logFC thresholds
        """

        fc_filt = self.df['fold_change'].abs() >= fc_thres
        fdr_filt = self.df['padj'] < fdr_thres

        return self.df[fc_filt.values & fdr_filt.values]

    def add_fold_change(self):
        """ Transform logFC in FC and add a column top df
        """

        self.df['fold_change'] = [math.copysign(2**x, x) for x in self.df['log2FoldChange']]


    def extract_gene_names_for_GO(self, df):
        """ Return a file with one Gene IDs per line for further GO analysis
        """
               
        outfile = os.path.splitext(self.name) + "geneIDs_list_for_GO.txt"
        with open(outfile, 'w') as f:
            for gene_ID in df.index:
                f.write(gene_ID + '\n')

    def describe(self):
        """ Describe the DR results
        """

        self.filter_fdr_logfc(self)
        

def get_all_results(deseq2_de_files, FDR_thres=0.05, logFC_thres=0):

    res_dict = {}
    
    for de_file in deseq2_de_files:
        df = read_de_res(de_file)
        df_filt = filter_FDR_logFC(df, FDR_thres=FDR_thres, logFC_thres=logFC_thres)
        res_dict[os.path.basename(df.index.name)] = describe_DEG(df_filt)
        extract_gene_names_for_GO(df_filt)

    return pd.DataFrame(res_dict).T
