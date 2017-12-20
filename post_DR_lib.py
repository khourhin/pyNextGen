import pandas as pd
import os
import math


def filter_FDR_logFC(df, FDR_thres=0.05, logFC_thres=0):
    """filter a DESeq2 result with FDR and logFC thresholds
    """

    logFC_filt =  df.filter(regex='.*\log2FoldChange').abs() >= logFC_thres
    FDR_filt = df.filter(regex='.*\padj') < FDR_thres

    # Verify there is only one column for each
    assert (logFC_filt.shape[1] == 1 & FDR_filt.shape[1] == 1)

    return df[logFC_filt.values & FDR_filt.values]

def describe_DEG(df):
    """Summarize a DESEQ2 DEG table
    """

    total_number_DE = int(df.shape[0])
    num_up_reg = int(((df.filter(regex='.*\log2FoldChange') > 0) == True).sum())
    num_down_reg = int(((df.filter(regex='.*\log2FoldChange') < 0) == True).sum())
    
    # Sanity check
    assert(total_number_DE == num_up_reg + num_down_reg)

    return({'total_num_DE':total_number_DE, 'num_down_reg':num_down_reg,
            'num_up_reg':num_up_reg})

def extract_gene_names_for_GO(df):
    """ Return a file with one Gene IDs per line for further GO analysis
    """

    assert(df.index.name)

    outfile = os.path.splitext(df.index.name)[0] + "geneIDs_list_for_GO.txt"
    with open(outfile, 'w') as f:
        for gene_ID in df.index:
            f.write(gene_ID + '\n')

def read_de_res(df_fname):

    df = pd.read_csv(df_fname, index_col=0)
    # Trick to keep the name of the file used to generate the df
    df.index.name = df_fname
    df.col

    return df

def add_fold_change(df):

    df = df['fold_change'] = [math.copysign(2**x, x) for x in df.filter(regex='.*\log2FoldChange')]
    

def get_all_results(deseq2_de_files, FDR_thres=0.05, logFC_thres=0):

    res_dict = {}
    
    for de_file in deseq2_de_files:
        df = read_de_res(de_file)
        df_filt = filter_FDR_logFC(df, FDR_thres=FDR_thres, logFC_thres=logFC_thres)
        res_dict[os.path.basename(df.index.name)] = describe_DEG(df_filt)
        extract_gene_names_for_GO(df_filt)

    return pd.DataFrame(res_dict).T
