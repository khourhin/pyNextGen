import pandas as pd
import os


def clean_sample_name(sample_path, extra_name_rm=[]):
    """ Clean sample names in feature count tables
    """

    new_name = os.path.basename(sample_path)
    new_name = os.path.splitext(new_name)[0]

    for pattern in extra_name_rm:
        new_name = new_name.replace(pattern, '')

    return new_name


def clean_feature_counts(feature_counts, extra_name_rm=[], drop_loc=True):
    """ Cleanup a featureCount count table
    - extra_name_rm: list of strings to remove from the header names
    - drop_loc: whether removing the 'Chr, Start, End, Strand, Length'
    columns
    """

    counts_df = pd.read_csv(feature_counts, sep='\t',
                            comment='#', index_col='Geneid')

    # Simplify the sample names
    counts_df.columns = [clean_sample_name(x, extra_name_rm=extra_name_rm)
                         for x in counts_df.columns]

    if drop_loc:
        counts_df.drop(['Chr', 'Start', 'End', 'Strand', 'Length'],
                       axis=1, inplace=True)
    
    counts_df.sort_index(inplace=True)

    return counts_df
    

    
