# Used for GPHN ribosome profiling analysis
# See http://localhost:8888/notebooks/code/allemand/gphn/10.5_last_Ribo.ipynb for usage

### SPECIFIC ###

import os
import pysam
import math
import numpy as np
import pandas as pd

from bokeh.plotting import figure, show
from bokeh import palettes

from qgrid import show_grid


def ribo_pause(df_row):
    """Infer ribosome pausing sites modified from Ishimura 2014
    Only change is adding longer reads (>31 bases)"""
    pos = df_row.ali_mid
    readlen = df_row.ali_len
    posframe = pos % 3
    lenrem = readlen % 6
    
    if readlen == 26 or readlen == 32:
        if posframe == 2 and lenrem == 2: offset = -3
        elif posframe == 0 and lenrem == 2: offset=-1
        elif posframe == 1 and lenrem == 2: offset=-2
        else: return -1    
        
    elif readlen == 27 or readlen == 33:
        if posframe == 2 and lenrem == 3: offset=0
        elif posframe == 0 and lenrem == 3: offset=-1
        elif posframe == 1 and lenrem == 3: offset=-2
        else: return -2

    elif readlen == 28 or readlen == 34:
        if posframe == 0 and lenrem == 4: offset=-1
        elif posframe == 1 and lenrem == 4: offset=-2
        elif posframe == 2 and lenrem == 4: offset=-3
        else: return -3
        
    elif readlen == 29 or readlen == 35:
        if posframe == 0 and lenrem == 5: offset=-1
        elif posframe == 1 and lenrem == 5: offset=-2
        elif posframe == 2 and lenrem == 5: offset=0
        else: return -4
        
    elif readlen == 30 or readlen == 36:
        if posframe == 1 and lenrem == 0: offset=-2
        elif posframe == 2 and lenrem == 0: offset=-3
        elif posframe == 0 and lenrem == 0: offset=-1
        else: return -5
        
    elif readlen == 31 or readlen == 25:
        if posframe == 1 and lenrem == 1: offset=-2
        elif posframe == 2 and lenrem == 1: offset=0
        elif posframe == 0 and lenrem == 1: offset=-1
        else: return -6
    
    else:
        newpos = np.nan
        return -7
            
    newpos = pos + offset
    return newpos


def get_ali_stats(ali):
    """Extract alignment information for all alignements"""
    
    ali = {'ref_start': ali.reference_start, 'ref_end': ali.reference_end,
           'read_len': ali.qlen, 'ali_len': ali.reference_length,
           'ali_mid': math.ceil((ali.reference_start + ali.reference_end)/2)}
    return ali


def get_all_ali_stats(bam_list):
    """Extract alignment information for all bams"""
    ribo_pos = {}
    for bam in bam_list:
        
        for ali in pysam.AlignmentFile(bam):
            ribo_pos[(bam, ali.reference_name, ali.qname)] = get_ali_stats(ali)
    
    ribo_df = pd.DataFrame(ribo_pos).T
    ribo_df.index.set_names(['bam', 'iso', 'read'], inplace=True)
    return ribo_df


def ribo_profiling_plot(pos_df, ref_fasta, bams, iso, orf_finder_csv=None, method='pos'):
    """Interactive plot of the coverage after Ishimura protocol method:
    can be either 'pos' ie number of reads at each position or
    'read_percentage' ie percentage of reads over the transcript at the particular position
    """

    tools = ['crosshair', 'wheel_zoom', 'box_zoom', 'reset']
    
    fas = {seq.name: seq.sequence for seq in pysam.FastxFile(ref_fasta)}
    iso_len = len(fas[iso])
    
    p1 = figure(plot_width=1550, plot_height=600, tools=tools)
    
    for bam, color in zip(bams, palettes.Category10[10][:len(bams)]):
        df = pos_df.loc[bam, iso].pos.value_counts()
        df = df.reindex(np.arange(iso_len)).fillna(0).reset_index()
        df['read_percentage'] = df['pos'] / sum(df['pos']) * 100
        
        # Avoiding "source" because its messing up the legend !
        p1.line(df['index'], df[method], color=color, muted_color=color,
                muted_alpha=0.1, legend=bam)
    
    p1.legend.click_policy = "mute"
    
    # Sam positions are 1 based (not 0 based). +1 to correct sequence position
    p1.text(x=np.arange(iso_len) + 1, y=np.repeat(-0.1, iso_len),
            text=list(fas[iso]), text_align='center')

    if orf_finder_csv:
        orf_df = parse_orf_finder(orf_finder_csv)
        orf_df = orf_df.loc[iso]

        p1.text(x=orf_df.start, y=np.repeat(10, len(orf_df.start)),
                text=["*" + str(x) for x in orf_df.start])
    
    show(p1, notebook_handle=True)
    return(p1)

### SPECIFIC
def parse_orf_finder(csv):
    """Parse orf finder orf summary csv
    """

    df = pd.read_csv(csv, index_col=[0,1])
    df['nuc_length'] = [len(seq) for seq in df['CDS']]
    df['start'] = [int(x.split(':')[-2]) + 718 for x in df['info']]
    df['end'] = [int(x.split(':')[-1]) + 718 for x in df['info']]
    df = df[df['start'] < df['end']]
    return df
