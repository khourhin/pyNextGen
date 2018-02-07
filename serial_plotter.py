#! /usr/bin/env python 

from bokeh.io import output_file, show
from bokeh.plotting import figure
import sys
import pandas as pd
import argparse
import numpy as np

# WARNINGS !!!
# NOW working with only one column (for hist)
# Assumes NO HEADER !


def summary(df, clim):
    "RUDIMENTARY SUMMARY FOR NOW"

    print(df.head())
    d = df['a'].value_counts().sort_values(ascending=False)
    print(d.head(min(d.shape[0], clim)))


def bar_hist(df):

    d = df['a'].value_counts().sort_values(ascending=False)

    f = figure(title="BLA",tools="save",
               background_fill_color="#E8DDCB")
    hist, edges = np.histogram(df['a'], bins=50)
    f.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
    show(f)

    print(d)

    
def get_input(sep):
    """
    Return a panda dataframe from the data passed by the stdin
    """
    
    df = pd.read_csv(sys.stdin, sep=sep, header=None)

    # Rename columns (necessary visibly to have column names for bokeh)
    df.columns = list('abcdefghijklmnopqrstuvwxyz'[:len(df.columns)])

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', '--sep',
                        help='The character used for field separator in the stdin table',
                        default='\t')
    parser.add_argument('--barhist',
                        help='Plot a bar plot histogram',
                        dest='action', action='store_const', const=bar_hist)
    parser.add_argument('-c', '--categories_lim',
                        help='Number of categories to print in the summary',
                        default=10, type=int)

    args = parser.parse_args()

    df = get_input(args.sep)
    args.action(df)
    summary(df, args.categories_lim)

# EXAMPLE USAGE
# tail -n +14 Arabidopsis_thaliana.TAIR10.34.gff3 |
# cut -f3,3 |
# python ~/Programming/pyNextGen/serial_plotter.py --barhist
