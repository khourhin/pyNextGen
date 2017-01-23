from bokeh.charts import Bar, Histogram, output_file, show
import sys
import pandas as pd
import argparse

# WARNINGS !!!
# NOW working with only one column (for hist)
# Assumes NO HEADER !

def hist(data):

    hist_plot = Histogram(data)

    output_file('/tmp/serial_plotter_plot.html', title='Serial Plotter plot')
    show(hist_plot)

def bar_hist(data):

    d = data['a'].value_counts().sort_values(ascending=False)
    print(d.head(5))
    hist_plot = Bar(d, legend=None)
    output_file('/tmp/serial_plotter_plot.html', title='Serial Plotter plot')
    show(hist_plot)
    
def get_input(sep):
    """
    Return a panda dataframe from the data passed by the stdin
    """
    
    df = pd.read_csv(sys.stdin, sep=sep, header=None)

    # Rename columns (necessary visibly to have column names for bokeh)
    df.columns = list('abcdefghijklmnopqrstuvwxyz'[:len(df.columns)])

    
    print(df.head())
    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', '--sep',
                        help='The character used for field separator in the stdin table',
                        default='\t')
    parser.add_argument('--hist',
                        help='Plot a classic histogram',
                        dest='action', action='store_const', const=hist)
    parser.add_argument('--barhist',
                        help='Plot a bar plot histogram',
                        dest='action', action='store_const', const=bar_hist)

    args = parser.parse_args()

    args.action(get_input(args.sep))

