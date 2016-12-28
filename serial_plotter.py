from bokeh.charts import Histogram, output_file, show
import sys
import pandas as pd
import argparse


def make_hist(data):
    hist1 = Histogram(data)

    output_file('/tmp/serial_plotter_plot.html', title='Serial Plotter plot')
    show(hist1)


def get_input(sep):
    """
    Return a panda dataframe from the data passed by the stdin
    """
    df = pd.read_csv(sys.stdin, sep=sep)
    print(df)
    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s', '--sep',
                        help='The character used for field separator in the stdin table',
                        default='\t')
    parser.add_argument('-p', '--plot',
                        help='The type of plot to use',
                        default='hist')
    args = parser.parse_args()

    make_hist(get_input(args.sep))

