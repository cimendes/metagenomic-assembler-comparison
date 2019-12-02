#!/usr/bin/env python3
"""

This script takes csv tables by species and plots them as series of scatter plots

The headers must be
"""

import pandas as pd
import sys
import os
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import plotly.express as px


def get_file_names(files):
    result = []
    for file in files:
        result.append(os.path.basename(file).split('.')[0].replace('_', ' '))
    return result


def main():

    try:
        csv_tables = sys.argv[1:]

    # If no arguments were given, just print the header line.
    except IndexError:
        print("No files provided. Exiting.. ")
        sys.exit(0)

    #fig = make_subplots(rows=2, cols=len(csv_tables), shared_yaxes=True, subplot_titles=get_file_names(csv_tables))

    for file in csv_tables:
        filename = os.path.basename(file).split('.')[0].replace('_', ' ')
        print("Processing... " + filename)

        data = pd.read_csv(file)

        position = 1 if csv_tables.index(file)/4 < 1 else 2
        fig_to_add = px.scatter(data, x="Contigs", y="Breadth of coverage", color="Assembler", title=filename)
        fig_to_add.show()

        #trace_to_add = fig_to_add['data'][0]
        #fig.add_trace(trace_to_add, row=1, col=csv_tables.index(file) + 1)

    #fig.show()


if __name__ == '__main__':
    main()
