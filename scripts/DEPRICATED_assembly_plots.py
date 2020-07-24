#!/usr/bin/env python3
"""

This script assesses contamination in

This script takes the following arguments (in this order):
  * list of alignment filenames PAF format (minimap2)

"""

import sys
import os
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import random
from itertools import groupby
from math import log
import plotly.express as px
import pandas as pd
import plotly.figure_factory as ff


def parse_fasta(fasta_file):

    size = []
    fh_reference = open(fasta_file, "r")

    entry = (x[1] for x in groupby(fh_reference, lambda line: line[0] == ">"))

    for header in entry:
        header_str = header.__next__()[1:].strip().split()[0]
        seq = "".join(s.strip() for s in entry.__next__())
        size.append(log(len(seq)))

    return size


def main():

    try:
        fasta_files = sys.argv[1:]

    # If no arguments were given, just print the header line.
    except IndexError:
        print("No files provided. Exiting.. ")
        sys.exit(0)
    fig = make_subplots(rows=1, cols=len(fasta_files), shared_yaxes=True,
                        subplot_titles=("MEGAHIT", "metaSPADES", "SKESA", "SPADES"))

    for fasta in fasta_files:
        data = parse_fasta(fasta)
        fig.add_trace(go.Violin(y=data, box_visible=True, meanline_visible=True), 1, fasta_files.index(fasta) + 1)

    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(showlegend=False)
    fig.show()


if __name__ == '__main__':
    main()



