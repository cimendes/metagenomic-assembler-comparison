#!/usr/bin/env python3
"""

This script assesses contamination in a metagenomic assemby. The reads must be identified with the
origin genome in their name ab«nd the alignment saved in PAF format

This script takes the following arguments (in this order):
  * list of alignment filenames PAF format (minimap2)

"""

import sys
import os
from collections import Counter
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import random


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_contig_info(paf_filename):

    contigs = {}

    with open(paf_filename, "r") as paf_file:

        for line in paf_file:
            line_fields = line.split("\t")
            # Skip lines with no ID -  where do they come from?
            if not is_number(line_fields[0].split('-')[0]):
                if line_fields[5] in contigs.keys():
                    contigs[line_fields[5]]["origin"].add(line_fields[0].split('-')[0].split('_')[0].split('.')[0])
                else:
                    contigs[line_fields[5]] = {"origin": {line_fields[0].split('-')[0].
                                                              split('_')[0].split('.')[0]}, "length": line_fields[6]}
    return contigs


def get_contaminated_contig_counts(filename, paf_contigs):

    data = {}

    # evaluate contamination
    data[filename] = {"contaminated": [], "good": []}
    to_report = []

    for contig in paf_contigs.keys():
        if len(paf_contigs[contig]["origin"]) > 1:  # more than one species mapped to a contig
            data[filename]["contaminated"].append(str(paf_contigs[contig]["length"]))
            to_report.append(str(paf_contigs[contig]["origin"]))
        else:
            data[filename]["good"].append(str(paf_contigs[contig]["length"]))

    print("Chimeric contigs: {} out of {} for {}".
          format(len(data[filename]["contaminated"]),
                 len(data[filename]["contaminated"]) + len(data[filename]["good"]),
                 filename))

    # Reporting....
    with open(filename + "_chimera.txt", "w") as outfile:
        counter = Counter(to_report)
        for key, value in counter.most_common():
            outfile.write(key + ':\t' + str(value) + '\n')

    return data


def main():
    try:
        paf_filename = sys.argv[1:]

    # If no arguments were given, just print the header line.
    except IndexError:
        print("No PAF mapping file provided. Exiting.. ")
        sys.exit(0)

    fig = make_subplots(rows=1, cols=len(paf_filename), shared_yaxes=True,
                        subplot_titles=("MEGAHIT", "metaSPADES", "SKESA", "SPADES"))

    for file in paf_filename:

        print("Parsing {}".format(file))
        paf_contigs = get_contig_info(file)
        filename = os.path.basename(file).split("_")[2]  # files are named reads_vs_{assembler}.paf

        data_dict = get_contaminated_contig_counts(filename, paf_contigs)

        fig.append_trace(go.Scatter(
                        x=[random.uniform(0, 1) for i in range(len(data_dict[filename]["contaminated"]))],
                        y=data_dict[filename]["contaminated"],
                        mode='markers',
                        name="{} chimera".format(filename),
                        opacity=0.7,
                    ), 1, paf_filename.index(file)+1)
        fig.append_trace(go.Scatter(
                        x=[random.uniform(0, 1) for i in range(len(data_dict[filename]["good"]))],
                        y=data_dict[filename]["good"],
                        mode='markers',
                        name=filename,
                        opacity=0.7
                    ), 1, paf_filename.index(file)+1)

    fig.update_xaxes(showticklabels=False)
    fig.show()


if __name__ == '__main__':
    main()
