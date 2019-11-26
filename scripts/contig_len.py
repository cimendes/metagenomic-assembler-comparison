#!/usr/bin/env python3
"""
TODO - Add description
     - Improve this....
"""

import sys
from itertools import groupby
from math import log
import plotly.graph_objs as go
from plotly.offline import plot


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()  # drop the ">"

        # join all sequence lines to one.
        try:
            seq = "".join(s.strip() for s in faiter.__next__())
        except StopIteration:  # some multifasta break here for some reason
            print(headerStr)
        yield (headerStr, seq)


def get_n50(contig_lengths):
    contig_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return n50


def main():
    try:
        assembly_filename = sys.argv[1:]

    except IndexError as e:
        print(e, file=sys.stderr)
        sys.exit(0)

    data_to_plot = []
    legends = []

    print(','.join(['Assembler', 'n Contigs', 'total bp', 'max contig size', 'n50', '% contigs>1000bp',
                    '% bp in contigs>1000bp', 'n50 in contigs>1000bp']))

    for file in assembly_filename:

        contig_len_log = []  # unused
        contig_len = []
        total_bp = 0
        good_contigs = []
        bp_in_good_contigs = 0

        fasta = fasta_iter(file)
        for header, seq in fasta:
            if len(seq) > 1000:
                good_contigs.append(len(seq))
                bp_in_good_contigs += len(seq)
            contig_len.append(len(seq))
            contig_len_log.append(log(len(seq)))  # log scale
            total_bp += len(seq)

        data_to_plot.append(contig_len)
        legends.append(file.split('.')[0].rsplit('_')[1])

        print(','.join([file.split('.')[0].rsplit('_')[1], f'{len(contig_len):.4f}', f'{total_bp}',
                        f'{max(contig_len):.4f}', f'{get_n50(contig_len)}',
                        f'{len(good_contigs)/len(contig_len):.4f}',
                        f'{bp_in_good_contigs/total_bp:.4f}', f'{get_n50(good_contigs)}']))

    # Create plot
    fig = go.Figure()

    for data, title in zip(data_to_plot, legends):
        fig.add_trace(go.Box(y=data, name=title, boxpoints='outliers', boxmean=True, marker_color='rgb(8, 81, 156)'))
    fig.update_layout(showlegend=False, yaxis_type="log", yaxis_title="Contig size (bp)",
                      title="Contig size distribution per assembler.")
    plot(fig)


if __name__ == '__main__':
    main()
