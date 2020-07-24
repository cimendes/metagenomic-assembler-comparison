#!/usr/bin/env python3
"""
TODO - Add description
     - Improve this....
"""

import sys
from itertools import groupby
import plotly.graph_objs as go
from plotly.offline import plot
import glob
import os
import fnmatch
import pandas as pd
import plotly.express as px

COLOURS = ['#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']
COLUMNS = ['Assembler','Contig', 'Contig Len', 'Mapped']


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip().split()[0]

        # join all sequence lines to one.
        try:
            seq = "".join(s.strip() for s in faiter.__next__())
        except StopIteration:
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

def get_mapped_contigs(paf_file):

    with open(paf_file) as f:
        mapped_contigs = [line.split()[0] for line in f]
    return mapped_contigs



def main():
    try:
        assemblies = glob.glob(sys.argv[1] + '/*')
        mappings = glob.glob(sys.argv[2] + '/*')

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)

    df = pd.DataFrame(columns=COLUMNS)

    data_to_plot_assembled = []
    legends = []
    data_to_plot_mapped = []
    data_to_plot_unmapped = []

    print(','.join(['Assembler', 'n Contigs', 'total bp', 'max contig size', 'n50', '% mapped contigs',
                    '% mapped bp']))

    for file in assemblies:

        filename = os.path.basename(file).split('.')[0].rsplit('_')[-1]
        mapped_contigs = get_mapped_contigs(fnmatch.filter(mappings, '*_'+filename+'.*')[0])

        contig_len = []
        total_bp = 0
        bp_in_mapped_contigs = 0
        mapped_contigs_to_plot = []
        unmapped_contigs_to_plot = []

        fasta = fasta_iter(file)
        for header, seq in fasta:
            contig_len.append(len(seq))
            total_bp += len(seq)
            if header in mapped_contigs:
                bp_in_mapped_contigs += len(seq)
                mapped_contigs_to_plot.append(len(seq))
                is_mapped = 'Mapped'
            else:
                unmapped_contigs_to_plot.append(len(seq))
                is_mapped = 'Unapped'

            df = df.append({'Assembler': filename, 'Contig': header, 'Contig Len': len(seq), 'Mapped': is_mapped},
                       ignore_index=True)

        data_to_plot_assembled.append(contig_len)
        legends.append(filename)
        data_to_plot_mapped.append(mapped_contigs_to_plot)
        data_to_plot_unmapped.append(unmapped_contigs_to_plot)

        print(','.join([filename, f'{len(contig_len)}', f'{total_bp}', f'{max(contig_len)}', f'{get_n50(contig_len)}',
                        f'{len(mapped_contigs)} ({(len(mapped_contigs)/len(contig_len))*100:.2f}%)',
                        f'{bp_in_mapped_contigs} ({(bp_in_mapped_contigs/total_bp)*100:.2f}%)']))

    df = df.reset_index()

    for assembler in df['Assembler'].unique():
        print(assembler)

    # Create plot
    # fig = go.Figure()
    #
    # for data_all, title, data_mapped, data_unmapped, colour in zip(data_to_plot_assembled, legends, data_to_plot_mapped, data_to_plot_unmapped, COLOURS):
    #     fig.add_trace(go.Box(y=data_all, name=title, boxpoints='outliers', boxmean=True, marker_color=colour))
    #     fig.add_trace(go.Box(y=data_mapped, name=title+'_mapped', boxpoints='outliers', boxmean=True, marker_color=colour))
    #     fig.add_trace(go.Box(y=data_unmapped, name=title+'_unmapped', boxpoints='outliers', boxmean=True, marker_color=colour))
    # fig.update_layout(showlegend=False, yaxis_type="log", yaxis_title="Contig size (bp)",
    #                   title="Contig size distribution per assembler (contigs over 1000 bp)")
    # plot(fig)

    fig = px.bar(df, x="Assembler", y="Contig Len", color="Mapped",
                 title="Contig size per assembler (contigs over 1000 bp)", yaxis_type="log")
    fig.show()




if __name__ == '__main__':
    main()