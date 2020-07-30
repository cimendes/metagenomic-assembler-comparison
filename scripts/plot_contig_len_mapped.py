#!/usr/bin/env python3
"""
Script to obtain information of percentage of mapping contigs/basepairs
and produce a descriptive boxplot.

Requires as input:
    - path for assembly data (filtered for over 1000bp ( files ending in *.fasta);
    - path for paf files (files ending in *.paf)
"""

import sys
from itertools import groupby
from plotly.offline import plot
import glob
import os
import fnmatch
import pandas as pd
import plotly.graph_objects as go

COLUMNS = ['Assembler','Contig', 'Contig Len', 'Mapped']  # columns for dataframe

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


def get_mapped_contigs(paf_file):
    """
    Gets list with the sizes of the mapped contigs.
    In the paf file, the first col is the contig name,
    the second is the contig size (excludes gaps)
    :param paf_file: path to the PAF file
    :return: list with contig sizes
    """
    with open(paf_file) as f:
        mapped_contigs = [line.split()[0] for line in f]
    return mapped_contigs


def parse_assemblies(assemblies, mappings):
    """
    Parses fastas and paf files and returns info on Assembler','Contig', 'Contig Len', 'Mapped' as dataframe
    :param assemblies: list of assembly files
    :param mappings: list of paf files
    :return: pandas dataframe
    """
    df = pd.DataFrame(columns=COLUMNS)

    for fasta_file in assemblies:

        filename = os.path.basename(fasta_file).split('.')[0].rsplit('_')[-1]
        mapped_contigs = get_mapped_contigs(fnmatch.filter(mappings, '*_' + filename + '.*')[0])

        fasta = fasta_iter(fasta_file)
        for header, seq in fasta:
            if header in mapped_contigs:
                is_mapped = 'Mapped'
            else:
                is_mapped = 'Unmapped'

            df = df.append({'Assembler': filename, 'Contig': header, 'Contig Len': len(seq), 'Mapped': is_mapped},
                           ignore_index=True)

    df = df.reset_index()

    return df


def save_unmapped_contigs(df, assembly_files):
    """

    :param df: dataframed
    :param assembly_files: list of fasta files
    :return:
    """
    for assembler in sorted(df['Assembler'].unique()):

        fasta = fasta_iter(fnmatch.filter(assembly_files, '*_' + assembler + '.*')[0])
        unmapped_contigs = list(df['Contig'][(df['Mapped'] == 'Unmapped') & (df['Assembler'] == assembler)])
        with open('unmapped_'+assembler+'.fasta', 'w') as fh:
            for header, seq in fasta:
                if header in unmapped_contigs:
                    fh.write(">" + header + "\n" + seq + "\n")


def main():
    try:
        assemblies = glob.glob(sys.argv[1] + '/*.fasta')
        mappings = glob.glob(sys.argv[2] + '/*.paf')

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)


    print(','.join(['Assembler', '% mapped contigs', '% mapped bp']))

    df = parse_assemblies(assemblies, mappings)

    save_unmapped_contigs(df, assemblies)

    # Create plot
    fig = go.Figure()

    for assembler in sorted(df['Assembler'].unique()):

        contigs = df['Contig Len'][df['Assembler'] == assembler]
        mapped_contigs = df['Contig Len'][(df['Mapped'] == 'Mapped') & (df['Assembler'] == assembler)]

        print(','.join([assembler, f'{len(mapped_contigs)} ({(len(mapped_contigs)/len(contigs))*100:.2f}%)',
                        f'{sum(mapped_contigs)} ({(sum(mapped_contigs)/sum(contigs))*100:.2f}%)']))

        fig.add_trace(go.Box(x=df['Contig Len'][df['Assembler'] == assembler], name=assembler, boxpoints='outliers',
                             boxmean=True, fillcolor='#D3D3D3', line=dict(color='#000000')))
        fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] == 'Unmapped') & (df['Assembler'] == assembler)],
                             name=assembler, boxpoints='all', pointpos=0, marker=dict(color='rgba(178,37,34,0.7)'),
                             line=dict(color='rgba(0,0,0,0)'), fillcolor='rgba(0,0,0,0)'))

    fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Contig size (Log bp)",
                      title="Contig size distribution per assembler (contigs over 1000 bp)",
                      plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))
    plot(fig)


if __name__ == '__main__':
    main()