#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information of percentage of mapping contigs/basepairs from the filtered assemblies
and produce a descriptive boxplot for mapped and unmapped contigs per assembler.

For each assembly, this script will output to the command line for each reference genome:
  * Assembler - assembler name (from fasta file name)
  * % mapped contigs - number of contigs that map to the reference files (and % over contigs > 1000bp)
  * % mapped bp - number of basepairs that map to the reference files (and % over total basepairs in
contigs > 1000bp)

Additionally, a Boxplot with boxes with size distribution of mapped contigs and scatter plot of the lenght o

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the unfiltered (raw) assembly files (ending in *.fasta)
  * Path to the mapped contings to the triple reference genomes (ending in *.paf)

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import sys
from plotly.offline import plot
import glob
import fnmatch
import plotly.graph_objects as go

#import commonly used functions from utils.py
import utils

def save_unmapped_contigs(df, assembly_files):
    """

    :param df: dataframed
    :param assembly_files: list of fasta files
    :return:
    """
    for assembler in sorted(df['Assembler'].unique()):

        fasta = utils.fasta_iter(fnmatch.filter(assembly_files, '*_' + assembler + '.*')[0])
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
    #add sanity check
    if len(assemblies) != len(mappings):
        print("Number of input files don't match.")
        sys.exit(0)



    print(','.join(['Assembler', '% mapped contigs', '% mapped bp']))

    # Dataframe with assembly info
    df = utils.parse_assemblies(assemblies, mappings)

    save_unmapped_contigs(df, assemblies)

    # Create plot
    fig = go.Figure()

    for assembler in sorted(df['Assembler'].unique()):

        contigs = df['Contig Len'][df['Assembler'] == assembler]
        mapped_contigs = df['Contig Len'][(df['Mapped'] == 'Mapped') & (df['Assembler'] == assembler)]

        print(','.join([assembler, f'{len(mapped_contigs)} ({(len(mapped_contigs)/len(contigs))*100:.2f}%)',
                        f'{sum(mapped_contigs)} ({(sum(mapped_contigs)/sum(contigs))*100:.2f}%)']))

        fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] == 'Mapped') & (df['Assembler'] == assembler)],
                            name=assembler, boxpoints='outliers',
                             boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))
        fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] == 'Unmapped') & (df['Assembler'] == assembler)],
                             name=assembler, boxpoints='all', pointpos=0, marker=dict(color='rgba(178,37,34,0.7)'),
                             line=dict(color='rgba(0,0,0,0)'), fillcolor='rgba(0,0,0,0)'))

    fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Contig size (Log bp)",
                      title="Contig size distribution per assembler (contigs over 1000 bp)",
                      plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))
    plot(fig)


if __name__ == '__main__':
    main()
