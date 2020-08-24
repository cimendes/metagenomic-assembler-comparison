#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information on the distribution of gap sizes filtered assemblies and respective mapping of the contigs
to the reference sequences.
Produces a descriptive boxplot for the gap size per assembler.

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the filtered (min length of 1000bp) assembly files (ending in *.fasta)
  * Path to the mapped contigs to the triple reference genomes (ending in *.paf)

The triple bacterial reference files for the zymos mock community are available at
"../../data/references/Zymos_Genomes_triple_chromosomes.fasta"

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import sys
from plotly.offline import plot
import pandas as pd
import glob
import fnmatch
import plotly.graph_objects as go
from itertools import groupby

#import commonly used functions from utils.py
import utils

REFERENCE_SEQUENCES = os.path.join(os.path.dirname(__file__),
                                   '..', '..', 'data', 'references', 'Zymos_Genomes_triple_chromosomes.fasta')

COLUMNS = ['Assembler', 'Gap size']  # columns for dataframe


def get_gaps(paf_file, ref_name, ref_len):
    """
    Function to process the mapping (*.paf) file for a given reference and output a list with gap sizes from
    the alignment.
    :param paf_file: tabular file with alignment information for an assembler
    :param ref_name: reference name to filter from the paf_filename
    :param ref_len: int with expected reference length
    :return: gap_sizes: list with gap sizes in the assembly for the ref_name reference
    """
    covered_bases_list = []

    with open(paf_file) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                start, end = int(parts[7]), int(parts[8])
                covered_bases_list.append([start, end])

    covered_bases = set()
    for item in sorted(covered_bases_list, key=lambda x: x[0]):
        start, stop = map(int, item[:])
        for base in range(start, stop):
            # Due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth.
            # Therefore, the coordinates are adjusted as follows:
            # [0; ref_len][ref_len+1; 2*ref_len][(2*ref_len)+1; 3*ref_len]
            if base <= ref_len:
                covered_bases.add(int(base))
            elif base <= 2*ref_len:
                covered_bases.add(int(base-ref_len))
            else:
                covered_bases.add(int(base-(2*ref_len)))

    covered_bases = list(covered_bases)
    gaps = [[s, e]for s, e in zip(covered_bases, covered_bases[1:]) if s+1 < e]  # get list of gap sizes coords
    gap_sizes = [coord[1]-coord[0] for coord in gaps]
    return gap_sizes


def gap_size_distribution(assemblies, mappings):
    """
    Parses paf files and returns info on 'Assembler' and 'Gap size' as dataframe
    :param assemblies: list of assembly files
    :param mappings: list of paf files
    :return: pandas dataframe with gap sizes for each assembler
    """
    df = pd.DataFrame(columns=COLUMNS)

    for assembly_file in sorted(assemblies):

        filename = utils.get_assember_name(assembly_file)

        # iterator for reference files (sequence length is needed)
        references = (x[1] for x in groupby(open(REFERENCE_SEQUENCES, "r"), lambda line: line[0] == ">"))

        paf_file = fnmatch.filter(mappings, '*_' + filename + '.*')[0]
        for header in references:
            header_str = header.__next__()[1:].strip().split()[0]
            seq = "".join(s.strip() for s in references.__next__())

            gaps = get_gaps(paf_file, header_str, len(seq)/3)
            for gap in gaps:
                df = df.append({'Assembler': filename, 'Gap size': gap}, ignore_index=True)

    return df


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

    df = gap_size_distribution(assemblies, mappings)
    # Create plot - gap size distribution per assembler
    fig = go.Figure()

    for assembler in sorted(df['Assembler'].unique()):
        fig.add_trace(go.Box(x=df['Gap size'][df['Assembler'] == assembler],
                             name=assembler, boxpoints='outliers',
                             boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))

    fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Contig size (Log bp)",
                      title="Gap size distribution per assembler (contigs over 1000 bp)",
                      plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))
    plot(fig)


if __name__ == '__main__':
    main()
