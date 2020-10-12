#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
For each contig, this script will evaluate if it's a missassembly and classy it into these main types:
    * Insertion
    * Deletion
    * Missense
    * Inversion
    * Chimera
    * Translocation
    * ...

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the metagenomic assembly files (ending in *.fasta)
  * Path to the mapped contigs to the triple reference genomes (ending in *.paf)

The triple bacterial reference files for the zymos mock community are available at
"../../data/references/Zymos_Genomes_triple_chromosomes.fasta"

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import random
import sys
from itertools import groupby
import glob
import os
import re
import fnmatch
import math
import pandas as pd
from plotly import subplots
from plotly.offline import plot
import plotly.graph_objects as go

#import commonly used functions from utils.py
import utils

REFERENCE_SEQUENCES = os.path.join(os.path.dirname(__file__),
                                   '..', '..', 'data', 'references', 'Zymos_Genomes_triple_chromosomes.fasta')


def add_matching_ref(df, mappings):
    """
    For each contig in the df, adds the correspondent reference if the contig is mapped. Drops the unmapped contigs.
    :param df: Pandas Dataframe with stats for each contig
    :param mappings: list of paf files
    :return: Pandas Dataframe with stats for each contig with reference info instead of 'Mapped' and rows with
    unmapped contigs removed
    """
    #paf_dict = {'Assembler': {'Contig': '', 'Data': [{}]}}

    for assembler in sorted(df['Assembler'].unique()):
        paf_file = fnmatch.filter(mappings, '*_' + assembler + '.*')[0]
        #paf_dict['Assembler'] = assembler
        with open(paf_file, 'r') as paf_fh:
            for line in paf_fh:
                line = line.split()
                if int(line[11]) == 0:
                    if line[1] != line[9]:  # number of residue matches
                        cigar = line[-1]
                        snp, indel, insertions, deletions = utils.parse_cs(cigar, 50)
                        print(cigar)
                        print('SNP: {}'.format(snp))
                        print('InDel: {}'.format(indel))
                        print('Insertions: {}'.format(len(insertions)))
                        print('Deletions: {}'.format(len(deletions)))

def main():
    try:
        assemblies = [sys.argv[1]]  # changed!
        mappings = [sys.argv[2]]

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)

    try:
        if sys.argv[3] == "--print-csv":
            print_csv = True
    except IndexError as e:
        print_csv = False

    # Dataframe with assembly info
    df = utils.parse_assemblies(assemblies, mappings)

    # Add correspondent reference to each dataframe contig
    df = add_matching_ref(df, mappings)

if __name__ == '__main__':
    main()