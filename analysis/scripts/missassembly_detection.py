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


def check_missassemblies(mappings):
    """
    :param df: Pandas Dataframe with stats for each contig
    :param mappings: list of paf files
    :return: dictionary
    """

    missmatch_dict = {}

    for paf_file in mappings:
        assembler = utils.get_assember_name(paf_file)

        with open(paf_file, 'r') as paf_fh:
            for line in paf_fh:
                line = line.split()

                contig, contig_len, query_start, query_end, strand = line[0:5]
                reference, reference_len, target_start, target_end = line[5:9]

                # a non-perfect alignment or a different number of residue matches
                if int(line[11]) != 0 or line[1] != line[9]:

                    cigar = line[-1]
                    exact_matches, snp, indel = utils.parse_cs(cigar)  # TODO - gap size to be adjusted by param

                    contig_dict = {'contig length': contig_len,
                                   'query start': int(query_start),
                                   'query end': int(query_end),
                                   'reference': reference,
                                   'reference length': int(reference_len)/3,
                                   'target start': utils.adjust_reference_coord(int(target_start),
                                                                                int(reference_len)/3),
                                   'target end': utils.adjust_reference_coord(int(target_end),
                                                                              int(reference_len)/3),
                                   'exact matches': exact_matches,
                                   'snp': snp,
                                   'indels': indel}

                    if assembler not in missmatch_dict.keys():
                        missmatch_dict[assembler] = {contig: [contig_dict]}
                    else:
                        if contig in missmatch_dict[assembler].keys():
                            missmatch_dict[assembler][contig].append(contig_dict)
                        else:
                            missmatch_dict[assembler][contig] = [contig_dict]

    return missmatch_dict


def evaluate_misassembled_contigs(dict):
    for assembler in dict.keys():
        for contig in dict[assembler].keys():
            if len(dict[assembler][contig]) > 1:  # why is contig 1 here?
                num_alignment_blocks = len(dict[assembler][contig])
                print(num_alignment_blocks)
                aligned_bases = 0
                contig_len = dict[assembler][contig][0]['contig length']
                reference = set()
                for alignment_block in dict[assembler][contig]:
                    print(alignment_block)
                    aligned_bases += (alignment_block['query end'] - alignment_block['query start'])




def main():
    try:
        mappings = [sys.argv[1]]  # TODO - to change

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)

    # Add correspondent reference to each dataframe contig
    misassembled_contigs = check_missassemblies(mappings)
    evaluate_misassembled_contigs(misassembled_contigs)


if __name__ == '__main__':
    main()