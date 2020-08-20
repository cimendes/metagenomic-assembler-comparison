#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assess assembly contamination/missassembly from a PAF file with the reads mapped to the assembly.
The reads must be identified with the origin genome in their name.

Purpose
-------
Get contamination/missassembly statistics for a directory  with the unmapped contigs of an assembly to the reference
genomes, and the mapping files (*.paf) of the original read data to the assembly files.

#TODO - Description of output

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the unmapped contigs of the assembly files (ending in *.fasta)
  * Path to the mapped reads to the complete assembly files (ending in *.paf)

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import sys
import glob
import fnmatch
from collections import Counter

#import commonly used functions from utils.py
import utils


def get_contig_info(paf_filename, fasta):
    """

    :param paf_filename:
    :param fasta:
    :return:
    """

    unmapped_contigs_iter = utils.fasta_iter(fasta)
    contigs = {c[0]: {'origin': "x", 'length': len(c[1])} for c in unmapped_contigs_iter}

    with open(paf_filename, "r") as paf_file:

        for line in paf_file:
            line_fields = line.split("\t")
            # Skip lines with no ID -  where do they come from?
            if not utils.is_number(line_fields[0].split('-')[0]):
                if line_fields[5] in contigs.keys():
                    if contigs[line_fields[5]]['origin'] == "x":
                        contigs[line_fields[5]]['origin'] = {line_fields[0].split('-')[0].split('_')[0].split('.')[0]}
                    else:
                        contigs[line_fields[5]]["origin"].add(line_fields[0].split('-')[0].split('_')[0].split('.')[0])

    return contigs


def get_plot_data(contigs_dict, filename):

    data = {}
    data[filename] = {"contaminated": [], "good": []}

    for contig in contigs_dict.keys():
        if len(contigs_dict[contig]["origin"]) > 1:  # more than one species mapped to a contig
            data[filename]["contaminated"].append(str(contigs_dict[contig]["length"]))
        else:
            data[filename]["good"].append(str(contigs_dict[contig]["length"]))

    print("Chimeric contigs: {} out of {} for {}".
          format(len(data[filename]["contaminated"]),
                 len(data[filename]["contaminated"]) + len(data[filename]["good"]),
                 filename))
    return data


def main():
    try:
        unmapped_contigs = glob.glob(sys.argv[1] + '/*.fasta')
        mappings = glob.glob(sys.argv[2] + '/*.paf')

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)
    # add sanity check
    if len(unmapped_contigs) != len(mappings):
        print("Number of input files don't match.")
        sys.exit(0)

    for assembly in unmapped_contigs:

        filename = utils.get_assember_name(assembly)
        paf_file = fnmatch.filter(mappings, '*_' + filename + '_*')[0]

        print("Parsing {}".format(filename))
        paf_contigs = get_contig_info(paf_file, assembly)

        cnt = Counter()

        for k, v in paf_contigs.items():
            cnt_ = Counter(paf_contigs[k]['origin'])
            cnt.update(cnt_)
        print(cnt)
        get_plot_data(paf_contigs, filename)


if __name__ == '__main__':
    main()
