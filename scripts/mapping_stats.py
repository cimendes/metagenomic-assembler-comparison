#!/usr/bin/env python3
"""
Adapted from https://raw.githubusercontent.com/rrwick/Long-read-assembler-comparison/master/scripts/assembly_stats.py

This script takes the following arguments (in this order):
    * Path with assemblies
    * Path with paf files

"""

import re
import sys
from itertools import groupby
import glob
import os
import fnmatch

REFERENCE_SEQUENCES = '../data/references/Zymos_Genomes_triple_chromosomes.fasta'


def fasta_iter(fasta_name):
    """
    Parse A Fasta File In Python
    :param fasta_name: fasta file name to be parsed
    :return: yield tuples of header, sequence
    """
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        try:
            seq = "".join(s.strip() for s in faiter.__next__())
        except StopIteration:
            print(headerStr)

        yield (headerStr, seq)


def get_NA50 (alignment_lengths):
    """

    :param alignment_lengths:
    :return:
    """
    sorted_lengths = sorted(alignment_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    na50 = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            na50 = contig_length
            break
    return na50


def get_covered_bases(covered_bases_list, ref_len):
    """

    :param covered_bases_list:
    :param ref_len:
    :return:
    """
    sorted_list = sorted(covered_bases_list, key=lambda x: x[0])

    covered_bases = set()

    for item in sorted_list:
        start, stop = map(int, item[:])

        # due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth
        # [0; ref_len][ref_len+1; 2*ref_len][(2*ref_len)+1; 3*ref_len]
        for base in range(start, stop):
            if base <= ref_len:
                covered_bases.add(base)
            elif base <= 2*ref_len:
                covered_bases.add(base-ref_len)
            else:
                covered_bases.add(base-(2*ref_len))
    return len(covered_bases) / ref_len


def get_lowest_window_identity(cigar, window_size):
    """

    :param cigar:
    :param window_size:
    :return:
    """
    lowest_window_id = float('inf')
    expanded_cigar = get_expanded_cigar(cigar)
    for i in range(0, len(expanded_cigar) - window_size):
        cigar_window = expanded_cigar[i:i+window_size]
        window_id = cigar_window.count('=') / window_size
        if window_id < lowest_window_id:
            lowest_window_id = window_id
    if lowest_window_id == float('inf'):
        return 0.0
    return lowest_window_id


def get_expanded_cigar(cigar):
    """

    :param cigar:
    :return:
    """
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def get_alignment_stats(paf_filename, ref_name, ref_length):
    """

    :param paf_filename:
    :param ref_name:
    :param ref_length:
    :return:
    """

    # Tracks the longest single alignment, in terms of the reference bases.
    longest_alignment = 0
    longest_alignment_id = 0
    longest_alignment_cigar = ''

    covered_bases = []  # dealing with overlapping regions :(

    # to calculate na50
    alignment_lengths = []

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                start, end = int(parts[7]), int(parts[8])
                alignment_lengths.append(end-start)

                matching_bases, total_bases = int(parts[9]), int(parts[10])
                cigar = [x for x in parts if x.startswith('cg:Z:')][0][5:]
                if end - start > longest_alignment:
                    longest_alignment = end - start
                    longest_alignment_id = matching_bases / total_bases
                    longest_alignment_cigar = cigar
                longest_alignment = max(longest_alignment, end - start)
                covered_bases.append([start, end])

    NA50 = get_NA50(alignment_lengths)

    relative_longest_alignment = longest_alignment / ref_length
    coverage = get_covered_bases(covered_bases, ref_length)
    lowest_window_id = get_lowest_window_identity(longest_alignment_cigar, 1000)

    return relative_longest_alignment, longest_alignment_id, lowest_window_id, coverage, NA50, len(alignment_lengths), \
           sum(alignment_lengths)


def main():
    """

    :return:
    """
    try:
        assemblies = sorted(glob.glob(sys.argv[1] + '/*'))
        mappings = glob.glob(sys.argv[2] + '/*')

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)

    for assembly_file in assemblies:
        filename = os.path.basename(assembly_file).split('.')[0].rsplit('_')[-1]

        print('\n\n------'+filename+'------\n')

        paf_file = fnmatch.filter(mappings, '*_' + filename + '.*')[0]

        # print specific file header to stdout.
        print(','.join(["reference", "reference length", "contiguity", "identity", "lowest identity",
                        "breadth of coverage", "aligned contigs", "aligned bp"]))

        reference = (x[1] for x in groupby(open(REFERENCE_SEQUENCES, "r"), lambda line: line[0] == ">"))

        for header in reference:
            header_str = header.__next__()[1:].strip().split()[0]
            seq = "".join(s.strip() for s in reference.__next__())

            contiguity, identity, lowest_window_identity, coverage, na50, aligned_seq, aligned_bp = \
                get_alignment_stats(paf_file, header_str, len(seq) / 3)  # divide by 3 for triple reference

            print(','.join([header_str, f'{len(seq)/3}', f'{contiguity:.2f}', f'{identity:.2f}',
                            f'{lowest_window_identity:.2f}', f'{coverage:.2f}', f'{aligned_seq}', f'{aligned_bp}']))


if __name__ == '__main__':
    main()