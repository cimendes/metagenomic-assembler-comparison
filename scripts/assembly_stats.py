#!/usr/bin/env python3
"""
Adapted from https://raw.githubusercontent.com/rrwick/Long-read-assembler-comparison/master/scripts/assembly_stats.py

This script takes the following arguments (in this order):
  * metagenomic assembly filename
  * metagenomic reads filename
  * alignment filenames to each sample in the metagenomic data
  * triple reference genome file

It outputs lots of information about the read set and the assembly. Run it
with no arguments to get the header line.
"""

import re
import sys
from itertools import groupby
from statistics import mean

def main():
    try:
        assembly_filename = sys.argv[1]
        sample_name = sys.argv[2]
        paf_filename = sys.argv[4]
        ref_sequence = sys.argv[5]
        assembler = sys.argv[6]

    # If no arguments were given, just print the header line.
    except IndexError:
        print('\t'.join(["name", "assembler", "contiguity", "identity", "lowest identity", "coverage"]))
        print('\t'.join(
            ["name", "assembler", "mean contiguity", " mean identity", " lowest identity", "mean coverage", "n contigs", "size",
             "n50"]), file=sys.stderr)
        sys.exit(0)

    fh_reference = open(ref_sequence, "r")

    entry = (x[1] for x in groupby(fh_reference, lambda line: line[0] == ">"))

    # get stats per reference
    total_length = 0
    contiguity_all = []
    identity_all = []
    lowest_identiy_all = []
    coverage_all = []

    for header in entry:
        header_str = header.__next__()[1:].strip().split()[0]
        seq = "".join(s.strip() for s in entry.__next__())
        total_length += len(seq)/3
        contiguity, identity, lowest_window_identity, coverage = get_alignment_stats(paf_filename, header_str, len(seq)/3)
        print('\t'.join([header_str, assembler, str(round(contiguity, 7)), str(round(identity, 7)), str(round(lowest_window_identity, 7)), str(round(coverage, 7))]))
        contiguity_all.append(contiguity)
        identity_all.append(identity)
        lowest_identiy_all.append(lowest_window_identity)
        coverage_all.append(coverage)

    contigs, size, n50 = get_assembly_stats(assembly_filename, total_length)

    result = [sample_name, assembler, f'{mean(contiguity_all):.7f}', f'{mean(identity_all):.7f}',
              f'{min(lowest_identiy_all):.7f}', f'{mean(coverage_all):.7f}', f'{contigs}', f'{size:.7f}', f'{n50:.7f}']

    print('\t'.join(result), file=sys.stderr)


def get_alignment_stats(paf_filename, ref_name, ref_length):

    # Tracks the longest single alignment, in terms of the reference bases. This can
    # exceed 100% in the case of overlapping ends.
    longest_alignment = 0
    longest_alignment_id = 0
    longest_alignment_cigar = ''

    # It also tracks the fraction of reference bases hit at least once (ranges from 0-100%).
    covered_bases = set()

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                start, end = int(parts[7]), int(parts[8])
                matching_bases, total_bases = int(parts[9]), int(parts[10])
                cigar = [x for x in parts if x.startswith('cg:Z:')][0][5:]
                if end - start > longest_alignment:
                    longest_alignment = end - start
                    longest_alignment_id = matching_bases / total_bases
                    longest_alignment_cigar = cigar
                longest_alignment = max(longest_alignment, end - start)
                for i in range(start, end):
                    covered_bases.add(i % ref_length)

    relative_longest_alignment = longest_alignment / ref_length
    coverage = len(covered_bases) / ref_length
    lowest_window_id = get_lowest_window_identity(longest_alignment_cigar, 1000)

    return relative_longest_alignment, longest_alignment_id, lowest_window_id, coverage


def get_lowest_window_identity(cigar, window_size):
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
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def get_assembly_stats(assembly_filename, ref_length):
    contig_lengths = sorted(get_contig_lengths(assembly_filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5 # why?
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return len(contig_lengths), total_length / ref_length, n50 / ref_length


def get_contig_lengths(filename):
    lengths = []
    with open(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


if __name__ == '__main__':
    main()
