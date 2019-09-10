#!/usr/bin/env python3
"""
Adapted from https://raw.githubusercontent.com/rrwick/Long-read-assembler-comparison/master/scripts/assembly_stats.py

This script takes the following arguments (in this order):
  * metagenomic assembly filename
  * sample name
  * alignment filename to the triple reference in PAF format (minimap2)
  * triple reference genome file
  * assembler name

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
        paf_filename = sys.argv[3]
        ref_sequence = sys.argv[4]
        assembler = sys.argv[5]

    # If no arguments were given, just print the header line.
    except IndexError:
        print('\t'.join(
            ["name", "assembler", "mean contiguity", " mean identity", " lowest identity", "mean coverage", "n contigs",
             "mean NA50", "mean aligned contigs", "mean aligned bp", "size", "n50"]), file=sys.stderr)
        sys.exit(0)

    fh_reference = open(ref_sequence, "r")

    entry = (x[1] for x in groupby(fh_reference, lambda line: line[0] == ">"))

    # get stats per reference
    contiguity_all = []
    identity_all = []
    lowest_identiy_all = []
    coverage_all = []
    na50_all = []
    aligned_contigs_all = []
    aligned_bp_all = []

    # print specific file header to stdout.
    print('\t'.join(["reference", "reference length", "contiguity", "identity", "lowest identity",
                     "breadth of coverage", "NA50", "aligned contigs", "aligned bp"]))

    for header in entry:

        header_str = header.__next__()[1:].strip().split()[0]
        seq = "".join(s.strip() for s in entry.__next__())

        contiguity, identity, lowest_window_identity, coverage, na50, aligned_seq, aligned_bp = get_alignment_stats(
            paf_filename, header_str, len(seq)/3)

        print('\t'.join([header_str, f'{len(seq)/3}', f'{contiguity:.4f}', f'{identity:.4f}', f'{lowest_window_identity:.4f}',
                         f'{coverage:.4f}', f'{na50}', f'{aligned_seq}', f'{aligned_bp}']))

        contiguity_all.append(contiguity)
        identity_all.append(identity)
        lowest_identiy_all.append(lowest_window_identity)
        coverage_all.append(coverage)
        na50_all.append(na50)
        aligned_contigs_all.append(aligned_seq)
        aligned_bp_all.append(aligned_bp)

    contigs, size, n50, min_contig_len, max_contig_len = get_assembly_stats(assembly_filename)

    print('\t'.join(
        ["assembler", "mean contiguity", " mean identity", "mean breadth of coverage",
         "mean NA50", "mean aligned contigs", "mean aligned bp", "contigs", "max contig length", "min contig lenght",
         "size bp", "n50"]), file=sys.stderr)

    result = [assembler, f'{mean(contiguity_all):.4f}', f'{mean(identity_all):.4f}',
              f'{mean(coverage_all):.4f}', f'{mean(na50_all):.4f}', f'{mean(aligned_contigs_all):.4f}',
              f'{mean(aligned_bp_all):.4f}' f'{contigs}', f'{max_contig_len}', f'{min_contig_len}', f'{size}', f'{n50}']

    print('\t'.join(result), file=sys.stderr)


def get_alignment_stats(paf_filename, ref_name, ref_length):

    # Tracks the longest single alignment, in terms of the reference bases.
    longest_alignment = 0
    longest_alignment_id = 0
    longest_alignment_cigar = ''

    covered_bases = 0

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
                covered_bases += end - start

    NA50 = get_NA50(alignment_lengths)

    relative_longest_alignment = longest_alignment / ref_length
    coverage = covered_bases / ref_length
    lowest_window_id = get_lowest_window_identity(longest_alignment_cigar, 1000)

    return relative_longest_alignment, longest_alignment_id, lowest_window_id, coverage, NA50, \
           len(alignment_lengths), sum(alignment_lengths)


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


def get_assembly_stats(assembly_filename):
    contig_lengths = sorted(get_contig_lengths(assembly_filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return len(contig_lengths), total_length, n50, max(contig_lengths), min(contig_lengths)


def get_contig_lengths(filename):

    lengths = []
    with open(filename, 'rt') as fasta_file:
        entry = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))

        for header in entry:
            header_str = header.__next__()[1:].strip().split()[0]
            seq = "".join(s.strip() for s in entry.__next__())
            lengths.append(len(seq))

    return lengths


def get_NA50 (alignment_lengths):
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


if __name__ == '__main__':
    main()
