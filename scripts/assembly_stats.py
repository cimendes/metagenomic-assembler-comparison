#!/usr/bin/env python3
"""
Get the basic statistics for a directory of RAW assemblies.
Outputs the statistics for all the contigs in the assembly and
for the contigs with over 1000pb
"""

import sys
import os
import glob
import fnmatch
from itertools import groupby


def get_assember_name(assembly_file):
    """
    get assembler name from filename. Expected format: `XX_<AssemblerName>.fasta`
    :param assembly_file: path
    :return: filename
    """
    return os.path.basename(assembly_file).split('.')[0].rsplit('_')[-1]


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


def get_N50(contig_lengths):
    """
    Calculate n50 from a list of contig lengths
    :param contig_lengths: list with contig lengths
    :return: int of the n50 of the contigs
    """
    target_length = sum(contig_lengths) * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return n50


def get_contig_lists(fasta):
    """
    Get basic assembly statistics from a fasta
    :param fasta: yield tuples of header, sequence
    :return:
    """

    contigs_len_over_1000 = []  # list of contigs with len > 1000
    contigs_len = []  # list with all contig lens

    for header, seq in fasta:
        if len(seq) > 1000:
            contigs_len_over_1000.append(len(seq))
        contigs_len.append(len(seq))

    return contigs_len, contigs_len_over_1000


def main():
    """
    in a directory with assemblies (ended in .fasta),
    calculate the assembly statistics (number of contigs, total number of basepairs, max contig size, n50
    for all contigs per assembly, including seperate stats for the contigs with over 1000bp
    """
    try:
        assemblies = sorted(glob.glob(sys.argv[1] + '/*.fasta'))
    except IndexError as e:
        print(e, "Directory not found.")
        sys.exit(0)

    print(','.join(['Assembler', 'Contigs', 'basepairs', 'Max contig size', 'n50', 'contigs>1000bp (%)',
                    ' bp in contigs>1000bp (%)','n50 in contigs>1000bp']))

    for assembly_file in assemblies:

        filename = get_assember_name(assembly_file)

        contigs, contigs_over_1000bp = get_contig_lists(fasta_iter(assembly_file))

        n50_contigs = get_N50(contigs)
        n50_contigs_over_1000bp = get_N50(contigs_over_1000bp)

        print(','.join([filename, f'{len(contigs)}', f'{sum(contigs)}', f'{n50_contigs}', f'{max(contigs)}',
                        f'{len(contigs_over_1000bp)} ({(len(contigs_over_1000bp)/len(contigs))*100:.2f}%)',
                        f'{sum(contigs_over_1000bp)} ({(sum(contigs_over_1000bp)/sum(contigs))*100:.2f}%)',
                        f'{n50_contigs_over_1000bp}']))


if __name__ == '__main__':
    main()
