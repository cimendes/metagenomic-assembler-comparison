#!/usr/bin/env python3
"""
Get the basic statistics for a directory of RAW assemblies.
Outputs the statistics for all the contigs in the assembly and
for the contigs with over 1000pb
"""

import sys
import glob

#import commonly used functions from utils.py
import utils


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

        filename = utils.get_assember_name(assembly_file)

        contigs, contigs_over_1000bp = get_contig_lists(utils.fasta_iter(assembly_file))

        n50_contigs = utils.get_N50(contigs)
        n50_contigs_over_1000bp = utils.get_N50(contigs_over_1000bp)

        print(','.join([filename, f'{len(contigs)}', f'{sum(contigs)}', f'{n50_contigs}', f'{max(contigs)}',
                        f'{len(contigs_over_1000bp)} ({(len(contigs_over_1000bp)/len(contigs))*100:.2f}%)',
                        f'{sum(contigs_over_1000bp)} ({(sum(contigs_over_1000bp)/sum(contigs))*100:.2f}%)',
                        f'{n50_contigs_over_1000bp}']))


if __name__ == '__main__':
    main()
