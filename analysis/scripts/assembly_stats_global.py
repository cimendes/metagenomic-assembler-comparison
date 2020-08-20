#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get the basic statistics for a directory of RAW assemblies.
Outputs the statistics for all the contigs in the assembly and for the contigs with over 1000pb

Purpose
-------
Get the basic statistics for a directory of RAW assembly files (*.fasta), unfiltered for length.

For each assembly, this script will output to the command line for each reference genome:
  * Assembler - assembler name (from fasta file name)
  * Contigs - number of contigs in the assembly
  * basepairs - total number of nucleotides in the assembly
  * Max contig size - size of the largest contig in the assembly
  * n50 - sequence length of the shortest contig at 50% of the total assembly length
  * contigs>1000bp (%) - number of contigs with size >= 1000 bp (and % over "Contigs")
  * bp in contigs>1000bp (%) - total number of nucleotides in contigs with size >= 1000 bp (and % over "basepairs")
  * n50 in contigs>1000bp - sequence length of the shortest contig at 50% of the total length of contigs with ´
size >= 1000 bp

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the unfiltered (raw) assembly files (ending in *.fasta)

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import sys
import glob

#import commonly used functions from utils.py
import utils


def get_contig_lists(fasta):
    """
    From a fasta iterator, get lists with contig lengths
    :param fasta: yield tuples of header, sequence
    :return:
        - contig_len: list with all contig lenghts in the assembly
        - contigs_len_over_1000: list with contig lenghts filtered for a minimum of 1000 nucleotides
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
    in a directory with assemblies (ended in "*.fasta"),
    calculate the assembly statistics (number of contigs, total number of basepairs, max contig size, n50)
    for all contigs per assembly, including separate stats for the contigs with over 1000bp
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
