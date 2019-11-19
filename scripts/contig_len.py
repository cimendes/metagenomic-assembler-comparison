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

import sys
from itertools import groupby
import matplotlib as mpl
## agg backend is used to create plot as a .png file
mpl.use('agg')
import matplotlib.pyplot as plt


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)


def get_n50(contig_lengths):
    contig_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return n50

def main():
    try:
        assembly_filename = sys.argv[1:]

    # If no arguments were given, just print the header line.
    except IndexError as e:
        print(e, file=sys.stderr)
        sys.exit(0)

    data_to_plot = []
    legends = []

    print(','.join(['Assembler', 'n Contigs', 'total bp', 'max contig size', 'n50', '% contigs over 1000bp']))

    for file in assembly_filename:

        contig_len = []
        total_bp = 0
        good_contigs = 0

        fasta = fasta_iter(file)
        for header, seq in fasta:
            if len(seq) > 1000:
                good_contigs += 1
            contig_len.append(len(seq))
            total_bp += len(seq)

        data_to_plot.append(contig_len)
        legends.append(file.split('.')[0].rsplit('_')[1])

        print(','.join([file.split('.')[0].rsplit('_')[1], f'{len(contig_len):.4f}', f'{total_bp}',
                        f'{max(contig_len):.4f}', f'{get_n50(contig_len):.4f}', f'{good_contigs/len(contig_len):.4f}']))

    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot, )

    ax.set_xticklabels(legends)

    # Save the figure
    fig.savefig('fig1.png', bbox_inches='tight')


if __name__ == '__main__':
    main()
