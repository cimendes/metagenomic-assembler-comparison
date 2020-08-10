#!/usr/bin/env python3
"""

"""
import sys
from itertools import groupby
import glob
import os
import re
import fnmatch

#import commonly used functions from utils.py
import utils

REFERENCE_SEQUENCES = os.path.join(os.path.dirname(__file__),
                                   '..', 'data', 'references', 'Zymos_Genomes_triple_chromosomes.fasta')


def get_C90(alignment_lengths, ref_len):
    """

    :param alignment_lengths: list with length of mapped contigs for the reference
    :param ref_len: int with the expected reference length
    :return: int with the number of contigs that represent
    """
    sorted_lengths = sorted(alignment_lengths, reverse=True)  # from longest to shortest
    target_length = ref_len * 0.9

    length_so_far = 0
    c90 = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            c90 += 1
    return c90


def get_covered_bases(covered_bases_list, ref_len):
    """
    Get ration of referee lengths (adjusted for triple reference) covered by mapping contigs
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


def get_expanded_cigar(cigar):
    """

    :param cigar: string with cigar values
    :return: string with expanded cigar values for alignment
    """
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def get_lowest_window_identity(cigar, window_size):
    """

    :param cigar: string with alignment cigar
    :param window_size: int with window size
    :return: float with lowest identity value for mapping contigs
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


def get_alignment_stats(paf_filename, ref_name, ref_length):
    """

    :param paf_filename:
    :param ref_name:
    :param ref_length:
    :return:
    """

    # Tracks the longest single alignment, in terms of the reference bases.
    longest_alignment = 0
    identity = 0
    longest_alignment_cigar = ''

    covered_bases = []

    # to calculate na50
    contigs = {}
    alignment_lengths = []

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                if parts[0] not in contigs.keys():
                    contigs[parts[0]] = int(parts[1])
                start, end = int(parts[7]), int(parts[8])
                alignment_lengths.append(end - start)
                matching_bases, total_bases = int(parts[9]), int(parts[10])
                cigar = [x for x in parts if x.startswith('cg:Z:')][0][5:]
                if end - start > longest_alignment:
                    longest_alignment = end - start
                    identity = matching_bases / total_bases
                    longest_alignment_cigar = cigar
                longest_alignment = max(longest_alignment, end - start)
                covered_bases.append([start, end])

    contiguity = longest_alignment / ref_length
    lowest_identity = get_lowest_window_identity(longest_alignment_cigar, 1000)

    coverage = get_covered_bases(covered_bases, ref_length)

    return contiguity, coverage, lowest_identity


def parse_paf_files(df, mappings):
    """
    Parses fasta, paf files references and returns info in dataframe.
    :param df: pandas DataFrame with assembly stats
    :param mappings: list of paf files
    :return: pandas Dataframe with columns in COLUMNS
    """

    for assembler in sorted(df['Assembler'].unique()):

        print('\n\n------' + assembler + '------\n')

        # Assembler df
        df_assembler = df[df['Assembler'] == assembler]

        # iterator for reference files (length is needed)
        references = (x[1] for x in groupby(open(REFERENCE_SEQUENCES, "r"), lambda line: line[0] == ">"))

        paf_file = fnmatch.filter(mappings, '*_' + assembler + '.*')[0]

        print(','.join(["reference", "reference length", "contiguity", "identity", "lowest identity",
                        "breadth of coverage", "c90", "aligned contigs", "na50", "aligned bp"]))

        for header in references:
            header_str = header.__next__()[1:].strip().split()[0]
            seq = "".join(s.strip() for s in references.__next__())

            df_assembler_reference = df_assembler[df_assembler['Mapped'] == header_str]

            mapped_contigs = df_assembler_reference['Contig Len'].astype('int').tolist()

            na50 = utils.get_N50(mapped_contigs)
            c90 = get_C90(mapped_contigs, len(seq)/3)  # adjust for triple reference

            contiguity, coverage, lowest_identity = get_alignment_stats(paf_file, header_str, len(seq)/3)

            mapped_contigs = df_assembler_reference['Contig Len']

            print(','.join([header_str, f'{len(seq)/3}', f'{contiguity:.2f}', f'{1:.2f}',
                            f'{lowest_identity:.2f}', f'{coverage:.2f}', f'{c90}',
                            f'{len(mapped_contigs)}',f'{na50}', f'{sum(mapped_contigs)}']))


def add_matching_ref(df, mappings):
    """
    For each contig in the df, adds the correspondent reference if the contig is mapped. Drops the unmapped contigs.
    :param df: Pandas Dataframe with stats for each contig
    :param mappings: list of paf files
    :return: Pandas Dataframe with stats for each contig with reference info instead of 'Mapped' and rows with
    unmapped contigs removed
    """
    for assembler in sorted(df['Assembler'].unique()):
        paf_file = fnmatch.filter(mappings, '*_' + assembler + '.*')[0]
        mapped_contigs = utils.get_mapped_contigs_with_ref(paf_file)  # dictionary with contigs
        
        for contig in df['Contig'][(df['Mapped'] == 'Mapped') & (df['Assembler'] == assembler)]:
            #get index.
            row_index = df[(df['Contig'] == contig) & (df['Mapped'] == 'Mapped') & (df['Assembler'] == assembler)]\
                .index.item()
            df.loc[row_index,'Mapped'] = mapped_contigs[contig] # update with reference

    #remove unmapped contigs from dataframe
    df = df.drop(df[df.Mapped == 'Unmapped'].index)
    return df


def main():
    try:
        assemblies = glob.glob(sys.argv[1] + '/*.fasta')
        mappings = glob.glob(sys.argv[2] + '/*.paf')

    except IndexError as e:
        print(e, "files not found")
        sys.exit(0)

    # Dataframe with assembly info
    df = utils.parse_assemblies(assemblies, mappings)

    # Add correspondent reference to each dataframe contig
    df = add_matching_ref(df, mappings)

    #for assembler in df['Assembler'].unique():
    #    print(assembler)
    #    print(len(df['Contig'][df['Assembler'] == assembler]))

    parse_paf_files(df, mappings)


if __name__ == '__main__':
    main()
