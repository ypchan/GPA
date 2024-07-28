#!/usr/bin/env python3

import os
import sys
import tqdm
import argparse
import fileinput

'''
concatenate_BUSCO_gene_msa.py -- concatenate trimmmed BUSCO gene alignment files.
Date:  2021-09-21
Bugs： Any bugs should be reported to yanpengch@qq.com
'''


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('alginment_lst',
                    metavar='<trimmed_alignment.list>',
                    type=str,
                    help='a list of trimmed BUSCO gene alignments. Stdin is allowed.')

parser.add_argument('-o', '--out',
                    metavar='<super_BUSCO_matrix.faa>',
                    type=str,
                    default='super_BUSCO_matrix.faa',
                    help='output filename. Default: super_BUSCO_matrix.faa')

args = parser.parse_args()


def read_all_alignments(fasta_file_list):
    '''read fasta into python dictionary
    '''
    all_alignment_dict = {}
    pbar = tqdm.tqdm([line for line in fileinput.input(fasta_file_list)])
    for line in pbar:
        alignment = line.rstrip('\n')
        pbar.set_description(f"Reading {alignment}")
        alignment_basename = os.path.basename(alignment).split('.')[0]
        all_alignment_dict[alignment_basename] = {}
        with open(alignment, 'rt') as infh:
            for line in infh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    fa_id = line.lstrip('>')
                    all_alignment_dict[alignment_basename][fa_id] = []
                else:
                    all_alignment_dict[alignment_basename][fa_id].append(line)

    for alignment_basename, fa_dict in all_alignment_dict.items():
        fd_dict = {k: ''.join(v) for k, v in fa_dict.items()}
        all_alignment_dict[alignment_basename] = fd_dict
    return all_alignment_dict


def alignment_length_and_taxa_list(all_alignment_dict):
    '''obtain the length of each alignment file
    '''
    alignement_len_dict = {}
    taxa_lst = []
    for alignment_basename, fa_dict in all_alignment_dict.items():
        n = 0
        for k, v in fa_dict.items():
            n += 1
            if n == 1:
                alignement_len_dict[alignment_basename] = len(v)
            taxa_lst.append(k)

    taxa_lst = list(set(taxa_lst))
    return alignement_len_dict, taxa_lst


def placeholder_for_missing_taxa(alignement_len_dict, taxa_lst, all_alignment_dict):
    '''add placeholder '---'for missing taxa
    '''
    for alignment_basename, fa_dict in all_alignment_dict.items():
        need_add_placeholder_lst = list(set(taxa_lst) - set(fa_dict.keys()))
        for taxa in need_add_placeholder_lst:
            all_alignment_dict[alignment_basename][taxa] = '-' * \
                alignement_len_dict[alignment_basename]
    return all_alignment_dict


def concatenate_alignment(all_alignment_dict, taxa_lst):
    '''concatenate all alignments by identical taxa label
    '''
    concatenated_dict = {}
    for taxa in taxa_lst:
        concatenated_dict[taxa] = []
        for _, fa_dict in all_alignment_dict.items():
            sequence = fa_dict[taxa]
            concatenated_dict[taxa].append(sequence)

    concatenated_dict = {k: ''.join(v) for k, v in concatenated_dict.items()}
    return concatenated_dict

def out_super_matrix(args_out, concatenated_dict):
    '''output concatenated super BUSCO matrix
    '''
    with open(args_out, 'wt') as ofh:
        for k, v in concatenated_dict.items():
            ofh.write(f'>{k}\n{v}\n')


if __name__ == '__main__':
    all_alignment_dict = read_all_alignments(args.alginment_lst)
    alignement_len_dict, taxa_lst = alignment_length_and_taxa_list(
        all_alignment_dict)
    all_alignment_dict = placeholder_for_missing_taxa(
        alignement_len_dict, taxa_lst, all_alignment_dict)
    concatenated_dict = concatenate_alignment(all_alignment_dict, taxa_lst)

    out_super_matrix(args.out, concatenated_dict)
    print('Done', file=sys.stdout, flush=True)
    sys.exit()
