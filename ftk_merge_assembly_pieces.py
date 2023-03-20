#!/usr/bin/env/python

'''
ftk_merge_assembly_pieces.py -- merge pieces of a genome into a single FASTA file.

Date: 2022-10-28
Bugs: Any bugs should be reported to yanpengch@qq.com

Usage:
    cat file_lst | ftk_merge_assembly_pieces.py - -o output_directory

Input file example:
GCA_021398005.1 ./GCA_021398005.1/ncbi_dataset/data/GCA_021398005.1/GCA_021398005.fna
GCA_021436885.1 ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr1.fna
GCA_021436885.1 ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr2.fna
GCA_021436885.1 ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr3.fna
GCA_021436885.1 ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr4.fna
GCA_021436885.1 ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/unplaced.scaf.fna
'''

import os
import sys
import argparse
import fileinput

from tqdm import tqdm

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('list',
                    type=str,
                    help='a two-column list, stdin is allowed. Input list is examplifed as above')
parser.add_argument('-o', '--outdir',
                    type=str,
                    default='.',
                    help='output directory (default: .). the first column is treated as the filename with .fna as the suffix')
args = parser.parse_args()

assblies_dict = dict()

with fileinput.input(files=args.list, mode='r') as infh:
    for line in infh:
        assembly, piece_path = line.rstrip('\n').split('\t')
        if assembly not in assblies_dict:
            assblies_dict[assembly] = [piece_path]
        else:
            assblies_dict[assembly].append(piece_path)

pbar = tqdm(list(assblies_dict.keys()))
for assembly in pbar:
    pbar.set_description("Merging %s" % assembly)
    piece_lst = assblies_dict[assembly]
    output_file = args.outdir + '/' + assembly + '.fna'
    command = f"cat {' '.join(piece_lst)} > {output_file}"
    res_code = os.system(command)
    if res_code != 0:
        sys.exit(f'Error: failed to call {command}')

sys.exit(0)
