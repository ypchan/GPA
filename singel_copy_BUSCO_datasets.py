#!/usr/bin/env python3
'''
BUSCO_gene_matrix_for_phylogenomic_analysis.py -- select BUSCO genes based on their taxa coverage.

Date:  2021-09-21
Bugs： Any bugs should be reported to yanpengch@qq.com

label_busco_full_table:                                                                                               
GCA_902806535.1    ./GCA_902806535.1_HR_busco/run_ascomycota_odb10/full_table.tsv
GCA_002246955.1    ./GCA_002246955.1_ASM224695v1_busco/run_ascomycota_odb10/full_table.tsv

--busco_desc:
262829at4890    Proteasome subunit alpha type                           https://www.orthodb.org/v10?query=262829at4890
331536at4890    Mediator of RNA polymerase II transcription subunit 6   https://www.orthodb.org/v10?query=331536at4890
'''
import os
import sys
import argparse
import fileinput
import textwrap
import tqdm
import pandas as pd

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('label_busco_full_table',
                    metavar='<label_full_table_path.txt>',
                    type=str,
                    help='label and path list of BUSCO assessment result files: full_table.tsv')

parser.add_argument('-B', '--busco_desc',
                    metavar='<BUSCO_gene_description.txt>',
                    type=str,
                    required=True,
                    help='table of BUSCO gene description including BUSCO ID, Ortho DB url and function')

parser.add_argument('-o', '--out_matrix',
                    metavar='<out_matrix.tsv>',
                    type=str,
                    default='busco_full_matrix.tsv',
                    help='specify the output filename. Default: busco_full_matrix.tsv')

parser.add_argument('-t', '--taxa_coverage',
                    metavar='<int>',
                    type=int,
                    choices=range(1, 101),
                    default=80,
                    help='taxa coverage = (No. of Taxa having this BUSCO gene / No. all Taxa) * 100. Default: 80')

parser.add_argument('-O', '--out_dir',
                    metavar='<out_directory>',
                    type=str,
                    required=False,
                    default='single_copy_BUSCO_dataset',
                    help='output directory name')


args = parser.parse_args()


def read_busco_results(args_busco_full_table):
    '''read all busco results into python ldictionary
    '''
    label_busco_result_dict = {}
    for line in fileinput.input(args_busco_full_table):
        label, full_table_path = line.strip('\n').split()
        label_busco_result_dict[label] = full_table_path
    return label_busco_result_dict


def read_busco_description(args_busco_desc):
    '''read busco gene description into python dictiondary
    '''
    busco_desc_dict = {}
    with open(args.busco_desc, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            line_lst = line.split('\t')
            busco_id = line_lst[0]
            url = line_lst[2]
            desc = line_lst[1]
            busco_desc_dict[busco_id] = [url, desc]
    return busco_desc_dict


def construct_matrix(busco_desc_dict, label_busco_result_dict):
    '''construct busco gene matrix, index represents busco gene id, colname represents taxa name.
    '''
    taxa_label_lst = list(label_busco_result_dict.keys())
    col_name_lst = ['OrthoDB_URL', 'Desc', 'No_taxa',
                    'Coverage%'] + taxa_label_lst

    row_name_lst = list(busco_desc_dict.keys())
    df = pd.DataFrame(columns=col_name_lst, index=row_name_lst)

    df['OrthoDB_URL'] = [lst[0] for k, lst in busco_desc_dict.items()]
    df['Desc'] = [lst[1] for k, lst in busco_desc_dict.items()]

    for label, full_table_path in label_busco_result_dict.items():
        buscoid_status = {}
        with open(full_table_path, 'rt') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue
                line_lst = line.split('\t')
                busco_id = line_lst[0]
                status = line_lst[1]
                if busco_id not in buscoid_status:
                    buscoid_status[busco_id] = status
        df[label] = [buscoid_status[k] for k in df.index.tolist()]
    return df


def filter_matrix(df, args_taxa_coverage, args_out_matrix):
    '''according given threshold of taxa coverage to remove unsatisfied BUSCO genes
    '''
    num_taxa = df.shape[1] - 4
    df['No_taxa'] = (df.iloc[:, 4:] == 'Complete').sum(axis=1)
    coverage_lst = df['No_taxa'] / num_taxa * 100
    df['Coverage%'] = [round(coverage, 2) for coverage in coverage_lst]

    df_filtered = df[df['Coverage%'] >= args_taxa_coverage]
    num_busco_filtered = df.shape[0] - df_filtered.shape[0]
    print(
        f'[INFO] Number of BUSCO genes with low taxa coverage: {num_busco_filtered}\n', file=sys.stdout, flush=True)

    df.to_csv(args_out_matrix, sep='\t', index=True)
    return df_filtered


def read_fasta(fasta_file):
    '''read fasta file into python list
    '''
    seq_lst = []
    with open(fasta_file, 'rt') as infh:
        for line in infh:
            line.rstrip('\n')
            if line.startswith('>'):
                continue
            seq_lst.append(line)
    seq = ''.join(seq_lst)
    return seq


def construct_busco_gene_dataset(label_busco_result_dict, df_filtered):
    '''obtain the path of single_copy_busco_sequences of each BUSCO results
    '''
    single_copy_busco_dict = {}

    busco_id_lst = df_filtered.index.tolist()
    pbar = tqdm.tqdm(label_busco_result_dict.keys())
    for label in pbar:
        pbar.set_description(f"Reading {label}")
        full_table_path = label_busco_result_dict[label]
        single_copy_sequence_path = os.path.dirname(
            full_table_path) + '/busco_sequences/single_copy_busco_sequences'

        for busco_id in busco_id_lst:
            if df_filtered.loc[busco_id, label] == 'Complete':
                busco_sequence_path = single_copy_sequence_path + '/' + busco_id + '.faa'
                busco_sequence = read_fasta(busco_sequence_path)
                if busco_id not in single_copy_busco_dict:
                    single_copy_busco_dict[busco_id] = {}
                single_copy_busco_dict[busco_id][label] = busco_sequence.upper()
    print('', file=sys.stdout, flush=True)
    return single_copy_busco_dict


def out_busco_single_copy_dataset(single_copy_busco_dict, args_out_dir):
    '''output single-copy BUSCO protein datasets
    '''
    pbar = tqdm.tqdm(single_copy_busco_dict.keys())
    for busco_id in pbar:
        pbar.set_description(f"Writing {busco_id}")
        taxa_seq_dict = single_copy_busco_dict[busco_id]
        with open(args_out_dir + '/' + busco_id + '.faa', 'wt') as ofh:
            for taxa, seq in taxa_seq_dict.items():
                wrapped_seq = textwrap.fill(seq, width=80)
                ofh.write(f'>{taxa}\n{wrapped_seq}\n')


if __name__ == '__main__':
    label_busco_result_dict = read_busco_results(args.label_busco_full_table)
    busco_desc_dict = read_busco_description(args.busco_desc)
    df = construct_matrix(busco_desc_dict, label_busco_result_dict)
    df_filtered = filter_matrix(df, args.taxa_coverage, args.out_matrix)
    single_copy_busco_dict = construct_busco_gene_dataset(
        label_busco_result_dict, df_filtered)
    out_busco_single_copy_dataset(single_copy_busco_dict, args.out_dir)
    sys.exit('Done!')
