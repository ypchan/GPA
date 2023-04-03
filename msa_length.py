#!/usr/bin/env python3
'''
msa_length.py -- get alignment length.

DATE:  2021-09-25
BUGS： Any bugs should be reported to yanpengch@qq.com
'''
import sys
import argparse

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 prog='msa_length.py')
parser.add_argument('msa',
                    nargs='+',
                    type=str,
                    help = 'one or multiple msa files')
args = parser.parse_args()

BUSCO_ID_2_length_dict = {}
for BUSCO_ID in args.msa:
    with open(BUSCO_ID) as fh:
        n = 0
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                n += 1
                fa_seq_lst = []
            else:
                fa_seq_lst.append(line)
            if n > 1:
                 break
            BUSCO_ID_2_length_dict[BUSCO_ID] = len(''.join(fa_seq_lst))

for k,v in BUSCO_ID_2_length_dict.items():
    print(k,v,sep='\t', file=sys.stdout, flush=True)
sys.exit(0)
