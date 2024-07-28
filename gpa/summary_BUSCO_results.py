#!/usr/bin/env python3
'''
summary_BUSCO_results.py -- summary multiple busco result files into a table

Date:
    2020-04-14
Bugs：
    Any bugs should be reported to yanpengch@qq.com
Usage:
    find . -name 'short_summary.specific*txt' | summary_BUSCO_results.py - > BUSCO_results_summary.txt
'''
import re
import os
import sys
import fileinput

def parse_busco_result2lst(busco_result_file):
    '''short_summary.specific*:
    # BUSCO version is: 4.1.4
    # The lineage dataset is: ascomycota_odb10 (Creation date: 2020-09-10, number of species: 365, number of BUSCOs: 1706)
    # Summarized benchmarking in BUSCO notation for file SRR10394974_contigs_200.faa
    # BUSCO was run in mode: proteins

            ***** Results: *****

            C:97.9%[S:97.6%,D:0.3%],F:0.7%,M:1.4%,n:1706
            1670    Complete BUSCOs (C)
            1665    Complete and single-copy BUSCOs (S)
            5       Complete and duplicated BUSCOs (D)
            12      Fragmented BUSCOs (F)
            24      Missing BUSCOs (M)
            1706    Total BUSCO groups searched
    '''
    with open(busco_result_file) as buscofh:
        for line in buscofh:
            line = line.strip()
            if line.startswith('# Summarized benchmarking in BUSCO notation for file'):
                genome_label = os.path.basename(line.split(' ')[-1])

            if line.startswith('C:'):
                index_lst = [1,3,5,8,10, 12]
                per_C, per_S, per_D, per_F, per_M, num_Total = [re.split(r'[:\[\],]', line)[pos] for pos in index_lst]
                continue
                
            if 'Complete BUSCOs (C)' in line:
                num_C = line.split()[0]
            if 'Complete and single-copy BUSCOs (S)' in line:
                num_S = line.split()[0]
            if 'Complete and duplicated BUSCOs (D)' in line:
                num_D = line.split()[0]
            if 'Fragmented BUSCOs (F)' in line:
                num_F = line.split()[0]
            if 'Missing BUSCOs (M)' in line:
                num_M = line.split()[0]
    return [genome_label, per_C, per_S, per_D, per_F, per_M, num_C, num_S, num_D, num_F, num_M, num_Total]

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(__doc__, file=sys.stderr, flush=True)
        sys.exit(1)
    if sys.argv[1] in ['-h', '--h', '--help']:
        print(__doc__, file=sys.stderr, flush=True)
        sys.exit(1)

    head_lst = ['genome_label', 'C', 'S', 'D', 'F', 'M', 'Complete BUSCOS (C)', 'Complete and single-copy BUSCOs (S)', 'Complete and duplicated BUSCOs (D)', 'Fragemented BUSCOs (F)', 'Missing BUSCOs (M)', 'Total BUSCO groups searched']
    print('\t'.join(head_lst), file=sys.stdout, flush=True)

   
    with fileinput.input(files=sys.argv[1]) as filefh:
        for file in filefh:
            out_values_lst = parse_busco_result2lst(file.strip('\n'))
            print('\t'.join(out_values_lst), file=sys.stdout, flush=True)
    sys.exit(0)
