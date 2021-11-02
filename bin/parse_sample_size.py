#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:36:06 2021

@author: madeline
"""

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Returns strain-specific sample size from a TSV file')

    parser.add_argument('--table', type=str, default=None,
                        help='Multi-strain TSV file generated in workflow that contains num_seqs column')
    parser.add_argument("--filename", default=None,
                        help='Entry in "file" column to parse (usually "strain" + "qc.fasta")')
    return parser.parse_args()


def find_sample_size(table, lineage):
    
    strain_tsv_df = pd.read_csv(table, header=0, delim_whitespace=True, usecols=['file', 'num_seqs'])  
    num_seqs = strain_tsv_df[strain_tsv_df['file']==lineage]['num_seqs'].values

    if len(num_seqs) == 0:
        return "n/a"
    else:
        return num_seqs[0]
    



if __name__ == '__main__':
    
    args = parse_args()
    sample_size = find_sample_size(args.table, args.filename)
    print(sample_size)
    
    
    
'''
Example use case:
sample_size=$(python parse_sample_size.py --lineage B.1.621.qc.fasta --table samples_stats.tsv)
echo $sample_size
'''

