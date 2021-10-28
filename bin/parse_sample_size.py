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
    parser.add_argument("--lineage", default='n/a',
                        help='lineage')
    return parser.parse_args()


def find_sample_size(strain_tsv, strain):
    
    strain_tsv_df = pd.read_csv(strain_tsv, header=0, delim_whitespace=True, usecols=['file', 'num_seqs'])  
    num_seqs = strain_tsv_df[strain_tsv_df['file'].str.contains(strain)]['num_seqs'].values[0]

    return num_seqs



if __name__ == '__main__':
    
    args = parse_args()
    
    sample_size = find_sample_size(args.table, args.lineage)
    
    print(sample_size)
    
    
    
'''
In bash:
sample_size=$(python parse_sample_size.py --lineage B.1.621 --table samples_stats.tsv)
echo $sample_size
'''

