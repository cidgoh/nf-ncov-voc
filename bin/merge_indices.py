#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:00:04 2024

@author: madeline

This script merges multiple mutation indices into one.
If an existing mutation index is provided, the new indices are added to it.
"""

import dask.dataframe as dd
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merges mutation indices into one.')
    parser.add_argument('--original_index', type=str, default=None,
                        help='Path to an existing TSV index of all mutations, \
                            to merge new indices to')
    parser.add_argument('--gvf_indices', type=str, default=None,
                        nargs='*', help='Paths to GVF indices to merge')
    parser.add_argument('--index_savefile', type=str,
                        default=None, help='TSV filename to save updated \
                            index of all mutations to')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    original_index = args.original_index
    gvf_indices_list = args.gvf_indices
    index_savefile = args.index_savefile

    # if an original index is provided, merge the new indices with it
    if original_index!=None:
        indices_to_merge = list(original_index) + gvf_indices_list
    # if an original index is not provided, make an entirely new index
    else:
        indices_to_merge = gvf_indices_list
    
    # read all indices into dask df
    ddf = dd.read_csv(indices_to_merge, sep='\t') 
    # fillna to make groupby() work
    ddf = ddf.fillna('n/a')
    # groupby all columns except 'lineages'...
    group_cols = [x for x in ddf.columns if x!='lineages']
    # ...and merge 'lineages' column into a comma-separated list
    d = {**dict.fromkeys(ddf.columns, 'first'), 'lineages': 'list'}
    '''
    # specify meta df: this throws an error, but the inferred meta works fine
    schema = {**dict.fromkeys(ddf.columns, 'str'), 'pos':'int'}
    meta_df = pd.DataFrame(columns=schema.keys()).astype(schema)
    '''
    # do groupby
    ddf = ddf.groupby(by=group_cols)['lineages'].apply(','.join).reset_index()
    # sort by 'pos'
    ddf = ddf.sort_values("pos")

    # save ddf as a single TSV
    ddf.to_csv(index_savefile, single_file=True, sep='\t', index=False)