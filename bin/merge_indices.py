#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:00:04 2024

@author: madeline

This script merges multiple mutation indices into one.
If an existing mutation index is provided, the new indices are added to it.
"""
import dask
dask.config.set({'dataframe.query-planning': True})
from dask.distributed import Client
import dask.dataframe as dd
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
    parser.add_argument('--cpus', type=int, default=1,
                        help='Number of CPUs to use')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    original_index = args.original_index
    gvf_indices_list = args.gvf_indices
    index_savefile = args.index_savefile
    cpus = args.cpus

    # Set up Dask distributed client with specified CPUs
    client = Client(n_workers=cpus)
    print(f"Dask client initialized with {cpus} workers. Dashboard available at {client.dashboard_link}")

    # if an original index is provided, merge the new indices with it
    if original_index is not None:
        indices_to_merge = [original_index] + gvf_indices_list
    else:
        indices_to_merge = gvf_indices_list

    # read all indices into dask df
    ddf = dd.read_csv(
        indices_to_merge,
        sep='\t',
        dtype={
            'alias': 'object',
            'alias_protein': 'object',
            'hgvs_aa_mutation': 'object',
            'hgvs_alias': 'object',
            'protein_name': 'object',
            'gene_name': 'object',
            'gene_symbol': 'object',
            'mat_pep': 'object',
            'protein_symbol': 'object',
            'pos': 'int64',  # Assuming 'pos' is an integer column
            'lineages': 'object'
        }
    )
    # fillna to make groupby() work
    ddf = ddf.fillna('n/a')
    # specify which columns to group by
    group_cols = [x for x in ddf.columns if x != 'lineages']
    # do groupby operation
    ddf = ddf.groupby(by=group_cols)['lineages'].apply(','.join).reset_index()
    # sort by 'pos'
    ddf = ddf.sort_values("pos")

    # save ddf as a single TSV
    ddf.to_csv(index_savefile, single_file=True, sep='\t', index=False)

    print("Index merging completed. File saved at:", index_savefile)