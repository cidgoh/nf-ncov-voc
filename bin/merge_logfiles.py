#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:00:04 2024

@author: madeline

This script merges multiple mutation logs into one.
If an existing mutation log is provided, the new logs are added to it.

--log_header is the path to a text file of the form:

Total sequences:	508910
New sequences:	171
New lineages:	['HH.1.1', 'FT.3.1.1']
"""
import dask
dask.config.set({'dataframe.query-planning-warning': False}) # this recommends dask-exp install which we don't want to use as it's unstable
import dask.dataframe as dd
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merges mutation logs into one.')
    parser.add_argument('--original_index', type=str, default=None,
                        help='Path to an existing TSV index of all mutations, \
                            used to determine which mutations to log as new')
    parser.add_argument('--log_header', type=str, default=None,
                        help='Path to a text file containing the log header')
    parser.add_argument('--gvf_logs', type=str, default=None,
                        nargs='*', help='Paths to GVF logs to merge')
    parser.add_argument('--log_savefile', type=str,
                        default=None, help='TSV filename to save updated \
                            log of all mutations to')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    original_index = args.original_index
    log_header = args.log_header
    gvf_logs_list = args.gvf_logs
    log_savefile = args.log_savefile

    # read all logs into dask df
    ddf = dd.read_csv(gvf_logs_list, sep='\t') 
    # fillna to make groupby() work
    ddf = ddf.fillna('n/a')
    # groupby all columns except 'lineages'...
    group_cols = [x for x in ddf.columns if x!='lineages']
    ddf = ddf.groupby(by=group_cols)['lineages'].apply(','.join).reset_index()

    # sort by 'pos'
    ddf = ddf.sort_values("pos")

    # if a mutation index is provided, remove
    # previously-identified mutations from the logfile
    if original_index!=None:

        # convert original_index to match the logfile columns
        # read index into pandas df
        og_index = pd.read_csv(original_index, sep='\t').fillna('n/a')
        # fill in 'new_mutations' column like: "gene:mutation"
        og_index['new_mutations'] = og_index["gene"] + ":" + og_index["mutation"]
        # for orf1ab mutations, fill in 'new_mutations' column like: "gene:mutation / nsp:alias"
        og_index.loc[og_index['alias']!='n/a', 'new_mutations'] = og_index['new_mutations'] + " / " + og_index["alias_protein"] + ":" + og_index["alias"]
        # drop duplicates (there shouldn't be any)
        og_index = og_index[['pos', 'new_mutations', 'lineages']].drop_duplicates()
        # drop any NaN rows
        og_index = og_index[og_index['pos'].notna()]
        # ensure pos is integer type
        og_index['pos'] = og_index['pos'].astype(int)
        # rename 'lineages' column entries
        og_index['lineages'] = 'OG_index'
        # concatenate original index to logfiles, with the og index first
        ddf = dd.concat([og_index, ddf], axis=0)

        # delete all duplicated rows, based on 'pos' and 'new_mutations'
        # by default, the first ('OG_index') duplicate is kept here
        ddf = ddf.drop_duplicates(subset=['pos', 'new_mutations']) 
        # remove 'OG_index' rows, thereby removing all duplicates
        ddf = ddf[ddf['lineages']!='OG_index']
        # sort by 'pos'
        ddf = ddf.sort_values("pos")
        ddf = ddf[['pos', 'new_mutations', 'lineages']]
        ddf['pos'] = ddf['pos'].astype(int)
        
    # print number of new mutations in >=5 lineages
    num_lineages = ddf['lineages'].compute()
    num_lineages = num_lineages.str.split(',').str.len()
    print("New mutations in >=5 lineages: ", (num_lineages>=5).sum())
    
    # add log header
    log_header_df = pd.read_csv(log_header, sep='\t', names=["pos", "new_mutations", "lineages"])
    log_header_df.loc[len(log_header_df)] = ["New mutations:", "", ""]
    ddf = dd.concat([log_header_df, ddf])
    # save ddf as a single TSV
    ddf.to_csv(log_savefile, single_file=True, sep='\t', index=False, header=False)
    
