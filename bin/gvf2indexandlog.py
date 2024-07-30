#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on Fri Feb 16 12:00:04 2024

@author: madeline

Given one GVF, creates a mutation index TSV and a logfile TSV.
'''

import pandas as pd
import numpy as np
import argparse
from functions import separate_attributes


def parse_args():
    parser = argparse.ArgumentParser(
        description='Creates an index TSV and a log TSV from one GVF')
    parser.add_argument('--gvf_file', type=str, default=None,
                        help='Path to one GVF file to process')
    parser.add_argument('--index_savefile', type=str,
                        default=None, help='Filename to save updated index of all mutations to')
    parser.add_argument('--log_savefile', type=str,
                        default=None, help='Filename to save log to')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    gvf_file = args.gvf_file
    index_savefile = args.index_savefile
    log_savefile = args.log_savefile


    # read in gvf
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                '#score', '#strand', '#phase', '#attributes']
    gvf = pd.read_csv(gvf_file, sep='\t', names=gvf_columns, usecols=['#start', '#seqid', '#attributes'])

    # remove pragmas and original header
    gvf = gvf[~gvf['#seqid'].str.contains("#")]

    # separate attributes
    gvf = separate_attributes(gvf)

    # create index from GVF
    # make empty index df
    index_cols=['pos', 'mutation', 'hgvs_aa_mutation', 'hgvs_nt_mutation', 'gene_name', 'gene_symbol', 'protein_name', 'protein_symbol', 'alias', 'hgvs_alias', 'alias_protein', 'Pokay_annotation', 'lineages']
    index = pd.DataFrame(np.empty((gvf.shape[0], len(index_cols))), columns=index_cols)
    # populate index df with gvf info
    index['pos'] = gvf['#start']
    index['mutation'] = gvf['Name'].str.replace("p.", "", regex=False)
    index['hgvs_aa_mutation'] = gvf['hgvs_aa']
    index['hgvs_nt_mutation'] = gvf['hgvs_nt']
    index['alias'] = gvf['alias']
    index['hgvs_alias'] = gvf['hgvs_alias']
    index['alias_protein'] = 'n/a'
    index.loc[index['alias']!='n/a', 'alias_protein'] = gvf['mat_pep']
    index['gene_name'] = gvf['gene_name']
    index['gene_symbol'] = gvf['gene_symbol']
    index['protein_name'] = gvf['protein_name']
    index['protein_symbol'] = gvf['protein_symbol']
    index['Pokay_annotation'] = gvf["function_description"].notna()
    index['lineages'] = gvf['viral_lineage']
    # tidying
    index = index.drop_duplicates()
    index = index.dropna(axis=0)
    # save index
    index.to_csv(index_savefile, sep='\t', header=True, index=False)

    # create log from index
    log = index.copy()
    # fill in 'new_mutations' column like: "gene:mutation"
    log['new_mutations'] = log["gene_symbol"] + ":" + log["mutation"]
    # for orf1ab mutations, fill in 'new_mutations' column like: "gene:mutation / nsp:alias"
    log.loc[log['alias']!='n/a', 'new_mutations'] = log['new_mutations'] + " / " + log["alias_protein"] + ":" + log["alias"]
    # drop duplicates (there shouldn't be any)
    log = log[['pos', 'new_mutations', 'lineages']].drop_duplicates()
    # drop any NaN rows
    log = log[log['pos'].notna()]
    # ensure pos is integer type
    log['pos'] = log['pos'].astype(int)
    # save log with column headers: ['pos', 'new_mutations', 'lineages']
    log.to_csv(log_savefile, sep='\t', header=True, index=False)

