#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 15:03:49 2023

@author: madeline

If args.ingvf, the attributes completed by this script are: 
['multi_aa_name', 'multiaa_comb_mutation']

This script is specific to Pokay functional annotations for
SARS-CoV-2, to harmonize Pokay mutation names with VCF mutation names.
"""

import argparse
import pandas as pd
import numpy as np
from functions import separate_attributes, rejoin_attributes, \
    split_names, unnest_multi
from functions import empty_attributes, gvf_columns

def parse_args():
    parser = argparse.ArgumentParser(
        description='Splits composite mutation names in specified files')
    
    parser.add_argument('--functional_annotations', type=str,
                        default=None, help='TSV file of functional '
                                           'annotations')
    parser.add_argument('--out_functions', type=str,
                        help='Filename for the output tsv file')
    
    parser.add_argument('--names_to_split', type=str,
                        default=None,
                        help='.tsv of multi-aa mutation names to '
                             'split up into individual aa names')
    return parser.parse_args()

    


if __name__ == '__main__':

    args = parse_args()
    
    # split names in functional annotations file
    
    # read in functional annotations file
    df = pd.read_csv(args.functional_annotations, sep='\t', header=0)
    df = df.drop_duplicates()
    
    dataFrame_cols = df.columns.tolist() + ['comb_mutation']

    # remove any leading/trailing spaces
    for column in df.columns:
        df[column] = df[column].astype(str).str.strip()

    # create "comb_mutation" column with the lists from 'original mutation description'
    df['comb_mutation'] = df['original mutation description']
    
    # unnest all these columns
    to_unnest = ['original mutation description', 'nucleotide position', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias']
    df['amino acid mutation alias'] = df['amino acid mutation alias'].replace('nan', np.nan, regex=False)

    for column in to_unnest:
        # convert not-null values to lists
        df.loc[df[column].notna(), column] = df[column].astype(str).str.split(',')
        #convert null values to empty lists
        df[column] = df[column].fillna("").apply(list)
        #df.loc[df[column].isnull(),[column]] = df.loc[df[column].isnull(),column].apply(lambda x: [])
        # make temporary column for lengths of lists
        len_count = column + ' length'
        df.loc[df[column].notna(), len_count] = df[column].str.len()

    # make sure all columns have the same number of elements; if they don't, throw an error
    for length_column in ['original mutation description length', 'nucleotide mutation length', 'amino acid mutation length', 'amino acid mutation alias length']:
        if not (df[length_column].equals(df['nucleotide position length'])):
            column = length_column.replace(' length', '')
            #print(column + " and nucleotide position columns are of unequal length")   
            # now fix it by padding the shorter list with nans
            # calculate how many elements to add to each list
            df['to_add'] = df['nucleotide position length'] - df[length_column]
            # convert columns to lists of lists
            A_list = df[column].tolist() # column that needs padding with NaNs, represented as a nested list
            add_list = df['to_add'].tolist() # list of integers representing how many elements to add to each list
            na_list = [np.nan]*len(add_list) # flat list of NaNs
            # create new list padded with NaNs
            nas_to_add_list = [[a]*b for a,b in zip(na_list, add_list)] # nested list of NaNs to add
            padded_list = [a + b for a,b in zip(A_list, nas_to_add_list)] # column to explode padded with NaNs
            # add padded column to df
            df[column] = padded_list
            #print("Fixed.")

    # explode list columns
    df = unnest_multi(df, to_unnest, reset_index=True)

    # keep only unique names in 'comb_mutation', no repeats
    comb_mutation_nested_list = df['comb_mutation'].str.split(',').tolist()
    unique_nested_list = [list(set(ls)) for ls in comb_mutation_nested_list]
    df['comb_mutation'] = unique_nested_list
    df['comb_mutation'] = df['comb_mutation'].str.join(",")

    # remove the one 'original mutation description' from the 'comb_mutation' entry
    df['comb_mutation'] = df.apply(lambda row : row['comb_mutation'].replace(str(row['original mutation description']), ''), axis=1)
    df['comb_mutation'] = df['comb_mutation'].str.strip(',')
    df['comb_mutation'] = df['comb_mutation'].str.replace(',,', ',')

    # data cleaning: fix parsing error
    #df['comb_mutation'] = df['comb_mutation'].str.replace(
    #    "B.1.617.2\\tT19R", "T19R", regex=False)
    # data cleaning: remove quotation marks, spaces
    
    # if names_to_split tsv is given, use it to split up the multi-amino acid names
    if args.names_to_split != None:
        # load names_to_split file
        names_to_split_df = pd.read_csv(args.names_to_split, sep='\t', header=0)
        # remove spaces and quotation marks from 'split_into' column
        for x in [' ', "'"]:
            names_to_split_df['split_into'] = names_to_split_df[
                'split_into'].str.replace(x, "")
        # make names_to_split_df into a dictionary
        split_dict = dict(zip(names_to_split_df.name, names_to_split_df.split_into))
        # do str.replace on 'comb_mutation' and 'original mutation description' columns
        for name in split_dict.keys():
            df['comb_mutation'] = df['comb_mutation'].str.replace(
                name, split_dict[name])
            df['original mutation description'] = df['original mutation description'].str.replace(
                name, split_dict[name])
        # unnest "mutation" column so each mutation name gets its own row
        df['original mutation description'] = df['original mutation description'].str.split(',')
        df = unnest_multi(df, ['original mutation description'], reset_index=True)
    
    # reset df columns
    df = df[dataFrame_cols]
    
    # fix column formats
    df['PMID'] = df['PMID'].str.replace('nan', '', regex=False)
    df['PMID'] = df['PMID'].str.replace('.0', '', regex=False)
    df['publication year'] = df['publication year'].str.replace('.0', '', regex=False)

    # drop duplicate rows
    df = df.drop_duplicates()

    # save modified file
    filepath = args.out_functions
    print("Saved as: ", filepath)
    df.to_csv(filepath, sep='\t', index=False, header=True)

