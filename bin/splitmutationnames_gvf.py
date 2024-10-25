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
    parser.add_argument('--ingvf', type=str, default=None,
                        help='Path to the GVF file output of vcf2gvf.py')
    parser.add_argument('--outgvf', type=str,
                        help='Filename for the output gvf file')
    parser.add_argument('--names_to_split', type=str,
                        default=None,
                        help='.tsv of multi-aa mutation names to '
                             'split up into individual aa names')
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    
    # split names in gvf file
    
    # read in gvf file
    gvf = pd.read_csv(args.ingvf, sep='\t', names=gvf_columns, index_col=False)

    # remove pragmas and original header row
    pragmas = gvf[gvf['#seqid'].astype(str).str.contains("##")]
    pragmas.columns = range(9)
    pragmas = pragmas.fillna('')
    gvf = gvf[~gvf['#seqid'].astype(str).str.contains("#")]

    # expand #attributes into columns to edit separately
    gvf = separate_attributes(gvf)

    # if names_to_split tsv is given, use it to split up the multi-amino acid names
    if args.names_to_split != None:
        # split names in "Names" attribute into separate rows
        gvf = split_names(args.names_to_split, gvf, col_to_split='Name')
    
    # rename IDs: rows with the same entry in 'Name'
    # get the same ID
    gvf['ID'] = 'ID_' + gvf.groupby('Name', sort=False).ngroup().astype(str)
    
    # merge attributes back into a single column
    gvf = rejoin_attributes(gvf, empty_attributes)

    # discard temporary columns
    gvf = gvf[gvf_columns]
    
    # add pragmas to gvf
    # columns are now 0, 1, ...
    final_gvf = pd.DataFrame(np.vstack([gvf.columns, gvf]))
    final_gvf = pragmas.append(final_gvf)
    
    # save modified file to .gvf
    filepath = args.outgvf
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)


