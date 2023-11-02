#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 08:25:48 2023

@author: madeline

This script annotates GVF files with the functional annotation.

The attributes completed by this script are: 
["function_category", "function_description", "source",
 "citation", "comb_mutation", "heterozygosity"]

"""
import argparse
import pandas as pd
import numpy as np
from functions import separate_attributes, rejoin_attributes
from functions import empty_attributes, gvf_columns, vcf_columns


def parse_args():
    parser = argparse.ArgumentParser(
        description='Adds functional annotation to a GVF file')
    parser.add_argument('--ingvf', type=str, default=None,
                        help='Path to a GVF file')
    parser.add_argument('--outgvf', type=str,
                        help='Filename for the output GVF file')
    parser.add_argument('--functional_annotations', type=str,
                        default=None, help='TSV file of functional '
                                           'annotations')
    parser.add_argument("--names", type=str, default='n/a',
    			help="Save mutation names without "
                              "functional annotations to "
                              "this .txt filename for "
                              "troubleshooting purposes")
    return parser.parse_args()
    

def add_pokay_annotations(gvf, annotation_file):
    
    # expand #attributes into columns to fill in separately
    gvf = separate_attributes(gvf)
    
    # drop columns that are going to be re-added in the merge
    functional_attributes = ["function_category", "function_description", 
                             "source", "citation", "comb_mutation", 
                             "heterozygosity"]
    gvf = gvf.drop(columns=functional_attributes)

    # load functional annotations spreadsheet
    df = pd.read_csv(annotation_file, sep='\t', header=0)
    # remove any leading/trailing spaces
    for column in df.columns:
        df[column] = df[column].str.strip()

    # merge annotated vcf and functional annotation files by 'Name' and 'alias'
    df = df.rename(columns={"mutation": "Name", "gene": "protein_symbol", "alias":"Pokay_alias"})
    merged_df = pd.merge(df, gvf, on=['Name', 'protein_symbol'], how='right') #, 'alias'

    # data cleaning
    merged_df['comb_mutation'] = merged_df['comb_mutation'].str.replace(
        "B.1.617.2\\tT19R", "T19R", regex=False)
    
    # update ID attribute based on mutation groups
    
    # collect all mutation groups (including reference mutation) in
    # merged_df["mutation_group"], sorted alphabetically
    
    # join columns with commas in between
    merged_df["mutation_group"] = \
        merged_df["Name"].astype(str) + "," + \
        merged_df["comb_mutation"].astype(str) + "," + \
        merged_df['multiaa_comb_mutation'].astype(str)
    # cleaning: remove nans, quotations marks, spaces
    for x in [' ', 'nan', "'", '"']:
        merged_df["mutation_group"] = merged_df[
            "mutation_group"].str.replace(x, '')
    # cleaning: remove extra commas
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.replace(',,', ',')
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.strip(',')
    # sort each row alphabetically
    # extract column into a list (each row is an element)
    unsorted_group_list = merged_df["mutation_group"].tolist()
    # convert each element of the list to a list, making a nested list
    nested_list = [x.split(",") for x in unsorted_group_list]
    # sort sublists alphabetically
    nested_list = [sorted(x) for x in nested_list]
    # return sublists to strings
    sorted_group_list = [",".join(x) for x in nested_list]
    # add sorted lists back to the df
    merged_df["mutation_group"] = pd.Series(sorted_group_list)

    # make another column to check if all members of the group are
    # represented individually in 'Name' (True/False)
    unique_Name_entries = set(merged_df['Name'].tolist())
    group_fully_represented = [set(x).issubset(unique_Name_entries)
                               for x in nested_list]
    merged_df['group_fully_represented'] = group_fully_represented
    # drop rows with mutation group members not found in 'Name',
    # leaving the index unchanged
    merged_df = merged_df[merged_df['group_fully_represented']==True]

    # update 'ID' attribute: now, rows with the same entry
    # in 'mutation_group' get the same ID
    merged_df['ID'] = 'ID_' + merged_df.groupby(
        'mutation_group', sort=False).ngroup().astype(str)
    
    # change semicolons in function descriptions to colons
    merged_df['function_description'] = merged_df[
        'function_description'].str.replace(';', ':')
    
    # change heterozygosity column to True/False
    merged_df['heterozygosity'] = merged_df['heterozygosity'] == \
                                  'heterozygous'
        
    # replace NaNs in df with empty string
    merged_df = merged_df.fillna('')

    # merge attributes back into a single column
    merged_df = rejoin_attributes(merged_df, empty_attributes)

    return merged_df[gvf_columns]


if __name__ == '__main__':

    args = parse_args()
    
    # read in gvf file
    gvf = pd.read_csv(args.ingvf, sep='\t', names=gvf_columns, index_col=False)

    # remove pragmas and original header row
    pragmas = gvf[gvf['#seqid'].astype(str).str.contains("##")]
    pragmas.columns = range(9)
    pragmas = pragmas.fillna('')
    gvf = gvf[~gvf['#seqid'].astype(str).str.contains("#")]

    # add functional annotations
    pokay_annotated_gvf = add_pokay_annotations(gvf, args.functional_annotations)

    # add pragmas to df, then save to .gvf
    # columns are now 0, 1, ...
    final_gvf = pd.DataFrame(np.vstack([pokay_annotated_gvf.columns,
                                            pokay_annotated_gvf]))
    final_gvf = pragmas.append(final_gvf)
    filepath = args.outgvf  # outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)


    # get name troubleshooting report
    if args.names!='n/a':
        # save unmatched names (in vcf/tsv but not in
        # functional_annotations) to a .tsv file
        
        # create mask to find which rows do not have a functional annotation
        notinPokay_mask = pokay_annotated_gvf["#attributes"].str.contains("function_category=;")
        # extract all mutation names from #attributes column
        names = pd.Series(pokay_annotated_gvf["#attributes"].str.findall('(?<=Name=)(.*?)(?=;)').str[0])
        # get unique mutation names not in Pokay
        unmatched_names = pd.Series(names[notinPokay_mask].unique())
        # save unmatched names to file
        if unmatched_names.shape[0] != 0:
            unmatched_names.to_csv(args.names, sep='\t',
                                   index=False, header=False)
            print("")
            print(str(unmatched_names.shape[0]) +
                  " mutation names not matched with functional annotations "
                  "file saved to " + leftover_names_filepath)

    print("")
    print("Processing complete.")
