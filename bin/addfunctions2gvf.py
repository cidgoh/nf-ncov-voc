#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 08:25:48 2023

@author: madeline

This script annotates GVF files with the functional annotation.

The attributes completed by this script are: 
["measured_variant_functional_effect", "variant_functional_effect_description", "source",
 "citation", "comb_mutation"]

"""

import argparse
import pandas as pd
import numpy as np
from functions import separate_attributes, rejoin_attributes
from functions import empty_attributes, gvf_columns, vcf_columns

# Function to parse command line arguments
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
    parser.add_argument('--functional_annotation_resource', type=str,
                        default='Pokay', help='Functional annotation file identifier; can be versioned (eg. "Pokay_v1.0")')
    parser.add_argument("--names", type=str, default='n/a',
    			help="Save mutation names without "
                              "functional annotations to "
                              "this .txt filename for "
                              "troubleshooting purposes")
    return parser.parse_args()

# Function to add template annotations to GVF file
def add_template_annotations(gvf, annotation_file, annotation_resource):
    
    # expand #attributes into columns to fill in separately
    gvf = separate_attributes(gvf)
    
    # drop columns that are going to be re-added in the merge
    functional_attributes = ["measured_variant_functional_effect", "variant_functional_effect_description", 
                            "URL", "citation", "organism", "reference_accession", "reference_database_name", 'CVX_code',
                            'DrugBank_Accession_Number', 'Antibody_Registry_ID', "author", "publication_year", "DOI", "PMID", "peer_review_status",
                            "curator", "mutation_functional_annotation_resource"] #"comb_mutation"
    gvf = gvf.drop(columns=functional_attributes)

    # load functional annotations spreadsheet
    df = pd.read_csv(annotation_file, sep='\t', header=0)
    # replace spaces in column names with underscores to match the gvf attributes
    template_columns = df.columns.tolist()
    underscore_columns = [col.replace(" ", "_") for col in template_columns]
    rename_dict = dict(zip(template_columns, underscore_columns))
    df = df.rename(columns=rename_dict)
    # if no author, fill with "UNKNOWN"
    df['author'] = df['author'].fillna('UNKNOWN')
    # remove any leading/trailing spaces
    for column in df.columns:
        df[column] = df[column].fillna('')
        df[column] = df[column].astype(str).str.strip()

    # merge annotated vcf and functional annotation files by 'original_mutation_description' and 'protein_symbol'
    df['citation'] = df['author'] + ' et al. (' + df['publication_year'].str.replace(".0", "", regex=False) + ')'
    df_columns = functional_attributes + ["original_mutation_description", "protein_symbol"]
    df = df[df_columns]

    merged_df = pd.merge(gvf, df, on=['original_mutation_description', 'protein_symbol'], how='left') #, 'alias'

    # data cleaning
    merged_df['comb_mutation'] = merged_df['comb_mutation'].str.replace(
        "B.1.617.2\\tT19R", "T19R", regex=False)
    
    # update ID attribute based on mutation groups
    
    # collect all mutation groups (including reference mutation) in
    # merged_df["mutation_group"], sorted alphabetically
    
    # join columns with commas in between
    merged_df["mutation_group"] = \
        merged_df["original_mutation_description"].astype(str) + "," + \
        merged_df["comb_mutation"].astype(str) + "," + \
        merged_df['multiaa_comb_mutation'].astype(str)
    # cleaning: remove nans, quotations marks, spaces
    for x in [' ', 'nan', "'", '"', 'n/a']:
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
    # represented individually in 'original_mutation_description' (True/False)
    unique_original_mutation_description_entries = set(merged_df['original_mutation_description'].tolist())
    group_fully_represented = [set(x).issubset(unique_original_mutation_description_entries)
                               for x in nested_list]
    merged_df['group_fully_represented'] = group_fully_represented
    # drop rows with mutation group members not found in 'original_mutation_description',
    # leaving the index unchanged
    merged_df = merged_df[merged_df['group_fully_represented']==True]

    # update 'ID' attribute: now, rows with the same entry
    # in 'mutation_group' get the same ID
    merged_df['ID'] = 'ID_' + merged_df.groupby(
        'mutation_group', sort=False).ngroup().astype(str)
    
    # change semicolons in function descriptions to colons
    merged_df['variant_functional_effect_description'] = merged_df[
        'variant_functional_effect_description'].str.replace(';', ':')
    
    # add functional_description_resource attribute
    merged_df['functional_annotation_resource'] = annotation_resource

    # make sure 'publication_year' is an integer
    merged_df['publication_year'] = merged_df['publication_year'].str.replace('.0', '', regex=False)

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
    template_annotated_gvf = add_template_annotations(gvf, args.functional_annotations, args.functional_annotation_resource)

    # add pragmas to df, then save to .gvf
    # columns are now 0, 1, ...
    final_gvf = pd.DataFrame(np.vstack([template_annotated_gvf.columns,
                                            template_annotated_gvf]))
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
        notintemplate_mask = template_annotated_gvf["#attributes"].str.contains("measured_variant_functional_effect=;")
        # extract all mutation names from #attributes column
        names = pd.Series(template_annotated_gvf["#attributes"].str.findall('(?<=original_mutation_description=)(.*?)(?=;)').str[0])
        # get unique mutation names not in template
        unmatched_names = pd.Series(names[notintemplate_mask].unique())
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
