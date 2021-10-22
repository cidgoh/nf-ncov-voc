#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:47:11 2021

@author: madeline
"""

'''
This script converts GVF files to TSVs, for later conversion to an HTML case report.
'''


import argparse
import pandas as pd


def parse_args():
    
    parser = argparse.ArgumentParser(
        description='Converts a GVF file to a TSV')
    parser.add_argument('--gvf_directory', type=str, default=None, nargs='*',
                        help='Paths to n GVF files, separated by a space')
    parser.add_argument('--outtsv', type=str, default=None,
                        help='Output filepath for finished .tsv')

    return parser.parse_args()


def gvf2tsv(gvf):
    #read in gvf
    gvf_columns = ['#seqid','#source','#type','#start','#end','#score','#strand','#phase','#attributes']

    df = pd.read_csv(gvf, sep='\t', names=gvf_columns)    
    df = df[~df['#seqid'].str.contains("#")] #remove pragmas and original header
    df = df.reset_index(drop=True) #restart index from 0

    #split #attributes column into separate columns for each tag
    attributes = df['#attributes'].str.split(pat=';').apply(pd.Series) #split at ;, form dataframe
    attributes = attributes.drop(labels=len(attributes.columns)-1, axis=1) #last column is a copy of the index so drop it
    
    for column in attributes.columns:
        split = attributes[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0].lower()
        content = split[1]
        attributes[column] = content #ignore "tag=" in column content
        attributes.rename(columns={column:title}, inplace=True) #make attribute tag as column label

    #replace attributes column in the original df with the new separated out attributes
    df = pd.concat((df, attributes), axis=1)
    df = df.drop(labels='#attributes', axis=1)
    
    #remove '#' from column names
    df.columns = df.columns.str.replace("#", "")
    
    #drop unwanted columns
    df = df.drop(labels=['source', 'seqid', 'type', 'end', 'strand', 'score', 'phase', 'id'], axis=1)

    #rename 'dp' column to 'sequence_depth', make 'viral_lineage' plural
    df = df.rename(columns={'dp': 'sequence_depth', 'viral_lineage': 'viral_lineages'})
    
    return df


if __name__ == '__main__':
    
    args = parse_args()
    
    filepath = args.outtsv
    gvf_files = args.gvf #this is a list of strings

    #convert all gvf files to tsv and concatenate them
    print("Processing:")
    print(gvf_files[0])
    tsv_df = gvf2tsv(gvf_files[0])
    for gvf in gvf_files[1:]:
        print(gvf)
        new_tsv_df = gvf2tsv(gvf)
        tsv_df = pd.concat([tsv_df, new_tsv_df], ignore_index=True)   
    
    
    #find identical rows across strains, and keep only one row.

    #change n/a to 0 in 'ao' for counting purposes
    tsv_df['ao'] = tsv_df['ao'].str.replace("n/a", "0")

    for colname in ['sequence_depth', 'ao', 'ro']:
        #split up at commas into new columns: make a new mini-df
        split_series = tsv_df[colname].str.split(pat=',').apply(pd.Series)
        #rename series columns to 'ao_0', 'a0_1', etc.
        split_series.columns = [colname + '_' + str(name) for name in split_series.columns.values]
        #ensure all counts are numeric
        for column in split_series.columns:
            split_series[column] = pd.to_numeric(split_series[column], errors='coerce')
        #append series to tsv_df
        tsv_df = pd.concat([tsv_df, split_series], axis=1)
        
    cols_to_check = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 'multiaa_comb_mutation', 'start', 'function_category', 'citation', 'comb_mutation', 'function_description', 'heterozygosity']

    agg_dict = dict((col,'first') for col in tsv_df.columns.values.tolist())
    agg_dict['viral_lineages'] = ', '.join
    agg_dict['clade_defining'] = ', '.join
    
    #sum split columns
    for string in ['sequence_depth_', 'ao_', 'ro_']:
        relevant_keys = [key for key, value in agg_dict.items() if string in key.lower()]
        for key in relevant_keys:
            agg_dict[key] = 'sum'
   
    final_df = tsv_df.groupby(cols_to_check).agg(agg_dict)
    
    #rejoin split columns (ao, sequence_depth, ro) with comma separation
    for string in ['sequence_depth_', 'ao_', 'ro_']:
        colnames = [i for i in tsv_df.columns.values.tolist() if string in i]
        final_df[string + 'combined'] = final_df[colnames].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
        #drop split columns
        final_df = final_df.drop(labels=colnames, axis=1)
 
    #replace ao, ro, sequence_depth with the added up columns
    final_df = final_df.drop(labels=['ao', 'ro', 'sequence_depth'], axis=1)
    final_df = final_df.rename(columns={'sequence_depth_combined': 'sequence_depth', 'ro_combined': 'ro', 'ao_combined': 'ao'})

    #reorder columns
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 
       'multiaa_comb_mutation', 'start', 'vcf_gene', 'chrom_region',
       'mutation_type', 'sequence_depth', 'ps_filter', 'ps_exc', 'mat_pep_id',
       'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao', 'reference_seq',
       'variant_seq', 'viral_lineages', 'function_category', 'citation',
       'comb_mutation', 'function_description', 'heterozygosity',
       'clade_defining', 'who_variant', 'status',
       'voi_designation_date', 'voc_designation_date',
       'vum_designation_date']
    final_df = final_df[cols]
    
    #save the final df to a .tsv
    final_df.to_csv(filepath, sep='\t', index=False)
    print("")
    print("Processing complete.")
    print("Saved as: " + filepath)

    