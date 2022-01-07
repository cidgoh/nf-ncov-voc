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
import os


def parse_args():
    
    parser = argparse.ArgumentParser(
        description='Converts GVF files to a TSV report')
    parser.add_argument('--gvf_files', type=str, default=None, nargs='*',
                        help='Paths to GVF files to process')
    parser.add_argument('--clades', type=str, default=None,
                        help='TSV file of WHO strain names and VOC/VOI status')
    parser.add_argument('--outtsv', type=str, default="surveillance_report",
                        help='Filepath for finished .tsv')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--specify_variants', type=str, default=None, nargs='*',
                        help='Name(s) of WHO variant to make report for.  Not case-sensitive.')
    group.add_argument('--all_variants', action="store_true",
                        help='Create reports for all variants, using all available reference lineage gvf files.  Not case-sensitive.')
    parser.add_argument('--table', type=str, default=None,
                        help='Multi-strain TSV file generated in workflow that contains num_seqs column')

    return parser.parse_args()



def match_gvfs_to_who_variant(pango_lineage_list, gvf_files_list):
    
    matched_files = []
    
    for pango_lineage in pango_lineage_list:
            for f in gvf_files_list:
                if f.startswith(pango_lineage.replace("*","")) and f.endswith(".gvf"):
                    matched_files.append(f)

    return matched_files
  
    

def find_variant_pop_size(table, pango_lineage_list):
    strain_tsv_df = pd.read_csv(table, header=0, delim_whitespace=True, usecols=['file', 'num_seqs'])  
    variant_num_seqs = []
    for lineage in pango_lineage_list:
        num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(lineage.replace("*",""))]['num_seqs'].values.tolist()
        num_seqs = [x.split(",") for x in num_seqs]
        num_seqs = [item for sublist in num_seqs for item in sublist]
        variant_num_seqs = variant_num_seqs + num_seqs
    variant_pop_size = sum(list(map(int, variant_num_seqs)))

    return variant_pop_size
    

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
    df = df.rename(columns={'sample_size':'obs_sample_size', 'viral_lineage': 'viral_lineages'})
        
    return df



def streamline_tsv(tsv_df):
    #find identical rows across strains, and keep only one row.

    #change n/a to 0 in 'ao' for counting purposes
    tsv_df['ao'] = tsv_df['ao'].str.replace("n/a", "0")

    for colname in ['dp', 'ao', 'ro']:
        #split up at commas into new columns: make a new mini-df
        split_series = tsv_df[colname].str.split(pat=',').apply(pd.Series)
        #rename series columns to 'ao_0', 'ao_1', etc.
        split_series.columns = [colname + '_' + str(name) for name in split_series.columns.values]
        #ensure all counts are numeric
        for column in split_series.columns:
            split_series[column] = pd.to_numeric(split_series[column], errors='coerce')
        #append series to tsv_df
        tsv_df = pd.concat([tsv_df, split_series], axis=1)
    
    #make sample size numeric
    tsv_df['obs_sample_size'] = pd.to_numeric(tsv_df['obs_sample_size'], errors='coerce')
    
    cols_to_check = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 'multiaa_comb_mutation', 'start', 'function_category', 'citation', 'comb_mutation', 'function_description', 'heterozygosity']

    agg_dict = dict((col,'first') for col in tsv_df.columns.values.tolist())
    agg_dict['viral_lineages'] = ', '.join
    agg_dict['clade_defining'] = ', '.join
    
    #sum split columns and sample_size
    for string in ['dp_', 'ao_', 'ro_', 'obs_sample_size']:
        relevant_keys = [key for key, value in agg_dict.items() if string in key.lower()]
        for key in relevant_keys:
            agg_dict[key] = 'sum'
   
    final_df = tsv_df.groupby(cols_to_check).agg(agg_dict)
    
    #rejoin split columns (ao, sequence_depth, ro) with comma separation
    for string in ['dp_', 'ao_', 'ro_']:
        colnames = [i for i in tsv_df.columns.values.tolist() if string in i]
        final_df[string + 'combined'] = final_df[colnames].apply(lambda row: ','.join(row.values.astype(str)), axis=1)
        #drop split columns
        final_df = final_df.drop(labels=colnames, axis=1)
 
    #replace ao, ro, sequence_depth with the added up columns; remove 'who_variant'
    final_df = final_df.drop(labels=['ao', 'ro', 'dp', 'who_variant'], axis=1)
    final_df = final_df.rename(columns={'dp_combined': 'dp', 'ro_combined': 'ro', 'ao_combined': 'ao', 'multiaa_comb_mutation': 'multiaa_mutation_split_names'})
    #remove trailing zeros and commas from 'ao'
    final_df.ao = final_df.ao.str.replace(',0.0','', regex=True)
    #make 'ao' into comma-separated integers
    final_df.ao = final_df.ao.str.replace('.0,',',', regex=True) 
    #remove any last trailing '.0'
    final_df.ao = final_df.ao.str.split('.').str[0]
    
    #add variant_pop_size
    final_df['variant_pop_size'] = variant_pop_size


    #combine viral_lineages and clade_defining into key-value pairs

    #split viral_lineages and clade_defining by ','
    split_lineages = final_df['viral_lineages'].str.split(pat=',').apply(pd.Series) #split at ,, form dataframe
    split_clade_defining = final_df['clade_defining'].str.split(pat=',').apply(pd.Series) #split at ,, form dataframe
    #go through and make key-value pairs of corresponding columns from each
    final_df['clade_defining_status'] = ''
    for col in split_clade_defining.columns:
        final_df['clade_defining_status'] = final_df['clade_defining_status'] + split_lineages[col].astype(str) + '=' + split_clade_defining[col].astype(str) + '; '
    #drop clade_defining status for n/a strains and empty nan=nan pairs
    final_df.clade_defining_status = final_df.clade_defining_status.str.replace('n/a=n/a; ', 'n/a; ')
    final_df.clade_defining_status = final_df.clade_defining_status.str.replace('= n/a', '=n/a')
    final_df.clade_defining_status = final_df.clade_defining_status.str.replace('nan=nan; ', '')
    #strip trailing spaces and semicolons
    final_df.clade_defining_status = final_df.clade_defining_status.str.rstrip("; ")
    #remove spaces after '='
    final_df.clade_defining_status = final_df.clade_defining_status.str.replace("= ", "=")
    #drop repeated key-value pairs in each row (find these rows as they contain spaces)
    for row in final_df['clade_defining_status']:
        if ' ' in row:
            mylist = row.split('; ')
            newlist = []
            for pair in mylist:
                pair = pair.replace(';', '')
                pair = pair.lstrip(' ')
                newlist.append(pair)
            mylist = list(set(newlist))
            row_str = ', '.join(str(e) for e in mylist)
            mask = final_df['clade_defining_status']==row
            final_df.loc[mask, 'clade_defining_status'] = row_str
            
    #drop repeated lineage names in each row of viral_lineages
    for row in final_df['viral_lineages']:
        if ' ' in row:
            lineage_list = row.split(', ')
            mylist = list(set(lineage_list))
            row_str = ', '.join(str(e) for e in mylist)
            mask = final_df['viral_lineages']==row
            final_df.loc[mask, 'viral_lineages'] = row_str

    
    
    #reorder columns
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 
       'multiaa_mutation_split_names', 'start', 'vcf_gene', 'chrom_region',
       'mutation_type', 'dp', 'obs_sample_size', 'variant_pop_size', 'ps_filter', 'ps_exc', 'mat_pep_id',
       'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao', 'reference_seq',
       'variant_seq', 'function_category', 'citation',
       'comb_mutation', 'function_description', 'heterozygosity',
       'viral_lineages', 'clade_defining_status', 'status',
       'voi_designation_date', 'voc_designation_date',
       'vum_designation_date']
    final_df = final_df[cols]
    
    return final_df




if __name__ == '__main__':
    
    args = parse_args()
    
    filepath = args.outtsv
    clade_file = args.clades
    gvf_files_list = args.gvf_files
    
    #read in WHO variant/PANGO lineage .tsv
    clades = pd.read_csv(clade_file, sep='\t', header=0, usecols=['who_variant', 'pango_lineage']) #load clade-defining mutations file
    
    #get lowercase WHO variant names    
    if args.specify_variants:
        who_variants_list = args.specify_variants 
    elif args.all_variants:
        who_variants_list = clades['who_variant'].tolist()

    #for each variant, create a surveillance report
    for who_variant in who_variants_list:
        who_variant = who_variant.capitalize()
        #get list of relevant pango lineages
        pango_lineages = clades[clades['who_variant']==who_variant]['pango_lineage'].values[0].split(',')  #list of pango lineages from that variant
    
        #get variant population size
        variant_pop_size = find_variant_pop_size(args.table, pango_lineages)
        
        #get list of gvf files pertaining to variant
        gvf_files = match_gvfs_to_who_variant(pango_lineages, gvf_files_list)
        print(str(len(gvf_files)) + " GVF files found for " + who_variant + " variant.")
    
        #if any GVF files are found, create a surveillance report
        if len(gvf_files) > 0:
    
            #convert all gvf files to tsv and concatenate them
            print("Processing:")
            print(gvf_files[0])
            tsv_df = gvf2tsv(gvf_files[0])
            for gvf in gvf_files[1:]:
                print(gvf)
                new_tsv_df = gvf2tsv(gvf)
                tsv_df = pd.concat([tsv_df, new_tsv_df], ignore_index=True)   
            
            #streamline final concatenated df, reorder/rename columns where needed
            final_df = streamline_tsv(tsv_df)
            
            #save report as a .tsv
            filename = filepath + '_' + who_variant + '.tsv'
            final_df.to_csv(filename, sep='\t', index=False)
            print("Processing complete.")
            print(who_variant + " surveillance report saved as: " + filename)
            print("")

    