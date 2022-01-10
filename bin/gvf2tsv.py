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
        num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(lineage.replace("*",""))]['num_seqs']
        num_seqs = num_seqs.str.replace(",","").values.tolist()
        variant_num_seqs = variant_num_seqs + num_seqs
    variant_pop_size = sum(list(map(int, variant_num_seqs)))

    return variant_pop_size


def add_ao_by_variant_seq(ao_str, variant_seq_str):
    '''
    Takes a pair of strings like ao_str="7,7,24" and variant_seq_str="T,T,T" or "8,21,3" and "T,T,A".
    Output should be a string "T=38" (first case) or "T=29, A=3" (second case).
    '''
    ao_list = ao_str.split(',')
    ao_list = [int(x) for x in ao_list]
    var_str_list = variant_seq_str.split(',')
    zipped_lists = list(zip(var_str_list, ao_list))
    
    ao_dict = dict()
    for pair in zipped_lists:
        #if a variant seq isn't in the dictionary, add it
        if pair[0] not in ao_dict:
            ao_dict[pair[0]] = pair[1]
        else:
            ao_dict[pair[0]] += pair[1]
    
    joined_string = ''
    for item in list(ao_dict.items()):
        joined_string = joined_string + item[0] + '=' + str(item[1]) + ', '
    joined_string = joined_string.rstrip(", ")
    
    return joined_string
    

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
    
    #make ro, dp, and obs_sample_size numeric
    for colname in ['ro', 'dp', 'obs_sample_size']:
        tsv_df[colname] = pd.to_numeric(tsv_df[colname], errors='coerce')
    
    cols_to_check = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 'multiaa_comb_mutation', 'start', 'function_category', 'citation', 'comb_mutation', 'function_description', 'heterozygosity']

    agg_dict = dict((col,'first') for col in tsv_df.columns.values.tolist())

    #join some columns with commas
    agg_dict['viral_lineages'] = ', '.join
    agg_dict['clade_defining'] = ','.join
    agg_dict['ao'] = ','.join
    agg_dict['variant_seq'] = ','.join
            
    #sum dp, ro and sample_size
    for key in ['dp', 'ro', 'obs_sample_size']:
        agg_dict[key] = 'sum'
   
    final_df = tsv_df.groupby(cols_to_check).agg(agg_dict)
   

    #add ao according to the heterogeneous mutations
    final_df['ao_by_var_seq'] = [add_ao_by_variant_seq(x, y) for x, y in zip(final_df['ao'], final_df['variant_seq'])]
    
    #remove 'who_variant'; rename 'multiaa_comb_mutation'
    final_df = final_df.drop(labels=['who_variant'], axis=1)
    final_df = final_df.rename(columns={'multiaa_comb_mutation': 'multiaa_mutation_split_names'})
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
    final_df.clade_defining_status = final_df.clade_defining_status.str.replace('nan=nan; ', '')
    #strip trailing spaces and semicolons
    final_df.clade_defining_status = final_df.clade_defining_status.str.rstrip("; ")

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
    #return an ordered list of lineages
    for row in final_df['viral_lineages']:
        if ' ' in row:
            lineage_list = row.split(', ')
            mylist = list(set(lineage_list))
            '''
            #order lineages alphanumerically
            #split each element of mylist into a sublist, split at the first '.'
            split_list = 
            '''
            row_str = ', '.join(str(e) for e in mylist)
            mask = final_df['viral_lineages']==row
            final_df.loc[mask, 'viral_lineages'] = row_str  
    
    #reorder columns
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name', 
       'multiaa_mutation_split_names', 'start', 'vcf_gene', 'chrom_region',
       'mutation_type', 'dp', 'obs_sample_size', 'variant_pop_size', 'ps_filter', 'ps_exc', 'mat_pep_id',
       'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao', 'ao_by_var_seq', 'reference_seq',
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

    