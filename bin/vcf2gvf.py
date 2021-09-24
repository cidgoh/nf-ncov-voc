#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 11:06:22 2021

@author: madeline
"""

'''
This script converts VCF files that have been annotated by snpEFF into GVF files, including the functional annotation.
Required user input is a VCF file.

test case: 
--vcffile /home/madeline/Downloads/B.1.525.variants.filtered.annotated.filtered.vcf --strain B.1.525 --outvcf b1525out.gvf
'''

import argparse
import pandas as pd
import re
import glob
import os
import numpy as np
import json
from natsort import index_natsorted

def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a snpEFF-annotated VCF file to a GVF file with functional annotation')
    parser.add_argument('--vcffile', type=str, default=None,
                        help='Path to a snpEFF-annotated VCF file')
    parser.add_argument('--pokay', type=str, default='../.github/data/functional_annotation/functional_annotation_V.0.3.tsv',
                        help='Anoosha\'s v.3 parsed pokay .tsv file')
    parser.add_argument('--clades', type=str, default='../.github/data/clade_defining/clade_defining_mutations.tsv',
                        help='.tsv of clade-defining mutations')
    parser.add_argument('--gene_positions', type=str,
                        default='../.github/data/gene_coordinates/gene_positions.json',
                        help='gene positions in json format')
    parser.add_argument('--names_to_split', type=str,
                        default='../.github/data/multi_aa_names/mutation_names_to_split.tsv',
                        help='.tsv of multi-aa mutation names to split up into individual aa names')
    parser.add_argument('--strain', type=str,
                        default=None,
                        help='lineage')
    parser.add_argument('--outvcf', type=str,
                        help='Filename for the output .gvf file')
    parser.add_argument("--single_genome", help="VCF file is of a single genome", action="store_true")
    parser.add_argument("--names", help="Save unmatched mutation names to .tsvs for troubleshooting naming formats", action="store_true")
    return parser.parse_args()


def map_pos_to_gene(pos, GENE_POSITIONS_DICT):
    """This function is inspired/lifted from Ivan's code.
    Map a series of nucleotide positions to SARS-CoV-2 genes.
    See https://www.ncbi.nlm.nih.gov/nuccore/MN908947.
    :param pos: Nucleotide position pandas series
    :type pos: int
    :return: series containing SARS-CoV-2 chromosome region names at each nucleotide position in ``pos``
    """
    gene_names = pos.astype(str) #make a series of the same size as pos to put gene names in
    for gene in GENE_POSITIONS_DICT:
        start = GENE_POSITIONS_DICT[gene]["start"]
        end = GENE_POSITIONS_DICT[gene]["end"]
        gene_mask = pos.between(start, end, inclusive=True)
        if gene == "Stem-loop":
            gene_names[gene_mask] = gene + ",3\' UTR"
        else:
            gene_names[gene_mask] = gene
    gene_names[gene_names.str.isnumeric()] = "intergenic"
    return gene_names


gvf_columns = ['#seqid','#source','#type','#start','#end','#score','#strand','#phase','#attributes']
vcf_colnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown']

def vcftogvf(var_data, strain, GENE_POSITIONS_DICT, names_to_split):
     
    df = pd.read_csv(var_data, sep='\t', names=vcf_colnames)    
    df = df[~df['#CHROM'].str.contains("#")] #remove pragmas
    df = df.reset_index(drop=True) #restart index from 0

    new_df = pd.DataFrame(index=range(0,len(df)),columns=gvf_columns)

    #fill in first 8 GVF columns
    new_df['#seqid'] = df['#CHROM']
    new_df['#source'] = '.'
    new_df['#type'] = '.' #info['type']
    new_df['#start'] = df['POS']
    new_df['#end'] = (df['POS'].astype(int) + df['ALT'].str.len() - 1).astype(str)  #this needs fixing
    new_df['#score'] = '.'
    new_df['#strand'] = '+'
    new_df['#phase'] = '.'

    #parse INFO column
    
    #sort out problematic sites tag formats
    df['INFO'] = df['INFO'].str.replace('ps_filter;','ps_filter=;')
    df['INFO'] = df['INFO'].str.replace('ps_exc;','ps_exc=;')
    df['INFO'] = df['INFO'].str.replace('=n/a','')
    
    #parse EFF entry in INFO    
    eff_info = df['INFO'].str.findall('\((.*?)\)') #series: extract everything between parentheses as elements of a list
    eff_info = eff_info.apply(pd.Series)[0] #take first element of list
    eff_info = eff_info.str.split(pat='|').apply(pd.Series) #split at pipe, form dataframe
    #hgvs names
    hgvs = eff_info[3].str.rsplit(pat='c.').apply(pd.Series)
    hgvs_protein = hgvs[0].str[:-1]
    hgvs_nucleotide = 'c.' + hgvs[1]
    
    new_df['nt_name'] = hgvs_nucleotide
    new_df['aa_name'] = hgvs_protein
    
    #change nucleotide names of the form "c.C*4378A" to c.C4378AN; change vcf_gene to "intergenic" here
    asterisk_mask = hgvs_nucleotide.str.contains('\*')
    hgvs_nucleotide[asterisk_mask] = 'c.' + df['REF'] + df['POS'] + df['ALT']
    eff_info[5][asterisk_mask] = "intergenic"
    
    #use nucleotide name where protein name doesn't exist (for 'Name' attribute)
    Names = hgvs[0].str[:-1]
    Names[~Names.str.contains("p.")] =  hgvs_nucleotide #fill in empty protein name spaces with nucleotide names ("c."...)

    new_df["Names"] = Names
    
    new_df['vcf_gene'] = eff_info[5]
    new_df['mutation_type'] = eff_info[1]
    
    new_df["multi_name"] = ''
    new_df["multiaa_comb_mutation"] = ''
    
    new_df['#attributes'] = 'chrom_region=' + map_pos_to_gene(df['POS'].astype(int), GENE_POSITIONS_DICT) + ';' #gene names including IGRs/UTRs

    #make 'INFO' column easier to extract attributes from
    info = df['INFO'].str.split(pat=';').apply(pd.Series) #split at ;, form dataframe
    for column in info.columns:
        split = info[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0].lower()
        content = split[1]
        info[column] = content #ignore "tag=" in column content
        info.rename(columns={column:title}, inplace=True) #make attribute tag as column label
    
    #add 'INFO' attributes by name
    for column in ['ps_filter', 'ps_exc', 'mat_pep_id', 'mat_pep_desc', 'mat_pep_acc']:
        info[column] = info[column].fillna('') #drop nans if they exist
        new_df['#attributes'] = new_df['#attributes'].astype(str) + column + '=' + info[column].astype(str) + ';'
       
    #add ro, ao, dp
    unknown = df['unknown'].str.split(pat=':').apply(pd.Series)
    new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ro=' + unknown[3].astype(str) + ';'
    if args.single_genome:
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ao=n/a;dp=1;'
    else:
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ao=' + unknown[5].astype(str) + ';'
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'dp=' + info['dp'].astype(str) + ';'

    #add columns copied straight from Zohaib's file
    for column in ['REF','ALT']:
        key = column.lower()
        if key=='ref':
            key = 'Reference_seq'
        elif key=='alt':
            key = 'Variant_seq'
        new_df['#attributes'] = new_df['#attributes'].astype(str) + key + '=' + df[column].astype(str) + ';'


    #split multi-aa names from the vcf into single-aa names (multi-row)
    multiaanames = pd.read_csv(names_to_split, sep='\t', header=0) #load names_to_split spreadsheet
    names_to_separate = np.intersect1d(multiaanames, new_df['Names'].str[2:]) #multi-aa names that are in the gvf (list form)
    #now split the rows apart in-place
    for multname in names_to_separate:
        splits = multiaanames[multiaanames['name']==multname]['split_into'].copy() #relevant part of 'split_into' column
        splits_list_1 = splits.str.split(pat=',').tolist()[0] #list of names to split into
        splits_list = [s.replace("'", '') for s in splits_list_1]
        split_index = new_df.index.get_loc(new_df.index[new_df['Names'].str[2:]==multname][0])  #new_df index containing multi-aa name
        seprows = new_df.loc[[split_index]].copy() #copy of rows to alter
        new_df = new_df.drop(split_index) #delete original combined mutation rows
    
        i=0
        for sepname in splits_list:
            seprows['Names'] = "p." + sepname #single-aa name
            seprows["multi_name"] = multname #original multi-aa name
            seprows["multiaa_comb_mutation"] = splits.tolist()[0].replace("'" + sepname + "'",'').replace(",,",',').replace("','","', '").strip(',') #other single-aa names corresponding to this multi-aa mutation
            new_df = pd.concat([new_df.loc[:split_index + i], seprows, new_df.loc[split_index + i:]]).reset_index(drop=True)
            i += 1

    #add attributes
    new_df['#attributes'] = 'Name=' + new_df["Names"] + ';' + new_df['#attributes'].astype(str)
    new_df['#attributes'] = new_df['#attributes'].astype(str) + 'nt_name=' + new_df['nt_name'] + ';'
    new_df['#attributes'] = new_df['#attributes'].astype(str) + 'aa_name=' + new_df['aa_name'] + ';'
    new_df['#attributes'] = new_df['#attributes'].astype(str) + 'vcf_gene=' + new_df['vcf_gene'] + ';' #gene names
    new_df['#attributes'] = new_df['#attributes'].astype(str) + 'mutation_type=' + new_df['mutation_type'] + ';' #mutation type 

    #add strain name, multi-aa notes
    new_df['#attributes'] = new_df['#attributes'] + 'viral_lineage=' + strain + ';'
    new_df['#attributes'] = new_df['#attributes'] + "multi_aa_name=" + new_df["multi_name"] + ';'
    new_df['#attributes'] = new_df['#attributes'] + "multiaa_comb_mutation=" + new_df["multiaa_comb_mutation"] + ';'    
      
    new_df = new_df[gvf_columns + ['multiaa_comb_mutation']] #only keep the columns needed for a gvf file, plus multiaa_comb_mutation to add to comb_mutation later
    #new_df.to_csv('new_df.tsv', sep='\t', index=False, header=False)
    return new_df



#takes 4 arguments: the output df of vcftogvf.py, Anoosha's annotation file from Pokay, the clade defining mutations tsv, the strain name, and the names_to_split tsv.
def add_functions(gvf, annotation_file, clade_file, strain):
    attributes = gvf["#attributes"].str.split(pat=';').apply(pd.Series)
    hgvs_protein = attributes[0].str.split(pat='=').apply(pd.Series)[1] #remember this includes nucleotide names where there are no protein names
    hgvs_nucleotide = attributes[1].str.split(pat='=').apply(pd.Series)[1]
    gvf["mutation"] = hgvs_protein.str[2:] #drop the prefix
    
    #merge annotated vcf and functional annotation files by 'mutation' column in the gvf
    df = pd.read_csv(annotation_file, sep='\t', header=0) #load functional annotations spreadsheet
   
    for column in df.columns:
        df[column] = df[column].str.lstrip()
    merged_df = pd.merge(df, gvf, on=['mutation'], how='right') #add functional annotations
    
    
    #collect all mutation groups (including reference mutation) in a column, sorted alphabetically
    #this is more roundabout than it needs to be; streamline with grouby() later
    merged_df["mutation_group"] = merged_df["comb_mutation"].astype(str) + ", '" + merged_df["mutation"].astype(str) + "'"
    merged_df["mutation_group_combadded"] = merged_df["comb_mutation"].astype(str) + ", '" + merged_df["mutation"].astype(str) + "'," + merged_df['multiaa_comb_mutation'].astype(str)
    mutation_groups = merged_df["mutation_group"].str.split(pat=',').apply(pd.Series)
    mutation_groups = mutation_groups.apply(lambda s:s.str.replace("'", ""))
    mutation_groups = mutation_groups.apply(lambda s:s.str.replace(" ", ""))
    mutation_groups = mutation_groups.transpose() 
    sorted_df = mutation_groups
    for column in mutation_groups.columns:
        sorted_df[column] = mutation_groups.sort_values(by=column, ignore_index=True)[column]
    sorted_df = sorted_df.transpose()
    
    #splitnames_df = merged_df[['mutation_group', 'multiaa_comb_mutation', 'comb_mutation','mutation_group_combadded']] #only keep the columns needed
    #splitnames_df.to_csv('splitnames_df.tsv', sep='\t', index=False, header=True)

    
    #since they're sorted, put everything back into a single cell, don't care about dropna
    df3 = sorted_df.apply(lambda x :','.join(x.astype(str)),axis=1)
    unique_groups = df3.drop_duplicates() 
    unique_groups_multicol = sorted_df.drop_duplicates() 
    merged_df["mutation_group_labeller"] = df3 #for sanity checking
    
    #make a unique id for mutation groups that have all members represented in the vcf
    #for groups with missing members, delete those functional annotations
    merged_df["id"] = 'NaN'
    id_num = 0
    for row in range(unique_groups.shape[0]):
        group_mutation_set = set(unique_groups_multicol.iloc[row])
        group_mutation_set = {x for x in group_mutation_set if (x==x and x!='nan')} #remove nan and 'nan' from set
        gvf_all_mutations = set(gvf['mutation'].unique())
        indices = merged_df[merged_df.mutation_group_labeller == unique_groups.iloc[row]].index.tolist()
        if group_mutation_set.issubset(gvf_all_mutations): #if all mutations in the group are in the vcf file, include those rows and give them an id
            merged_df.loc[merged_df.mutation_group_labeller == unique_groups.iloc[row], "id"] = "ID_" + str(id_num)
            id_num += 1
        else:
            merged_df = merged_df.drop(indices) #if not, drop group rows, leaving the remaining indices unchanged

    #change semicolons in function descriptions to colons
    merged_df['function_description'] = merged_df['function_description'].str.replace(';',':')
    #change heteozygosity column to True/False
    merged_df['heterozygosity'] = merged_df['heterozygosity']=='heterozygous'
    merged_df['citation'] = merged_df['citation'].str.strip() #remove trailing spaces from citation 
    #add key-value pairs to attributes column
    for column in ['function_category', 'source', 'citation', 'comb_mutation', 'function_description', 'heterozygosity']:
        key = column.lower()
        merged_df[column] = merged_df[column].fillna('') #replace NaNs with empty string
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + key + '=' + merged_df[column].astype(str) + ';'


    #if strain is in clades file, merge that too
    clades = pd.read_csv(clade_file, sep='\t', header=0) #load clade-defining mutations file
    available_strains = clades['strain'].drop_duplicates().tolist() 
    if strain in available_strains:
        clades = clades.loc[clades.strain == strain] #only look at the relevant part of that file
        clades = clades.replace(np.nan, '', regex=True)
        #extract status and WHO strain name from clades file
        variant = clades['variant'].drop_duplicates().tolist()[0]
        who_label = clades['who_name'].drop_duplicates().tolist()[0]
        variant_status = clades['variant_status'].drop_duplicates().tolist()[0]
        voi_designation_date = clades['voi_designation_date'].drop_duplicates().tolist()[0]
        voc_designation_date = clades['voc_designation_date'].drop_duplicates().tolist()[0]
        alert_designation_date = clades['alert_designation_date'].drop_duplicates().tolist()[0]
        #merge clades with function-annotated dataframe
        merged_df = pd.merge(clades, merged_df, on=['mutation'], how='right') #add clade-defining mutations
        #change clade-defining attribute to True/False depending on content of 'strain' column
        merged_df.loc[merged_df.strain == strain, "#attributes"] = merged_df.loc[merged_df.strain == strain, "#attributes"].astype(str)  + "clade_defining=True;"
        merged_df.loc[merged_df.strain != strain, "#attributes"] = merged_df.loc[merged_df.strain != strain, "#attributes"].astype(str)  + "clade_defining=False;"
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "who_label=" + who_label + ';'
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "variant=" + variant + ';'
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "variant_status=" + variant_status + ';'
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "voi_designation_date=" + voi_designation_date + ';'
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "voc_designation_date=" + voc_designation_date + ';'
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "alert_designation_date=" + alert_designation_date + ';'
        
        
    #add ID to attributes
    merged_df["#attributes"] = 'ID=' + merged_df['id'].astype(str) + ';' + merged_df["#attributes"].astype(str)
    
    if args.names:
        #get list of names in tsv but not in functional annotations, and vice versa, saved as a .tsv
        tsv_names = gvf["mutation"].unique()
        pokay_names = df["mutation"].unique()
        print(str(np.setdiff1d(tsv_names, pokay_names).shape[0]) + "/" + str(tsv_names.shape[0]) + " mutation names were not found in pokay")
        in_pokay_only = pd.DataFrame({'in_pokay_only':np.setdiff1d(pokay_names, tsv_names)})
        in_tsv_only = pd.DataFrame({'in_tsv_only':np.setdiff1d(tsv_names, pokay_names)})
        leftover_names = in_tsv_only
        leftover_names["strain"] = strain
    
        clade_names = clades["mutation"].unique()
        leftover_clade_names = pd.DataFrame({'unmatched_clade_names':np.setdiff1d(clade_names, tsv_names)})
        leftover_clade_names["strain"] = strain
    
        return merged_df[gvf_columns], leftover_names, gvf["mutation"].tolist(), leftover_clade_names
    
    else:
        return merged_df[gvf_columns]
    
    
      

if __name__ == '__main__':
    
    args = parse_args()
    with open(args.gene_positions) as fp:
        GENE_POSITIONS_DICT = json.load(fp)
    annotation_file = args.pokay
    clade_file = args.clades

    #make empty list in which to store mutation names from all strains in the folder together
    all_strains_mutations = []
    leftover_df = pd.DataFrame() #empty dataframe to hold unmatched names
    unmatched_clade_names = pd.DataFrame() #empty dataframe to hold unmatched clade-defining mutation names
    pragmas = pd.DataFrame([['##gff-version 3'], ['##gvf-version 1.10'], ['##species NCBI_Taxonomy_URI=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049']]) #pragmas are in column 0

    file = args.vcffile

    print("Processing: " + file)
        
    #create gvf from annotated vcf (ignoring pragmas for now)
    gvf = vcftogvf(file, args.strain, GENE_POSITIONS_DICT, args.names_to_split)
    #add functional annotations
    if args.names:
        annotated_gvf, leftover_names, mutations, \
        leftover_clade_names = add_functions(gvf,
                                             annotation_file,
                                             clade_file, args.strain)
    else:
        annotated_gvf = add_functions(gvf, annotation_file,
                                      clade_file, args.strain)
    #add pragmas to df, then save to .gvf
    annotated_gvf = pd.DataFrame(np.vstack([annotated_gvf.columns, annotated_gvf])) #columns are now 0, 1, ...
    final_gvf = pragmas.append(annotated_gvf)
    filepath = args.outvcf #outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)
    
    
    #get name troubleshooting reports
    if args.names:        
        all_strains_mutations.append(mutations)
        leftover_df = leftover_df.append(leftover_names)
        unmatched_clade_names = unmatched_clade_names.append(leftover_clade_names)
        
        #save unmatched names (in tsv but not in Pokay) across all strains to a .tsv file
        leftover_names_filepath = "_leftover_names.tsv"
        leftover_df.to_csv(leftover_names_filepath, sep='\t', index=False)
        print("")
        print("Mutation names not found in Pokay saved to " + leftover_names_filepath)
    
        #save unmatched clade-defining mutation names across all strains to a .tsv file
        leftover_clade_names_filepath = "leftover_clade_defining_names.tsv"
        unmatched_clade_names.to_csv(leftover_clade_names_filepath, sep='\t', index=False)
        print("Clade-defining mutation names not found in the annotated VCFs saved to " + leftover_clade_names_filepath)

        #print number of unique mutations across all strains    
        flattened = [val for sublist in all_strains_mutations for val in sublist]
        arr = np.array(flattened)
        if args.vcffile:
            print("# unique mutations in " + strain + " VCF file: ", np.unique(arr).shape[0])
        if args.vcfdir:
            print("# unique mutations across all processed VCFs: ", np.unique(arr).shape[0])

    print("")        
    print("Processing complete.")
        