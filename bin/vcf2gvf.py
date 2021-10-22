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
import numpy as np
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a snpEFF-annotated VCF file to a GVF file with functional annotation')
    parser.add_argument('--vcffile', type=str, default=None,
                        help='Path to a snpEFF-annotated VCF file')
    parser.add_argument('--functional_annotations', type=str, default=None,
                        help='TSV file of functional annotations')
    parser.add_argument('--clades', type=str, default=None,
                        help='TSV file of outbreak.info clade-defining mutations')
    parser.add_argument('--clades_threshold', type=float,
                        default=0.75,
                        help='Alternate frequency cutoff for clade-defining mutations')
    parser.add_argument('--gene_positions', type=str,
                        default=None,
                        help='gene positions in JSON format')
    parser.add_argument('--names_to_split', type=str,
                        default=None,
                        help='.tsv of multi-aa mutation names to split up into individual aa names')
    parser.add_argument('--strain', type=str,
                        default='n/a',
                        help='lineage')
    parser.add_argument('--outvcf', type=str,
                        help='Filename for the output GVF file')
    parser.add_argument("--single_genome", help="VCF file is of a single genome", action="store_true")
    parser.add_argument("--names", help="Save mutation names without functional annotations to TSV files for troubleshooting purposes", action="store_true")
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


def clade_defining_threshold(threshold, df):
    """Specifies the clade_defining attribute as True if AF > threshold, False if AF <= threshold, and n/a if the VCF is for a single genome"""
    if args.single_genome:
        df["#attributes"] = df["#attributes"].astype(str)  + "clade_defining=n/a;"
    else:
        df.loc[df.AF > threshold, "#attributes"] = df.loc[df.AF > threshold, "#attributes"].astype(str)  + "clade_defining=True;"
        df.loc[df.AF <= threshold, "#attributes"] = df.loc[df.AF <= threshold, "#attributes"].astype(str)  + "clade_defining=False;"
    return df


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
        title = split[0].drop_duplicates().tolist()[0]
        if isinstance(title, str):
            title = title.lower()
        content = split[1]
        info[column] = content #ignore "tag=" in column content
        info.rename(columns={column:title}, inplace=True) #make attribute tag as column label
    
    #add 'INFO' attributes by name
    for column in ['ps_filter', 'ps_exc', 'mat_pep_id', 'mat_pep_desc', 'mat_pep_acc']:
        info[column] = info[column].fillna('') #drop nans if they exist
        new_df['#attributes'] = new_df['#attributes'].astype(str) + column + '=' + info[column].astype(str) + ';'
       
    #add ro, ao, dp
    unknown = df['unknown'].str.split(pat=':').apply(pd.Series)
    if args.single_genome:
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ro=n/a;ao=n/a;dp=1;'
    else:
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ro=' + unknown[3].astype(str) + ';'
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'ao=' + unknown[5].astype(str) + ';'
        new_df['#attributes'] = new_df['#attributes'].astype(str) + 'dp=' + info['dp'].astype(str) + ';'

    #add af column for clade-defining cutoff (af=ao/dp)
    new_df['AF'] =  unknown[5].astype(int) / info['dp'].astype(int)
        
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
      
    new_df = new_df[gvf_columns + ['multiaa_comb_mutation', 'AF']] #only keep the columns needed for a gvf file, plus multiaa_comb_mutation to add to comb_mutation later
    #new_df.to_csv('new_df.tsv', sep='\t', index=False, header=False)
    return new_df



#takes 4 arguments: the output df of vcftogvf.py, the functional annotation file, the clade defining mutations tsv, the strain name, and the names_to_split tsv.
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
    merged_df["mutation_group"] = merged_df["comb_mutation"].astype(str) + ", '" + merged_df["mutation"].astype(str) + "', " + merged_df['multiaa_comb_mutation'].astype(str)
    merged_df["mutation_group"] = merged_df["mutation_group"].str.replace("nan, ", "")
    merged_df["mutation_group"] = merged_df["mutation_group"].str.rstrip(' ').str.rstrip(',')

    #separate the mutation_group column into its own df with one mutation per column
    mutation_groups = merged_df["mutation_group"].str.split(pat=',').apply(pd.Series)
    mutation_groups = mutation_groups.apply(lambda s:s.str.replace("'", ""))
    mutation_groups = mutation_groups.apply(lambda s:s.str.replace(" ", ""))
    mutation_groups = mutation_groups.transpose() #now each mutation has a column instead
    #sort each column alphabetically
    sorted_df = mutation_groups

    for column in mutation_groups.columns:
        sorted_df[column] = mutation_groups.sort_values(by=column, ignore_index=True)[column]
    sorted_df = sorted_df.transpose()
    
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


    #get clade_defining status, and then info from clades file
    clades = pd.read_csv(clade_file, sep='\t', header=0) #load clade-defining mutations file

    #find the relevant pango_lineage line in the clade file that matches args.strain
    cladefile_strain = 'None'
    available_strains = clades['pango_lineage'].tolist()
    for strain in available_strains:
        for pango_strain in strain.replace("*","").split(','):
            if args.strain.startswith(pango_strain):
                cladefile_strain = strain
    
    #if strain in available_strains:
    if cladefile_strain != 'None':
        clades = clades.loc[clades.pango_lineage == cladefile_strain] #only look at the relevant part of that file
        clades = clades.replace(np.nan, '', regex=True)

        #extract status, WHO strain name, etc. from clades file
        who_variant = clades['who_variant'] 
        who_variant = clades.iloc[0]['who_variant']
        status = clades.iloc[0]['status']
        voi_designation_date = clades.iloc[0]['voi_designation_date']
        voc_designation_date = clades.iloc[0]['voc_designation_date']
        vum_designation_date = clades.iloc[0]['vum_designation_date']
      
        #get True/False/n/a designation for clade-defining status
        merged_df = clade_defining_threshold(args.clades_threshold, merged_df)
        
        #add remaining attributes from clades file
        merged_df["#attributes"] = merged_df["#attributes"].astype(str) + "who_variant=" + who_variant + ';' + "status=" + status + ';' + "voi_designation_date=" + voi_designation_date + ';' + "voc_designation_date=" + voc_designation_date + ';' + "vum_designation_date=" + vum_designation_date + ';'
    else:
        merged_df["#attributes"] = merged_df["#attributes"].astype(str)  + "clade_defining=n/a;" + "who_variant=n/a;" + "status=n/a;" + "voi_designation_date=n/a;" + "voc_designation_date=n/a;" + "vum_designation_date=n/a;"

        
    #add ID to attributes
    merged_df["#attributes"] = 'ID=' + merged_df['id'].astype(str) + ';' + merged_df["#attributes"].astype(str)
    
    if args.names:
        #get list of names in tsv but not in functional annotations, and vice versa, saved as a .tsv
        tsv_names = gvf["mutation"].unique()
        functional_annotation_names = df["mutation"].unique()
        print(str(np.setdiff1d(tsv_names, functional_annotation_names).shape[0]) + "/" + str(tsv_names.shape[0]) + " mutation names were not found in functional_annotations")
        leftover_names = pd.DataFrame({'in_tsv_only':np.setdiff1d(tsv_names, functional_annotation_names)})
        leftover_names["strain"] = strain
        '''
        clade_names = clades["mutation"].unique()
        leftover_clade_names = pd.DataFrame({'unmatched_clade_names':np.setdiff1d(clade_names, tsv_names)})
        leftover_clade_names["strain"] = strain
        '''
        return merged_df[gvf_columns], leftover_names, gvf["mutation"].tolist() #, leftover_clade_names
    
    else:
        return merged_df[gvf_columns]
    
    
      

if __name__ == '__main__':
    
    args = parse_args()
    with open(args.gene_positions) as fp:
        GENE_POSITIONS_DICT = json.load(fp)
    annotation_file = args.functional_annotations
    clade_file = args.clades

    #make empty list in which to store mutation names from all strains in the folder together
    all_strains_mutations = []
    leftover_df = pd.DataFrame() #empty dataframe to hold unmatched names
    #unmatched_clade_names = pd.DataFrame() #empty dataframe to hold unmatched clade-defining mutation names
    pragmas = pd.DataFrame([['##gff-version 3'], ['##gvf-version 1.10'], ['##species NCBI_Taxonomy_URI=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049']]) #pragmas are in column 0

    file = args.vcffile


    print("Processing: " + file)
        
    #create gvf from annotated vcf (ignoring pragmas for now)
    gvf = vcftogvf(file, args.strain, GENE_POSITIONS_DICT, args.names_to_split)
    #add functional annotations
    if args.names:
        '''
        annotated_gvf, leftover_names, mutations, \
        leftover_clade_names = add_functions(gvf,
                                             annotation_file,
                                             clade_file, args.strain)
        '''
        annotated_gvf, leftover_names, mutations = add_functions(gvf,
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
        #unmatched_clade_names = unmatched_clade_names.append(leftover_clade_names)
        #save unmatched names (in tsv but not in functional_annotations) across all strains to a .tsv file
        leftover_names_filepath = "leftover_names.tsv"
        leftover_df.to_csv(leftover_names_filepath, sep='\t', index=False)
        print("")
        print("Mutation names not found in functional annotations file saved to " + leftover_names_filepath)
        '''
        #save unmatched clade-defining mutation names to a .tsv file
        leftover_clade_names_filepath = "leftover_clade_defining_names.tsv"
        unmatched_clade_names.to_csv(leftover_clade_names_filepath, sep='\t', index=False)
        print("Clade-defining mutation names not found in the annotated VCF saved to " + leftover_clade_names_filepath)
        '''
        #print number of unique mutations across all strains    
        flattened = [val for sublist in all_strains_mutations for val in sublist]
        arr = np.array(flattened)
        print("# unique mutations in VCF file: ", np.unique(arr).shape[0])

    print("")        
    print("Processing complete.")
        