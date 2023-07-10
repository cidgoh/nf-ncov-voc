#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: madeline

This script converts VCF files that have been annotated into GVF
files, including the functional annotation. Required user
input is a VCF file.
    
The attributes completed by this script are: 
['ID', 'Name', 'chrom_region', 'protein', 'ps_filter', 'ps_exc', 'mat_pep_id',
'mat_pep_desc','mat_pep_acc', 'ro', 'ao', 'dp', 'sample_size', 'Reference_seq',
'Variant_seq', 'nt_name', 'aa_name', 'vcf_gene', 'mutation_type',
'viral_lineage', 'multi_aa_name', 'multiaa_comb_mutation',
'alternate_frequency']
"""

import argparse
import pandas as pd
import numpy as np
import json
from functions import parse_INFO, find_sample_size, \
    unnest_multi, get_unknown_labels, separate_attributes, rejoin_attributes, \
        clade_defining_threshold, split_names, map_pos_to_gene_protein
from functions import empty_attributes, gvf_columns, vcf_columns, pragmas


def vcftogvf(vcf, strain, GENE_PROTEIN_POSITIONS_DICT, names_to_split, sample_size):
    vcf_df = pd.read_csv(vcf, sep='\t', names=vcf_columns)
    # get variant-calling source
    var_cols = get_unknown_labels(vcf_df)
    # remove pragmas
    vcf_df = vcf_df[~vcf_df['#CHROM'].str.contains("#")]
    # restart index from 0
    vcf_df = vcf_df.reset_index(drop=True)
    
    # expand INFO column into multiple columns
    vcf_df = parse_INFO(vcf_df, var_cols)

    # create an empty df to make the new GVF in
    new_gvf = pd.DataFrame(index=range(0, len(vcf_df)), columns=gvf_columns)

    # fill in GVF columns from VCF
    new_gvf['#seqid'] = vcf_df['#CHROM']
    new_gvf['#source'] = '.'
    new_gvf['#start'] = vcf_df['POS']
    # 'end' is not used in creating the COVID-MVP heatmap,
    # but is a required GVF column, so use 'POS' for 'end' as well
    new_gvf['#end'] = vcf_df['POS'] 
    new_gvf['#score'] = '.'
    new_gvf['#strand'] = '+'
    new_gvf['#phase'] = '.'
    if "type" in vcf_df.columns: # "type" is not an attribute of INFO for wastewater
        new_gvf['#type'] = vcf_df['type']
    else:
        new_gvf['#type'] = '.'
    # fill '#attributes' column with empty key-value pairs to fill in later
    new_gvf['#attributes'] = empty_attributes
            
    # expand #attributes into columns to fill in separately
    new_gvf = separate_attributes(new_gvf)

    # fill in attributes from vcf_df columns by name if they exist
    vcf_df_cols_to_add = ['nt_name', 'aa_name', 'vcf_gene', 'mutation_type',
                        'ps_filter', 'ps_exc', 'mat_pep_id','mat_pep_desc',
                        'mat_pep_acc', 'Reference_seq', 'Variant_seq',
                        "dp", "ro", "ao"]
    for column in list(set(vcf_df.columns) & set(vcf_df_cols_to_add)):
        # drop nans if they exist
        vcf_df[column] = vcf_df[column].fillna('')
        new_gvf[column] = vcf_df[column]

    # add other attributes
    new_gvf['sample_size'] = sample_size
    new_gvf['Name'] = vcf_df["Names"]
    new_gvf['viral_lineage'] = strain
    new_gvf['alternate_frequency'] = vcf_df["AF"]
    
    # add chrom_region and protein attributes
    gene_names, protein_names = map_pos_to_gene_protein(
        vcf_df['POS'].astype(int), GENE_PROTEIN_POSITIONS_DICT)
    new_gvf['chrom_region'] = gene_names
    new_gvf['protein'] = protein_names
    
    # add clade_defining attribute
    new_gvf = clade_defining_threshold(args.clades_threshold,
                                             new_gvf, sample_size)
    
    # split up composite mutation names into separate rows
    if names_to_split != 'n/a':
        new_gvf = split_names(names_to_split, new_gvf)
        
    # add 'ID' attribute: here, rows with the same entry in 'Name'
    # get the same ID (should all be different)
    new_gvf['ID'] = 'ID_' + new_gvf.groupby('Name', sort=False).ngroup().astype(str)
    
    # merge attributes back into a single column
    new_gvf = rejoin_attributes(new_gvf, empty_attributes)
    
    return new_gvf



def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a annotated VCF file to a GVF '
                    'file with functional annotation')
    parser.add_argument('--vcffile', type=str, default=None,
                        help='Path to a snpEFF-annotated VCF file')
    parser.add_argument('--size_stats', type=str, default='n/a',
                        help='Statistics file for for size extraction')
    parser.add_argument('--clades_threshold', type=float,
                        default=0.75,
                        help='Alternate frequency cutoff for '
                             'clade-defining mutations')
    parser.add_argument('--gene_positions', type=str,
                        default=None,
                        help='gene positions in JSON format')
    # --names_to_split needs updating: 13 January, 2023
    parser.add_argument('--names_to_split', type=str,
                        default='n/a',
                        help='.tsv of multi-aa mutation names to '
                             'split up into individual aa names')
    parser.add_argument('--strain', type=str,
                        default='n/a',
                        help='Lineage; user mode is if strain="n/a"')
    parser.add_argument("--wastewater", help="Activate wastewater data mode",
                        action="store_true")
    parser.add_argument('--outgvf', type=str,
                        help='Filename for the output GVF file')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    
    # Reading the gene & proetin coordinates of SARS-CoV-2 genome
    with open(args.gene_positions) as fp:
        GENE_PROTEIN_POSITIONS_DICT = json.load(fp)
    # Assigning the vcf file to a variable
    vcf_file = args.vcffile

    # print("Processing: " + vcf_file)

    sample_size = find_sample_size(args.size_stats, args.strain, vcf_file, args.wastewater)
    
    # create gvf from annotated vcf (ignoring pragmas for now)
    gvf = vcftogvf(vcf_file, args.strain, GENE_PROTEIN_POSITIONS_DICT,
                   args.names_to_split, sample_size)
    
    # add pragmas to df, then save to .gvf
    # columns are now 0, 1, ...
    final_gvf = pd.DataFrame(np.vstack([gvf.columns, gvf]))
    final_gvf = pragmas.append(final_gvf)
    filepath = args.outgvf  # outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)

    print("")
    print("Processing complete.")
