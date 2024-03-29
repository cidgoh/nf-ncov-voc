#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: madeline

This script converts VCF files that have been annotated into GVF
files. Required user input is a VCF file.
    
The attributes completed by this script are: 
['ID', 'Name', 'gene', 'protein_name', 'protein_symbol', 'protein_id', 'ps_filter', 'ps_exc', 'mat_pep',
'mat_pep_desc','mat_pep_acc', 'ro', 'ao', 'dp', 'sample_size', 'Reference_seq',
'Variant_seq', 'nt_name', 'aa_name', 'hgvs_nt', 'hgvs_aa', 'hgvs_alias', 'vcf_gene', 'mutation_type',
'viral_lineage', 'alternate_frequency']
"""

import argparse
import pandas as pd
import numpy as np
import json
from functions import parse_INFO, find_sample_size, \
    unnest_multi, get_unknown_labels, separate_attributes, rejoin_attributes, \
    clade_defining_threshold, map_pos_to_gene_protein, add_alias_names, \
    convert_amino_acid_codes, rewrite_nt_snps_as_hgvs, remove_nts_from_nt_name, \
    add_hgvs_names
from functions import empty_attributes, gvf_columns, vcf_columns, pragmas


def vcftogvf(vcf, strain, GENE_PROTEIN_POSITIONS_DICT, sample_size):
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
    vcf_df_cols_to_add = ['vcf_gene', 'mutation_type',
                        'ps_filter', 'ps_exc', 'mat_pep','mat_pep_desc',
                        'mat_pep_acc', 'Reference_seq', 'Variant_seq',
                        "dp", "ro", "ao"]
    for column in list(set(vcf_df.columns) & set(vcf_df_cols_to_add)):
        # drop nans if they exist
        vcf_df[column] = vcf_df[column].fillna('')
        new_gvf[column] = vcf_df[column]
    
    for column in ['nt_name', 'aa_name']:
        vcf_df[column] = vcf_df[column].fillna('n/a')
        new_gvf[column] = vcf_df[column]   

    # add other attributes
    new_gvf['sample_size'] = sample_size
    new_gvf['Name'] = vcf_df["Names"]
    new_gvf['viral_lineage'] = strain
    new_gvf['alternate_frequency'] = vcf_df["AF"]
    
    # add gene and protein attributes from JSON
    json_df = map_pos_to_gene_protein(
        vcf_df['POS'].astype(int), GENE_PROTEIN_POSITIONS_DICT)
    new_gvf["gene"] = json_df["gene"]
    new_gvf["protein_name"] = json_df["protein_name"]
    new_gvf["protein_symbol"] = json_df["protein_symbol"]
    new_gvf["protein_id"] = json_df["protein_id"]
    
    # add 'alias' column for ORF1a/b mutations
    new_gvf = add_alias_names(new_gvf, GENE_PROTEIN_POSITIONS_DICT)

    # add clade_defining attribute
    new_gvf = clade_defining_threshold(args.clades_threshold,
                                             new_gvf, sample_size)
    
    # add HGVS names columns: 'hgvs_nt', 'hgvs_aa', 'hgvs_alias'
    new_gvf = add_hgvs_names(new_gvf)
    
    # save HGVS names for troubleshooting
    #new_gvf[['nt_name', 'hgvs_nt', 'aa_name', 'hgvs_aa', 'alias', 'hgvs_alias']].to_csv('hgvs_troubleshooting.tsv', sep='\t')
    
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
    parser.add_argument('--vcffile', type=str, default=None, required=True,
                        help='Path to a snpEFF-annotated VCF file')
    parser.add_argument('--grouping_criteria', type=str, default=None, required=True,
                        choices=['time', 'wastewater', 'lineage']
                        help="The sample group type")
    parser.add_argument('--sample_group', type=str, default=None, required=True,
                        help='sample group name, ie. lineage, date range')
    parser.add_argument('--size_stats', type=str, default=None,
                        help='Statistics file for for size extraction')
    parser.add_argument('--clades_threshold', type=float,
                        default=0.75,
                        help='Alternate frequency cutoff for '
                             'clade-defining mutations')
    parser.add_argument('--gene_positions', type=str,
                        default=None,
                        help='gene positions in JSON format')
    parser.add_argument('--strain', type=str,
                        default=None,
                        help='Lineage; user mode is if strain="n/a"')
    parser.add_argument("--wastewater", help="Activate wastewater data mode",
                        action="store_true")
    parser.add_argument('--outgvf', type=str, required=True,
                        help='Filename for the output GVF file')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    
    # Reading the gene & protein coordinates of SARS-CoV-2 genome
    with open(args.gene_positions) as fp:
        GENE_PROTEIN_POSITIONS_DICT = json.load(fp)
    
    # Assigning variables
    vcf_file = args.vcffile
    grouping_criteria = args.grouping_criteria
    sample_group = args.sample_group

    # If the strain and/or stats file are None, set them as 'n/a'
    size_stats = args.size_stats
    strain = args.strain
    
    if size_stats == None:
            size_stats='n/a'
    if strain == None:
            strain='n/a'
    

    # print("Processing: " + vcf_file)

    sample_size = find_sample_size(size_stats, strain, vcf_file, args.wastewater)
    
    # create gvf from annotated vcf (ignoring pragmas for now)
    gvf = vcftogvf(vcf_file, strain, GENE_PROTEIN_POSITIONS_DICT,
                   sample_size)
    
    # add species to pragmas
    species = GENE_PROTEIN_POSITIONS_DICT['Src']['species']
    pragmas[0] = pragmas[0].str.replace("##species", "##species " + str(species))
    # temporary pragma, subject to change
    pragmas[0] = pragmas[0].str.replace("##sample-description", "##sample-description " + 'grouping_criteria=' + str(grouping_criteria) +';' + 'sample_group=' + str(sample_group) + ';')

    # combine pragmas, header, GVF contents
    final_gvf = pd.DataFrame(np.vstack([gvf.columns, gvf]))
    final_gvf = pragmas.append(final_gvf)
    
    # save GVF
    filepath = args.outgvf  # outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)

    print("")
    print("Processing complete.")
