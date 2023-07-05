#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 08:25:48 2023

@author: madeline

This script annotates GVF files with the functional annotation.

"""
import argparse
import pandas as pd
import numpy as np
from functions import separate_attributes, rejoin_attributes


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
    parser.add_argument('--strain', type=str,
                        default='n/a',
                        help='Lineage; user mode is if strain="n/a"')
    parser.add_argument("--names", help="Save mutation names without "
                                        "functional annotations to "
                                        "TSV files for "
                                        "troubleshooting purposes",
                        action="store_true")
    return parser.parse_args()


def add_pokay_annotations(gvf, annotation_file, strain):
    
    # expand #attributes into columns to fill in separately
    gvf = separate_attributes(gvf)
    
    # drop columns that are going to be re-added in the merge
    functional_attributes = ["function_category", "function_description", 
                             "source", "citation", "comb_mutation", 
                             "heterozygosity"]
    gvf = gvf.drop(columns=functional_attributes)

    # remember this includes nucleotide names where there are no
    # protein names
    gvf["mutation"] = gvf["Name"]

    # merge annotated vcf and functional annotation files by
    # 'mutation' column in the gvf
    # load functional annotations spreadsheet
    df = pd.read_csv(annotation_file, sep='\t', header=0)

    for column in df.columns:
        df[column] = df[column].str.lstrip()
        # add functional annotations
    merged_df = pd.merge(df, gvf, on=['mutation'], how='right')

    # collect all mutation groups (including reference mutation) in a
    # column, sorted alphabetically
    # this is more roundabout than it needs to be; streamline with
    # grouby() later
    merged_df["mutation_group"] = merged_df["comb_mutation"].astype(
        str) + ", '" + merged_df["mutation"].astype(str) + "', " + \
                                  merged_df[
                                      'multiaa_comb_mutation'].astype(
                                      str)
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.replace("nan, ", "")
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.rstrip(' ').str.rstrip(',')

    # separate the mutation_group column into its own df with one
    # mutation per column
    mutation_groups = merged_df["mutation_group"].str.split(
        pat=',').apply(pd.Series)
    mutation_groups = mutation_groups.apply(
        lambda s: s.str.replace("'", ""))
    mutation_groups = mutation_groups.apply(
        lambda s: s.str.replace(" ", ""))
    # now each mutation has a column instead
    mutation_groups = mutation_groups.transpose()
    # sort each column alphabetically
    sorted_df = mutation_groups

    for column in mutation_groups.columns:
        sorted_df[column] = mutation_groups.sort_values(by=column,
                                                        ignore_index=True)[
            column]
    sorted_df = sorted_df.transpose()

    # since they're sorted, put everything back into a single cell,
    # don't care about dropna
    df3 = sorted_df.apply(lambda x: ','.join(x.astype(str)), axis=1)
    unique_groups = df3.drop_duplicates()
    unique_groups_multicol = sorted_df.drop_duplicates()
    # for sanity checking
    merged_df["mutation_group_labeller"] = df3

    # make a unique id for mutation groups that have all members
    # represented in the vcf
    # for groups with missing members, delete those functional
    # annotations
    merged_df["id"] = 'NaN'
    id_num = 0
    for row in range(unique_groups.shape[0]):
        group_mutation_set = set(unique_groups_multicol.iloc[row])
        # remove nan and 'nan' from set
        group_mutation_set = {x for x in group_mutation_set if (x == x
                                                                and x != 'nan')}
        gvf_all_mutations = set(gvf['mutation'].unique())
        indices = merged_df[merged_df.mutation_group_labeller ==
                            unique_groups.iloc[row]].index.tolist()
        # if all mutations in the group are in the vcf file, include
        # those rows and give them an id
        if group_mutation_set.issubset(gvf_all_mutations):
            merged_df.loc[merged_df.mutation_group_labeller ==
                          unique_groups.iloc[row], "id"] = "ID_" + \
                                                           str(id_num)
            id_num += 1
        else:
            # if not, drop group rows, leaving the remaining indices
            # unchanged
            merged_df = merged_df.drop(indices)

            # change semicolons in function descriptions to colons
    merged_df['function_description'] = merged_df[
        'function_description'].str.replace(';', ':')
    # change heteozygosity column to True/False
    merged_df['heterozygosity'] = merged_df['heterozygosity'] == \
                                  'heterozygous'
    # remove trailing spaces from citation
    merged_df['citation'] = merged_df['citation'].str.strip()
    # add key-value pairs to attributes column
    for column in ['function_category', 'source', 'citation',
                   'comb_mutation', 'function_description',
                   'heterozygosity']:
        key = column.lower()
        # replace NaNs with empty string
        merged_df[column] = merged_df[column].fillna('')
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + key + '=' + merged_df[column].astype(str) + ';'
    
    # add ID to attributes
    merged_df["ID"] = merged_df['id'].astype(str)
    
    # merge attributes back into a single column
    merged_df = rejoin_attributes(merged_df, empty_attributes)

    return merged_df[gvf_columns]


if __name__ == '__main__':

    args = parse_args()
    
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']
    
    pragmas = pd.DataFrame([['##gff-version 3'], ['##gvf-version '
                                                 '1.10'], [
                                '##species NCBI_Taxonomy_URI=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049']])  # pragmas are in column 0

    empty_attributes = 'ID=;Name=;chrom_region=;protein=;ps_filter=;ps_exc=; \
        mat_pep_id=;mat_pep_desc=;mat_pep_acc=; ro=;ao=;dp=;sample_size=; \
        Reference_seq=;Variant_seq=;nt_name=;aa_name=;vcf_gene=; \
        mutation_type=; viral_lineage=;multi_aa_name=;multiaa_comb_mutation=; \
        alternate_frequency=;function_category=;source=; citation=; \
        comb_mutation=;function_description=;heterozygosity=;clade_defining=; \
        variant=;variant_type=;voi_designation_date=;voc_designation_date=; \
        vum_designation_date=;status=;'
    empty_attributes = empty_attributes.replace(" ", "")                                              


    # read in gvf file
    gvf = pd.read_csv(args.ingvf, sep='\t', names=gvf_columns, index_col=False)

    # remove pragmas and original header row
    gvf = gvf[~gvf['#seqid'].astype(str).str.contains("#")]

    # separate pragmas
    #pragmas = gvf[gvf['#seqid'].str.contains("##")]
        
    # add functional annotations
    pokay_annotated_gvf = add_pokay_annotations(gvf, args.functional_annotations, args.strain)

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
    if args.names:
        # save unmatched names (in vcf/tsv but not in
        # functional_annotations) to a .tsv file
        
        # create mask to find which rows do not have a functional annotation
        notinPokay_mask = pokay_annotated_gvf["#attributes"].str.contains("function_category=;")
        # extract all mutation names from #attributes column
        names = pd.Series(pokay_annotated_gvf["#attributes"].str.findall('(?<=Name=)(.*?)(?=;)').str[0])
        # get unique mutation names not in Pokay
        unmatched_names = pd.Series(names[notinPokay_mask].unique())
        # save unmatched names to TSV
        if unmatched_names.shape[0] != 0:
            leftover_names_filepath = "unmatched_names.tsv"
            unmatched_names.to_csv(leftover_names_filepath, sep='\t',
                                   index=False, header=False)
            print("")
            print(str(unmatched_names.shape[0]) +
                  " mutation names not matched with functional annotations "
                  "file saved to " + leftover_names_filepath)

    print("")
    print("Processing complete.")
