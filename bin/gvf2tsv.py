#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os

"""
Created on Mon Aug  9 10:47:11 2021

@author: madeline
"""

'''This script converts GVF files to TSVs, for later conversion to an 
HTML case report. '''


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a GVF file to a TSV')
    parser.add_argument('--gvf_files', type=str, default=None,
                        nargs='*',
                        help='Path to GVF-containing directory')
    parser.add_argument('--clades', type=str, default=None,
                        help='TSV file of WHO strain names and '
                             'VOC/VOI status')
    parser.add_argument('--outtsv', type=str,
                        default="surveillance_report",
                        help='Filepath for finished .tsv')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--specify_variants', type=str, default=None,
                       nargs='*', help='Name(s) of WHO variant to '
                                       'make report for (Not '
                                       'case-sensitive).')
    group.add_argument('--all_variants', action="store_true",
                       help='Create reports for all variants, using '
                            'all available reference lineage gvf '
                            'files')

    return parser.parse_args()


def find_gvfs(lineages, gvf_list):
    """This function searches for gvf files in a directory which have
    filenames that contain specified lineage names.
    :param lineages: Pangolin lineages to include in the report
    :type lineages: list of strings
    :param gvf_list: directory to search for gvf files
    :type gvf_list: string
    :return: list of gvf filenames to process for the report
    """
    file_list = []
    if not "." in lineages:
        for pango_lineage in lineages:
            for f in gvf_list:
                if f.startswith(pango_lineage.replace("*", "")) and \
                        f.endswith(".gvf"):
                    file_list.append(f)
    else:
        for f in gvf_list:
            if lineages in f and \
                    f.endswith(".gvf"):
                file_list.append(f)

    return file_list


def gvf2tsv(gvf):
    # read in gvf
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']

    df = pd.read_csv(gvf, sep='\t', names=gvf_columns)
    # remove pragmas and original header
    df = df[~df['#seqid'].str.contains("#")]
    # restart index from 0
    df = df.reset_index(drop=True)

    # split #attributes column into separate columns for each tag
    # split at ;, form dataframe
    attributes = df['#attributes'].str.split(pat=';').apply(pd.Series)
    # last column is a copy of the index so drop it
    attributes = attributes.drop(labels=len(attributes.columns) - 1,
                                 axis=1)

    for column in attributes.columns:
        split = attributes[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0].lower()
        content = split[1]
        # ignore "tag=" in column content
        attributes[column] = content
        # make attribute tag as column label
        attributes.rename(columns={column: title}, inplace=True)

    # replace attributes column in the original df with the new
    # separated out attributes
    df = pd.concat((df, attributes), axis=1)
    df = df.drop(labels='#attributes', axis=1)

    # remove '#' from column names
    df.columns = df.columns.str.replace("#", "")

    # drop unwanted columns
    df = df.drop(labels=['source', 'seqid', 'type', 'end', 'strand',
                         'score', 'phase', 'id'], axis=1)

    # rename 'dp' column to 'sequence_depth', make 'viral_lineage'
    # plural
    df = df.rename(columns={'dp': 'sequence_depth',
                            'viral_lineage': 'viral_lineages'})

    return df


def streamline_tsv(df):
    # find identical rows across strains, and keep only one row.

    # change n/a to 0 in 'ao' for counting purposes
    df['ao'] = df['ao'].str.replace("n/a", "0")

    for colname in ['sequence_depth', 'ao', 'ro']:
        # split up at commas into new columns: make a new mini-df
        split_series = df[colname].str.split(pat=',').apply(
            pd.Series)
        # rename series columns to 'ao_0', 'a0_1', etc.
        split_series.columns = [colname + '_' + str(name) for name in
                                split_series.columns.values]
        # ensure all counts are numeric
        for column in split_series.columns:
            split_series[column] = pd.to_numeric(split_series[
                                                     column],
                                                 errors='coerce')
        # append series to df
        df = pd.concat([df, split_series], axis=1)

    cols_to_check = ['name', 'nt_name', 'aa_name', 'multi_aa_name',
                     'multiaa_comb_mutation', 'start',
                     'function_category', 'citation',
                     'comb_mutation', 'function_description',
                     'heterozygosity']

    agg_dict = dict((col, 'first') for col in
                    df.columns.values.tolist())
    agg_dict['viral_lineages'] = ', '.join
    agg_dict['clade_defining'] = ', '.join

    # sum split columns
    for string in ['sequence_depth_', 'ao_', 'ro_']:
        relevant_keys = [key for key, value in agg_dict.items() if
                         string in key.lower()]
        for key in relevant_keys:
            agg_dict[key] = 'sum'

    final_df = df.groupby(cols_to_check).agg(agg_dict)

    # rejoin split columns (ao, sequence_depth, ro) with comma
    # separation
    for string in ['sequence_depth_', 'ao_', 'ro_']:
        colnames = [i for i in df.columns.values.tolist() if
                    string in i]
        final_df[string + 'combined'] = final_df[colnames].apply(
            lambda row: ','.join(row.values.astype(str)), axis=1)
        # drop split columns
        final_df = final_df.drop(labels=colnames, axis=1)

    # replace ao, ro, sequence_depth with the added up columns
    final_df = final_df.drop(labels=['ao', 'ro', 'sequence_depth'],
                             axis=1)
    final_df = final_df.rename(columns={'sequence_depth_combined':
                                            'sequence_depth',
                                        'ro_combined': 'ro',
                                        'ao_combined': 'ao'})

    # reorder columns
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name',
            'multiaa_comb_mutation', 'start', 'vcf_gene',
            'chrom_region',
            'mutation_type', 'sequence_depth', 'sample_size',
            'ps_filter', 'ps_exc', 'mat_pep_id',
            'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao', 'reference_seq',
            'variant_seq', 'viral_lineages', 'function_category',
            'citation',
            'comb_mutation', 'function_description', 'heterozygosity',
            'clade_defining', 'who_variant', 'status',
            'voi_designation_date', 'voc_designation_date',
            'vum_designation_date']
    final_df = final_df[cols]

    return final_df


if __name__ == '__main__':

    args = parse_args()

    filepath = args.outtsv
    clade_file = args.clades

    # gvf_directory = args.gvf_directory #directory to search
    gvf_files_list = args.gvf_files

    # read in WHO variant/PANGO lineage .tsv
    clades = pd.read_csv(clade_file, sep='\t', header=0,
                         usecols=['who_variant', 'pango_lineage'])

    # get lowercase WHO variant names
    who_variants_list = []
    if args.specify_variants:
        who_variants_list = args.specify_variants
    elif args.all_variants:
        for i in range(0, len(clades['who_variant'])):
            if not clades.loc[i, 'who_variant'] == "Unnamed":
                who_variants_list.append(clades.loc[i, 'who_variant'])
            else:
                who_variants_list.append(clades.loc[i, 'pango_lineage'])

    print(who_variants_list)

    # for each variant called, create a surveillance report
    for who_variant in sorted(who_variants_list):
        who_variant = who_variant.capitalize()
        # get list of relevant pango lineages
        # list of pango lineages from that variant
        if not "." in who_variant:
            pango_lineages = \
                clades[clades['who_variant'] == who_variant][
                    'pango_lineage'].values[0].split(',')
        else:
            pango_lineages = who_variant

        # get list of gvf files pertaining to variant
        gvf_files = find_gvfs(lineages=pango_lineages,
                              gvf_list=gvf_files_list)
        print(str(len(gvf_files)) + " GVF files found for " +
              who_variant + " variant.")

        # if any GVF files are found, create a surveillance report
        if len(gvf_files) > 0:

            # convert all gvf files to tsv and concatenate them
            print("Processing:")
            print(gvf_files[0])
            tsv_df = gvf2tsv(gvf_files[0])
            for gvf_file in gvf_files[1:]:
                print(gvf_file)
                new_tsv_df = gvf2tsv(gvf=gvf_file)
                tsv_df = pd.concat([tsv_df, new_tsv_df],
                                   ignore_index=True)

            # streamline final concatenated df, reorder/rename
            # columns where needed
            surveillance_df = streamline_tsv(df=tsv_df)

            # save report as a .tsv
            filename = filepath + '_' + who_variant + '.tsv'
            surveillance_df.to_csv(filename, sep='\t', index=False)
            print("Processing complete.")
            print(who_variant + " surveillance report saved as: " +
                  filename)
            print("")
