#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:47:11 2021

@authors: madeline & zohaib

This script converts GVF files to TSVs, for later conversion to a PDF
case report.

"""
import argparse
import pandas as pd
import os
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts GVF files to a TSV report')
    parser.add_argument('--gvf_files', type=str, default=None,
                        nargs='*', help='Paths to GVF files to process')
    parser.add_argument('--clades', type=str, default=None,
                        help='TSV file of WHO strain names and '
                             'VOC/VOI status')
    parser.add_argument('--outtsv', type=str,
                        default="surveillance_report.tsv",
                        help='Filepath for finished .tsv')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--specify_variants', type=str, default=None,
                       nargs='*', help='Name(s) of WHO variant(s) to '
                                       'make report for. Not '
                                       'case-sensitive.  Space separated.')
    group.add_argument('--all_variants', action="store_true",
                       help='Create reports for all variants, '
                            'using all available reference lineage '
                            'gvf files. Not case-sensitive.')
    parser.add_argument('--table', type=str, default=None,
                        help='Multi-strain TSV file generated in '
                             'workflow that contains num_seqs column')
    parser.add_argument('--user', action="store_true",
                        help='Use user-uploaded file')

    return parser.parse_args()


def match_gvfs_to_who_variant(pango_lineage_list, gvf_files_list):
    matched_files = []
    if len(pango_lineage_list) > 1:
        for lineage in pango_lineage_list:
            if "*" in lineage:
                lineage = lineage.replace("*", "")
            matched_files.extend([i for i in gvf_files_list if
                                  i.startswith(lineage)])

            matched_files = sorted(set(matched_files))
    else:
        for lineage in pango_lineage_list:
        #    matched_files.extend([i for i in gvf_files_list if
        #                          i.startswith(lineage)])
            matched_files = ([i for i in gvf_files_list if
                            i[:i.find("_")]==lineage])

    return matched_files


def find_variant_pop_size(table, pango_lineage_list):
    strain_tsv_df = pd.read_csv(table, header=0,
                                delim_whitespace=True, thousands=r',',
                                usecols=['file', 'num_seqs'])

    files = match_gvfs_to_who_variant(
        pango_lineage_list=pango_lineage_list,
        gvf_files_list=strain_tsv_df['file'].tolist())
    # print(files)

    pop_size = strain_tsv_df.loc[strain_tsv_df['file'].isin(
        files), 'num_seqs'].sum()
    return pop_size


def add_ao_by_variant_seq(ao_str, variant_seq_str):
    """
    Takes a pair of strings like ao_str="7,7,24" and
    variant_seq_str="T,T,T" or "8,21,3" and "T,T,A".
    Output should be a string "T=38" (first case) or "T=29, A=3" (
    second case).
    """
    ao_list = ao_str.split(',')
    ao_list = [int(x) for x in ao_list]
    var_str_list = variant_seq_str.split(',')
    zipped_lists = list(zip(var_str_list, ao_list))

    ao_dict = dict()
    for pair in zipped_lists:
        # if a variant seq isn't in the dictionary, add it
        if pair[0] not in ao_dict:
            ao_dict[pair[0]] = pair[1]
        # if it's already there, add this ao to the existing key
        else:
            ao_dict[pair[0]] += pair[1]

    # create 3 strings: one for ao, one for var_seq, one for both
    # joined together with an '=' between
    joined_string = ''
    ao_string = ''
    var_string = ''
    for item in list(ao_dict.items()):
        joined_string = joined_string + item[0] + '=' + str(item[1]) \
                        + ', '
        ao_string = ao_string + str(item[1]) + ','
        var_string = var_string + item[0] + ','
    joined_string = joined_string.rstrip(", ")
    ao_string = ao_string.rstrip(", ")
    var_string = var_string.rstrip(", ")

    return joined_string, ao_string, var_string


def add_one_from_each_lineage(count_str, lineage_str, mode):
    """
    Takes a pair of strings like count_str="7,7,24" and
    lineage_str="Q.3,Q.3,Q.1".
    Output should be a string "31", adding unique lineage
    values only, if mode='add'.
    If mode='comma', returns a comma-separated set of values,
    like "7,24".
    """
    count_list = count_str.split(';')
    if mode == 'add':
        count_list = [int(x) for x in count_list]
    lineage_str_list = lineage_str.replace(" ", "")
    lineage_str_list = lineage_str.split(',')
    zipped_lists = list(zip(lineage_str_list, count_list))

    count_dict = dict()
    for pair in zipped_lists:
        # if a lineage isn't in the dictionary, add it
        if pair[0] not in count_dict:
            count_dict[pair[0]] = pair[1]
        # if it's already there, do nothing

    # if mode=='comma', return a comma-separated string of values
    if mode == 'comma':
        return_str = ','.join(list(count_dict.values()))
    # elif mode=='add', add up all the keys and return the value as a
    # string
    elif mode == 'add':
        return_str = str(sum(list(count_dict.values())))

    return return_str


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
    df = df.drop(labels=['#source', '#attributes'], axis=1)

    # remove '#' from column names
    df.columns = df.columns.str.replace("#", "")

    # drop unwanted columns
    df = df.drop(labels=['seqid', 'type', 'end',
                         'strand', 'score', 'phase', 'id'], axis=1)

    # rename 'dp' column to 'sequence_depth', make 'viral_lineage'
    # plural
    df = df.rename(columns={'sample_size': 'obs_sample_size',
                            'viral_lineage': 'viral_lineages',
                            'source': 'citation_url'})

    return df


def streamline_user(tsv_df):
    tsv_df = tsv_df.rename(columns={'multiaa_comb_mutation':
                                        'multiaa_mutation_split_names'})
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name',
            'multiaa_mutation_split_names', 'start', 'vcf_gene',
            'chrom_region', 'mutation_type', 'dp', 'obs_sample_size',
            'ps_filter', 'ps_exc', 'mat_pep_id',
            'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao',
            'variant_seq', 'reference_seq', 'function_category',
            'citation', 'citation_url', 'comb_mutation',
            'function_description',
            'heterozygosity']
    tsv_df = tsv_df[cols]

    return tsv_df


def streamline_tsv(tsv_df):
    # find identical rows across strains, and keep only one row.
    # change n/a to 0 in 'ao' for counting purposes
    tsv_df['ao'] = tsv_df['ao'].str.replace("n/a", "0")

    '''
    # make ro, dp, and obs_sample_size numeric
    for colname in ['ro', 'dp', 'obs_sample_size']:
        tsv_df[colname] = pd.to_numeric(tsv_df[colname],
                                        errors='coerce')
    '''
    # make obs_sample_size numeric
    tsv_df['obs_sample_size'] = pd.to_numeric(tsv_df['obs_sample_size'],
                                              errors='coerce')

    agg_dict = dict((col, 'first') for col in
                    tsv_df.columns.values.tolist())

    agg_dict['obs_sample_size'] = 'sum'
    # join some columns with commas
    agg_dict['viral_lineages'] = ', '.join
    agg_dict['clade_defining'] = ','.join
    agg_dict['ao'] = ';'.join
    agg_dict['dp'] = ';'.join
    agg_dict['ro'] = ';'.join
    agg_dict['variant_seq'] = ','.join

    cols_to_check = ['name', 'nt_name', 'aa_name', 'multi_aa_name',
                     'multiaa_comb_mutation', 'start',
                     'function_category', 'citation',
                     'comb_mutation', 'function_description',
                     'heterozygosity']

    final_df = tsv_df.groupby(cols_to_check).agg(agg_dict)
    final_df = final_df.rename(columns={'ao': 'ao_all',
                                        'variant_seq':
                                            'variant_seq_all',
                                        'multiaa_comb_mutation':
                                            'multiaa_mutation_split_names'})

    # add dp, ro per mutation
    for colname in ['dp', 'ro']:
        final_df[colname] = [add_one_from_each_lineage(x, y, mode='add')
                             for x, y in
                             zip(final_df[colname],
                                 final_df['viral_lineages'])]

    # add ao, variant_seq per mutation
    for colname in ['ao_all', 'variant_seq_all']:
        final_df[colname] = [
            add_one_from_each_lineage(x, y, mode='comma') for x, y in
            zip(final_df[colname],
                final_df['viral_lineages'])]

    # add ao according to the heterogeneous mutations
    final_df['ao_by_var_seq'] = [add_ao_by_variant_seq(x, y)[0] for x, y
                                 in zip(final_df['ao_all'],
                                        final_df['variant_seq_all'])]
    final_df['ao'] = [add_ao_by_variant_seq(x, y)[1] for x, y in
                      zip(final_df['ao_all'],
                          final_df['variant_seq_all'])]
    final_df['variant_seq'] = [add_ao_by_variant_seq(x, y)[2] for x, y
                               in zip(final_df['ao_all'],
                                      final_df['variant_seq_all'])]

    # remove 'who_variant'; rename 'multiaa_comb_mutation'
    final_df = final_df.drop(labels=['variant'], axis=1)
    # add variant_pop_size
    final_df['variant_pop_size'] = variant_pop_size

    # combine viral_lineages and clade_defining into key-value pairs

    # split viral_lineages and clade_defining by ','
    split_lineages = final_df['viral_lineages'].str.split(
        pat=',').apply(pd.Series)  # split at ,, form dataframe
    split_clade_defining = final_df['clade_defining'].str.split(
        pat=',').apply(pd.Series)  # split at ,, form dataframe
    # go through and make key-value pairs of corresponding columns
    # from each
    final_df['clade_defining_status'] = ''
    for col in split_clade_defining.columns:
        final_df['clade_defining_status'] = final_df[
                                                'clade_defining_status'] + \
                                            split_lineages[col].astype(
                                                str) + '=' + \
                                            split_clade_defining[
                                                col].astype(str) + '; '
    # drop clade_defining status for n/a strains and empty nan=nan pairs
    final_df.clade_defining_status = \
        final_df.clade_defining_status.str.replace('n/a=n/a; ',
                                                   'n/a; ')
    final_df.clade_defining_status = \
        final_df.clade_defining_status.str.replace('nan=nan; ', '')
    # strip trailing spaces and semicolons
    final_df.clade_defining_status = \
        final_df.clade_defining_status.str.rstrip("; ")

    # drop repeated key-value pairs in each row (find these rows as
    # they contain spaces)
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
            mask = final_df['clade_defining_status'] == row
            final_df.loc[mask, 'clade_defining_status'] = row_str

    # drop repeated lineage names in each row of viral_lineages
    # return an ordered list of lineages
    for row in final_df['viral_lineages']:
        if ' ' in row:
            lineage_list = row.split(', ')
            mylist = list(set(lineage_list))
            '''
            #order lineages alphanumerically
            #split each element of mylist into a sublist, split at
            the first '.'
            split_list =
            '''
            row_str = ', '.join(str(e) for e in mylist)
            mask = final_df['viral_lineages'] == row
            final_df.loc[mask, 'viral_lineages'] = row_str

    # reorder columns
    cols = ['name', 'nt_name', 'aa_name', 'multi_aa_name',
            'multiaa_mutation_split_names', 'start', 'vcf_gene',
            'chrom_region', 'mutation_type', 'dp', 'obs_sample_size',
            'variant_pop_size', 'ps_filter', 'ps_exc', 'mat_pep_id',
            'mat_pep_desc', 'mat_pep_acc', 'ro', 'ao_by_var_seq', 'ao',
            'variant_seq', 'reference_seq', 'function_category',
            'citation',
            'citation_url', 'comb_mutation', 'function_description',
            'heterozygosity', 'viral_lineages',
            'clade_defining_status', 'status',
            'voi_designation_date', 'voc_designation_date',
            'vum_designation_date']
    final_df = final_df[cols]

    return final_df


if __name__ == '__main__':

    args = parse_args()
    outfile = args.outtsv
    gvf_list = args.gvf_files



    if not args.user:
        clade_file = args.clades
        # read in WHO variant/PANGO lineage .tsv
        clades = pd.read_csv(clade_file, sep='\t', header=0, usecols=[
            'variant', 'pango_lineage'])

        # get lowercase WHO variant names
        who_variants_list = []
        if args.specify_variants:
            who_variants_list = args.specify_variants
        elif args.all_variants:
            # if all variants, read the file and add from who_variant
            # column, in case of Variant under monitoring, lineage will
            # be added as the variant name to avoid several lineages
            # being combined as one variant
            for i in range(0, len(clades['variant'])):
                if clades.loc[i, 'variant'] == "Unnamed":
                    who_variants_list.append(clades.loc[i,
                                                        'pango_lineage'])

                else:
                    who_variants_list.append(clades.loc[i,
                                                        'variant'])

        # for each variant, create a surveillance report
        for who_variant in who_variants_list:
            if "_" in who_variant:
                who_variant = who_variant[0:who_variant.find(
                "_")].capitalize() + who_variant[who_variant.find(
                "_"):]
            else:
                who_variant = who_variant.capitalize()
            
            # get list of relevant pango lineages
            pango_lineages = []
            for var in clades[clades['variant']==who_variant]['pango_lineage'].tolist():
                if "," in var:
                    for temp in var.split(","):
                        if "[" not in var:
                            pango_lineages.append(temp)
                        else:
                            parent = temp[0]
                            child = temp[2:-3].split("|")
                            for c in child:
                                pango_lineages.append(parent + str(c))
                                pango_lineages.append(parent + str(c) + ".*")
                else:
                    pango_lineages.append(var)

            #print(pango_lineages)

            # get list of gvf files pertaining to variant
            gvf_files = []
            gvf_files = match_gvfs_to_who_variant(
                pango_lineage_list=pango_lineages,
                gvf_files_list=gvf_list)
            print(str(len(gvf_files)) + " GVF files found for " +
                  who_variant + " variant.")


            if len(gvf_files) > 0:
                # convert all gvf files to tsv and concatenate them
                print("")
                print("Processing:")
                print(gvf_files[0])
                gvf_df = gvf2tsv(gvf=gvf_files[0])
                if len(gvf_files) > 1:
                    for gvf in gvf_files[1:]:
                        print(gvf)
                        new_gvf_df = gvf2tsv(gvf=gvf)
                        gvf_df = pd.concat([gvf_df, new_gvf_df],
                                           ignore_index=True)


                if args.table:
                    # get variant population size
                    variant_pop_size = find_variant_pop_size(table=args.table,
                                                             pango_lineage_list=
                                                             pango_lineages)

                out_df = streamline_tsv(tsv_df=gvf_df)
                filename = who_variant + '_' + outfile
                out_df.to_csv(filename, sep='\t', index=False)
                print("Processing complete.")
                print(who_variant + " surveillance report saved as: " +
                      filename)
                print("")


    # if user-provided, who_variant is the provided filename
    else:
        gvf_files = gvf_list
        if "_qc" not in gvf_list[0]:
            who_variant = gvf_list[0].replace(
                ".filtered.SNPEFF.annotated.gvf", "")
        else:
            who_variant = gvf_list[0].replace(
                "_qc.sorted.variants.normalized.filtered.SNPEFF"
                ".annotated.gvf", "")

        if args.table:
            # get variant population size
            variant_pop_size = find_variant_pop_size(table=args.table,
                                                     pango_lineage_list=
                                                     pango_lineages)
        else:
            variant_pop_size = "n/a"

    # if any GVF files are found, create a surveillance report


    # convert all gvf files to tsv and concatenate them
        print("Processing:")
        print(gvf_files[0])
        gvf_df = gvf2tsv(gvf=gvf_files[0])

        # streamline final concatenated df, reorder/rename
        # columns where needed
        out_df = streamline_user(tsv_df=gvf_df)

            # save report as a .tsv
        filename = who_variant + '_' + outfile
        out_df.to_csv(filename, sep='\t', index=False)
        print("Processing complete.")
        print(who_variant + " surveillance report saved as: " +
              filename)
        print("")
