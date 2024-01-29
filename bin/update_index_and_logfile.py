#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Given a new batch of GVFs and a new partial logfile, adds the new mutation and lineage info to the mutation index and produces a new logfile.
If no mutation index is provided, a new one is created from scratch using the GVFs.

--partial_logfile is the path to a text file of the form:

Total sequences:	508910
New sequences:	171
New lineages:	['HH.1.1', 'FT.3.1.1']
'''

import pandas as pd
import numpy as np
import argparse
import glob
import os
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Creates a list of all unique mutation names across all GVFs')
    parser.add_argument('--mutation_index', type=str, default=None,
                        help='Path to the TSV index of all mutations')
    parser.add_argument('--gvf_files', type=str, default=None,
                        nargs='*', help='Paths to GVF files to process')
    parser.add_argument('--index_savefile', type=str,
                        default=None, help='Filename to save updated index of all mutations to')
    parser.add_argument('--partial_logfile', type=str,
                        default=None, help='Path to partial log file, to append new mutations to')
    parser.add_argument('--log_savefile', type=str,
                        default=None, help='Filename to save log to')

    return parser.parse_args()


def gvf2df(gvf):

    # read in gvf
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']
    gvf = pd.read_csv(gvf, sep='\t', names=gvf_columns, usecols=['#start', '#seqid', '#attributes'])
    # remove pragmas and original header
    gvf = gvf[~gvf['#seqid'].str.contains("#")]

    # create new dataframe in the format of the index, but just for this one gvf
    df = pd.DataFrame(np.empty((gvf.shape[0], 9)), columns=index_cols)
    df['pos'] = gvf['#start'].tolist()
    df['mutation'] = gvf["#attributes"].str.findall('(?<=Name=)(.*?)(?=;)').str[0].tolist()
    df['mutation'] = df['mutation'].str.replace("p.", "", regex=False)
    df['hgvs_aa_mutation'] = gvf["#attributes"].str.findall('(?<=hgvs_aa=)(.*?)(?=;)').str[0].tolist()
    df['hgvs_nt_mutation'] = gvf["#attributes"].str.findall('(?<=hgvs_nt=)(.*?)(?=;)').str[0].tolist()
    #df['alias'] = gvf["#attributes"].str.findall('(?<=alias=)(.*?)(?=;)').str[0].tolist()
    df['hgvs_alias'] = gvf["#attributes"].str.findall('(?<=hgvs_alias=)(.*?)(?=;)').str[0].tolist()
    df['gene'] = gvf["#attributes"].str.findall('(?<=gene=)(.*?)(?=;)').str[0].tolist()
    df['protein_name'] = gvf["#attributes"].str.findall('(?<=protein_name=)(.*?)(?=;)').str[0].tolist()
    df['Pokay_annotation'] = (~gvf["#attributes"].str.contains("function_description=;")).tolist()
    df['lineages'] = gvf["#attributes"].str.findall('(?<=lineage=)(.*?)(?=;)').str[0].tolist()
    lineage = df['lineages'][0]
    df = df.astype(str)
    df = df.drop_duplicates()

    return df, lineage

if __name__ == '__main__':

    args = parse_args()
    gvf_list = args.gvf_files
    mutation_index_path = args.mutation_index
    index_savefile = args.index_savefile
    partial_logfile = args.partial_logfile
    log_savefile = args.log_savefile
    
    # set empty structures for creating the logfile at the end
    lineages = []
    logfile_df = pd.DataFrame(np.empty((0, 4)), columns=['pos', 'alias', 'new_mutations', 'lineages'])

    # set mutation index columns
    index_cols=['pos', 'mutation', 'hgvs_aa_mutation', 'hgvs_nt_mutation', 'hgvs_alias', 'gene', 'protein_name', 'Pokay_annotation', 'lineages']

    # open the mutation index if the path was provided
    if mutation_index_path!=None:
        mutation_index = pd.read_csv(mutation_index_path, sep='\t', dtype='str')
    # if no path was given, create a new mutation index from scratch
    else:
        mutation_index = pd.DataFrame(np.empty((0, 9)), columns=index_cols)
    mutation_index['pos'] = mutation_index['pos'].astype(int)


    # for troubleshooting: make practice mutation_index with all mentions of ['HH.1.1', 'FT.3.1.1'] removed
    #for lineage in ['HH.1.1', 'FT.3.1.1']:
    #    mutation_index = mutation_index[mutation_index['lineage']!=lineage]
    #    mutation_index['lineage'] = mutation_index['lineage'].str.replace(lineage, "")


    # iterate through gvfs in the list argument and get set from each
    for file in gvf_list:
        print("Processing: " + file)

        # open gvf and reformat to match the mutation index
        df, lineage = gvf2df(file)
        df['pos'] = df['pos'].astype(int)
        lineages.append(lineage)
        # append the new gvf df to the index and use groupby() to add the lineage where 
        mutation_index = pd.concat([mutation_index, df])
        mutation_index['hgvs_alias'] = mutation_index['hgvs_alias'].astype(str)
        # groupby group_cols, adding new lineages to "lineages" in a list
        group_cols = ['pos', 'mutation', 'hgvs_aa_mutation', 'hgvs_nt_mutation', 'hgvs_alias', 'gene', 'protein_name', 'Pokay_annotation']

        mutation_index = mutation_index.groupby(by=group_cols, as_index=False).agg(list)
        # convert lineage column from lists back to strings
        mutation_index["lineages"] = [','.join(map(str, l)) for l in mutation_index["lineages"]]

        # extract novel mutations to a new dataframe to save to the logfile
        # novel mutations are found in the rows where "lineage"==gvf_lineage 
        lineage_mask = ((mutation_index['lineages']==lineage) | mutation_index['lineages'].str.contains("," + lineage))
        lineage_chunk = mutation_index.loc[lineage_mask,:].copy()
        lineage_chunk['new_mutations'] = lineage_chunk["gene"] + ":" + lineage_chunk["mutation"]
        # for ORF1ab mutations, add the alias
        orf1ab_mask = lineage_chunk['gene'].astype(str).str.contains("ORF1ab")
        lineage_chunk.loc[orf1ab_mask, 'new_mutations'] = lineage_chunk['new_mutations'] + " / " + lineage_chunk["protein_name"] + ":" + lineage_chunk["hgvs_alias"]
        # for intergenic mutations, replace nan with "intergenic"
        intergenic_mask = lineage_chunk['protein_name'].astype(str).str.contains("nan")
        lineage_chunk.loc[intergenic_mask, 'new_mutations'] =  lineage_chunk.loc[intergenic_mask, 'new_mutations'].str.replace("nan", "intergenic", regex=True)
        # drop duplicates (there shouldn't be any)
        lineage_chunk = lineage_chunk[['pos', 'hgvs_alias', 'new_mutations', 'lineages']].drop_duplicates()
        # drop NaN rows
        lineage_chunk = lineage_chunk[lineage_chunk['pos'].notna()]
        # append to logfile_df
        logfile_df = pd.concat([logfile_df, lineage_chunk])
        logfile_df['pos'] = logfile_df['pos'].astype(int)
            

    # save updated mutation index
    mutation_index.to_csv(index_savefile, sep='\t', header=True, index=False)

    # clean up logfile_df: restrict df to rows with new mutations, added to the index from the GVFs added above
    gvf_lineage_set = lineages
    # get list of lists for "lineages" column
    lineages_col = logfile_df["lineages"].str.split(",").tolist()
    # keep only rows that have only lineages from the new GVFs
    new_mutation_check = [set(subl).issubset(gvf_lineage_set) for subl in lineages_col]
    logfile_df["new_mutation_check"] = new_mutation_check
    logfile_df = logfile_df[logfile_df['new_mutation_check']==True]
    logfile_df = logfile_df[['new_mutations', 'lineages']]

    # read in partial logfile, append the mutations list, and save to the new name
    partial_logfile_df = pd.read_csv(partial_logfile, sep='\t', names=["new_mutations","lineages"])
    partial_logfile_df.loc[len(partial_logfile_df)] = ["New mutations:", ""]
    logfile_df = pd.concat([partial_logfile_df, logfile_df])

    # save logfile
    logfile_df.to_csv(log_savefile, sep='\t', header=False, index=False)
