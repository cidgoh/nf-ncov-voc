#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: anoosha & madeline

This script converts Pokay .txt files into a functional annotation
key (tsv) files that can be used by nf-ncov-voc workflow for
annotation of variant called files.

"""

import os
import argparse
from pathlib import Path
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script produces a TSV file from TXT files '
                    'in POKAY  '
                    'https://github.com/nodrogluap/pokay/tree/master'
                    '/data for annotating SARS-COV-2 mutations')
    parser.add_argument('--inputdir', type=str, default=None,
                        help='directory path for input files')
    parser.add_argument('--accession', type=str, default=None,
                        help='versioned reference accession from RefSeq')
    parser.add_argument('--mutation_index', type=str, default=None,
                        help='index of all mutations (.TSV)')
    parser.add_argument('--outputfile', type=str, default=None,
                        help='output file (.TSV) format')
    return parser.parse_args()


def combination_mutation(c_mutations, c_mutation):
    comb_mutation = []
    for y in range(len(c_mutations)):
        if c_mutations[y] != c_mutation:
            comb_mutation.append(c_mutations[y])
    return comb_mutation


def data_cleanup(dframe):
    dframe['mutation functional effect description'] = dframe[
        'mutation functional effect description'].apply(
        lambda x: ','.join(map(str, x)))

    dframe['mutation functional effect description'] = dframe[
        'mutation functional effect description'].str.replace(',#', '')
    dframe['mutation functional effect description'] = dframe[
        'mutation functional effect description'].str.replace('#', '')
    dframe['mutation functional effect description'] = dframe[
        'mutation functional effect description'].str.strip()
    
    dframe['comb_mutation'] = dframe['comb_mutation'].apply(
        lambda x: x[1:-1])

    return dframe


def extract_source_citation(dframe):
    '''
    Fill in 'publication year', 'author', 'peer review status', 'DOI', and 'URL' columns
    '''
    if '|' in dframe['url']:
        dframe['citation'] = dframe.url.str.split(
            "|").str[0] + '|' + dframe.source.str.split('|')[1]
        dframe['URL'] = dframe.url.str.split(
            "|").str[1]

    else:
        dframe['citation'] = dframe.url.str.split(
            "http").str[0]
        dframe['citation'] = dframe['citation'].str.replace('#',
                                                            '')
        dframe['URL'] = dframe.url.str.split(
            ")").str[1]
        
    dframe['publication year'] = dframe['citation'].str.extract('.*\((.*)\).*')
    dframe['author'] = dframe['citation'].str.extract('^(.*?)et al') #.str.strip()
    dframe['peer review status'] = dframe['URL'].str.extract('(?<=\[)(.*)')
    dframe['URL'] = dframe['URL'].str.extract('^(.*?)\[') #.str.strip()
    dframe['DOI'] = dframe['URL'].str.extract('(?<=doi.org/)(.*)') #.str.strip()

    return dframe


def extract_metadata(inp_file, chunk, df):
    mutation_name = chunk[-1].strip()
    check_combination = 0
    if ";" in mutation_name:
        mutation_name = mutation_name.split(";")
        check_combination = 1
    elif "," in mutation_name:
        mutation_name = mutation_name.split(",")
    else:
        mutation_name = [mutation_name]
    for x in mutation_name:
        heterozygosity = ""
        if x.startswith("[") or x.endswith("]"):
            heterozygosity = "heterozygous"
            x = x.replace("[", "")
            x = x.replace("]", "")
        if check_combination == 1:
            comb_mutation = combination_mutation(
                c_mutations=mutation_name, c_mutation=x)
        else:
            comb_mutation = ""
        gene_name = inp_file.split('_')[0]
        function_category = Path(inp_file).stem
        function_category = function_category.split('_', 1)[1].replace(
            '_', ' ')
        url = []
        for i in range(0, len(chunk)):
            chunk[i] = chunk[i].strip("\n")
            if "http" in chunk[i]:
                url.append(i)
        function = {}
        for index_url in range(0, len(url)):
            if index_url == 0:
                function[chunk[url[index_url]]] = chunk[
                                                  0:url[index_url]]
            else:
                if len(chunk[
                       url[index_url - 1] + 1: url[index_url]]) > 0:
                    function[chunk[url[index_url]]] = chunk[
                                                      url[index_url -
                                                          1] + 1: url[
                                                          index_url]]
                else:
                    new_key = str(chunk[url[index_url - 1]]).strip() \
                              + " | " + chunk[url[index_url]]
                    function[new_key] = chunk[
                                        url[index_url - 2] + 1: url[
                                            index_url - 1]]
                    del function[chunk[url[index_url - 1]]]

        df_func = pd.DataFrame(function.items(), columns=['url',
                                                          'mutation functional effect description'])

        df_list = [x, gene_name,
                   function_category, comb_mutation] #, heterozygosity]
        # print(df_list)
        df1 = pd.DataFrame(
            columns=['original mutation description', 'protein symbol', 'mutation functional effect category',
                     'comb_mutation']) #, 'heterozygosity'])
        df1.loc[len(df1)] = df_list

        df_func['original mutation description'] = df1['original mutation description'].iloc[0]
        df_func['protein symbol'] = df1['protein symbol'].iloc[0]
        df_func['mutation functional effect category'] = df1['mutation functional effect category'].iloc[0]
        df_func['comb_mutation'] = str(df1['comb_mutation'].iloc[0])
        #df_func['heterozygosity'] = str(df1['heterozygosity'].iloc[0])

        df_func = data_cleanup(dframe=df_func)
        df_func = extract_source_citation(dframe=df_func)

        df = pd.concat([df, df_func], ignore_index=True)
        df = df.drop(labels='url', axis=1)
    return df


def write_tsv(dframe):
    dframe.to_csv(args.outputfile, sep="\t", escapechar='|',
                  quoting=csv.QUOTE_ALL, index=False, header=True)


if __name__ == '__main__':

    args = parse_args()

    # Folder Path
    path = args.inputdir

    # Change the directory
    # os.chdir(path)


    dataFrame_cols = ['organism', 'reference accession', 'reference database name', 'nucleotide position',
'original mutation description', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias',
'gene name', 'gene symbol', 'protein name', 'protein symbol', 'assay', 'mutation functional effect category',
'mutation functional effect description', 'author', 'publication year', 'URL', 'DOI', 'PMID',
'peer review status', 'curator', 'mutation functional annotation resource']
    
    dataFrame = pd.DataFrame(columns=dataFrame_cols)

    for file in os.listdir(path):
        if file.endswith("binding.txt") and "_" in file:
            file_path = os.path.join(path, file)
            print(file)
            f = open(file_path, 'r')
            lines = f.readlines()
            mutations = []
            for i in range(0, len(lines)):
                if not lines[i].startswith("#") and lines[i] != "\n":
                    mutations.append(i)
            func_chunk = []
            for index in range(0, len(mutations)):
                # fetching function if there is only one mutation
                if index == 0:
                    func_chunk = lines[0:mutations[index] + 1]
                    # print(function)
                    dataFrame = extract_metadata(inp_file=file,
                                                 chunk=func_chunk,
                                                 df=dataFrame)
                else:
                    func_chunk = lines[mutations[index - 1] +
                                       2:mutations[index] + 1]
                    dataFrame = extract_metadata(inp_file=file,
                                                 chunk=func_chunk,
                                                 df=dataFrame)

    # fill in first three columns
    reference_accession = args.accession
    dataFrame['reference accession'] = reference_accession
    dataFrame['reference database name'] = 'RefSeq'
    if reference_accession=='NC_045512.2':
        dataFrame['organism'] = 'Severe acute respiratory syndrome coronavirus 2'
    elif reference_accession=='NC_063383.1':
        dataFrame['organism'] = 'Monkeypox virus'
    else:
        dataFrame['organism'] = 'unknown'

    dataFrame['mutation functional annotation resource'] = 'Pokay'

    # add HGVS mutations names and nucleotide positions from mutation index
    mutation_index = pd.read_csv(args.mutation_index, sep='\t')
    # merge on "mutation" and "gene" in index ("original mutation description" and "gene symbol" in functional annotation file)
    mutation_index = mutation_index.rename(columns={"mutation": "original mutation description", 'pos':'nucleotide position',
                                                    'hgvs_aa_mutation':'amino acid mutation','hgvs_nt_mutation':'nucleotide mutation',
                                                    'gene':'protein symbol','hgvs_alias':'amino acid mutation alias'})
    # remove doubled columns from dataFrame
    index_cols_to_use = ['nucleotide position','nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias']
    dataFrame = dataFrame.drop(columns=index_cols_to_use)
    merged_dataFrame = pd.merge(dataFrame, mutation_index, on=['original mutation description', 'protein symbol'], how='left') #, 'alias'
    merged_dataFrame = merged_dataFrame[dataFrame_cols]

    # tidying
    merged_dataFrame = merged_dataFrame[merged_dataFrame["nucleotide position"].notna()] # names that aren't matched in the index are dropped here
    merged_dataFrame["nucleotide position"] = merged_dataFrame["nucleotide position"].astype(int)
    merged_dataFrame["author"] = merged_dataFrame["author"].str.strip()
    merged_dataFrame["DOI"] = merged_dataFrame["DOI"].str.strip()
    merged_dataFrame["URL"] = merged_dataFrame["URL"].str.strip()
    merged_dataFrame["peer review status"] = ~merged_dataFrame["peer review status"].str.contains("Preprint")

    write_tsv(dframe=merged_dataFrame)
