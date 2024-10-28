#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: anoosha & madeline

This script converts Pokay .txt files into a functional annotation
key (tsv) files that can be used by nf-ncov-voc workflow for
annotation of variant called files.

"""

import json
import os
import argparse
from pathlib import Path
import pandas as pd
import csv
from functions import map_pos_to_gene_protein, unnest_multi


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script produces a TSV file from TXT files '
                    'in POKAY  '
                    'https://github.com/nodrogluap/pokay/tree/master'
                    '/data for annotating SARS-COV-2 mutations')
    parser.add_argument('--inputdir', type=str, required=True,
                        help='directory path for input files')
    parser.add_argument('--accession', type=str, required=True,
                        help='versioned reference accession from RefSeq')
    parser.add_argument('--mutation_index', type=str, default='n/a',
                        help='index of all mutations (.TSV)')
    parser.add_argument('--gene_positions', type=str,
                        default=None, required=True,
                        help='gene positions in JSON format')
    parser.add_argument('--outputfile', type=str, required=True,
                        help='output file (.TSV) format')
    parser.add_argument('--save_dois', type=str, default=None,
                        help='output file (.TXT) format')
    return parser.parse_args()


def combination_mutation(c_mutations, c_mutation):
    comb_mutation = []
    for y in range(len(c_mutations)):
        if c_mutations[y] != c_mutation:
            comb_mutation.append(c_mutations[y])
    return comb_mutation


def data_cleanup(dframe):
    dframe['variant functional effect description'] = dframe[
        'variant functional effect description'].apply(
        lambda x: ','.join(map(str, x)))

    dframe['variant functional effect description'] = dframe[
        'variant functional effect description'].str.replace(',#', '')
    dframe['variant functional effect description'] = dframe[
        'variant functional effect description'].str.replace('#', '')
    dframe['variant functional effect description'] = dframe[
        'variant functional effect description'].str.strip()
    
    #dframe['comb_mutation'] = dframe['comb_mutation'].apply(
    #    lambda x: x[1:-1])

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
    dframe['peer review status'] = dframe['URL'].astype(str).str.extract('(?<=\[)(.*)')
    dframe['URL'] = dframe['URL'].astype(str).str.extract('^(.*?)\[') #.str.strip()
    dframe['DOI'] = dframe['URL'].astype(str).str.extract('(?<=doi.org/)(.*)') #.str.strip()

    return dframe


def extract_metadata(inp_file, chunk, df):
    mutation_name = chunk[-1].strip()
    '''
    ## MRI: for now, removing 'comb_mutation' column and having 'mutation_name' comma-separated
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
    '''
    mutation_name = mutation_name.replace(';', ',')
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
                                                          'variant functional effect description'])

        df_list = [mutation_name, gene_name,
                   function_category] #, comb_mutation, heterozygosity]
        # print(df_list)
        df1 = pd.DataFrame(
            columns=['original mutation description', 'pokay_id', 'measured variant functional effect'])
                     #'comb_mutation', 'heterozygosity'])
        df1.loc[len(df1)] = df_list

        df_func['original mutation description'] = df1['original mutation description'].iloc[0]
        df_func['pokay_id'] = df1['pokay_id'].iloc[0]
        df_func['measured variant functional effect'] = df1['measured variant functional effect'].iloc[0]
        #df_func['comb_mutation'] = str(df1['comb_mutation'].iloc[0])
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

    if args.accession=='NC_063383.1': # add gene orientation and strand orientation for MPOX
        dataFrame_cols = ['organism', 'reference accession', 'reference database name', 'nucleotide position',
'original mutation description', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias',
'gene name', 'gene symbol', 'gene orientation', 'strand orientation', 'protein name', 'protein symbol', 'measured variant functional effect', 'inferred variant functional effect', 'viral life cycle functional effect',
'variant functional effect description', 'CVX code', 'DrugBank Accession Number', 'Antibody Registry ID', 'author', 'publication year', 'URL', 'DOI', 'PMID',
'peer review status', 'curator', 'mutation functional annotation resource']
    
    else:
        dataFrame_cols = ['organism', 'reference accession', 'reference database name', 'nucleotide position',
'original mutation description', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias',
'gene name', 'gene symbol', 'protein name', 'protein symbol', 'measured variant functional effect', 'inferred variant functional effect', 'viral life cycle functional effect',
'variant functional effect description', 'CVX code', 'DrugBank Accession Number', 'Antibody Registry ID', 'author', 'publication year', 'URL', 'DOI', 'PMID',
'peer review status', 'curator', 'mutation functional annotation resource']
    
    dataFrame = pd.DataFrame(columns=dataFrame_cols)

    for file in os.listdir(path):
        if file.endswith(".txt") and "_" in file:
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

    # strip whitespace
    dataFrame["author"] = dataFrame["author"].str.strip()
    dataFrame["DOI"] = dataFrame["DOI"].str.strip()
    dataFrame["URL"] = dataFrame["URL"].str.strip()

    # add BioRegistry prefix for doi
    dataFrame["DOI"] = "doi:" + dataFrame["DOI"]

    # add temporary curator name
    dataFrame["curator"] = "Paul Gordon" ## for now: need to go through

    # clean up 'peer review status'
    dataFrame.loc[dataFrame['peer review status'].isna(), 'peer review status'] = '' ##up for debate
    dataFrame.loc[dataFrame['peer review status'].str.contains("Journal"), 'peer review status'] = 'true'
    dataFrame.loc[dataFrame['peer review status'].str.contains("Preprint"), 'peer review status'] = 'false'
    dataFrame.loc[dataFrame['peer review status'].str.contains("Grey literature"), 'peer review status'] = 'false' ##up for debate

    # map 'pokay_id' to 'protein_symbol' and 'mat_pep' based on JSON
    pokay_id_set = set(dataFrame['pokay_id'].tolist())
    #pokay_id_set = {'RdRp', 'nsp4', 'M', 'ORF8', '3CL', 'nsp1', 'E', 'S', 'nsp8', 'ORF6', 'N', 'PLpro', 'nsp7', 'ORF3a', 'nsp9', 'nsp6', 'nsp2'}
    dataFrame['gene name'] = ''
    dataFrame['gene symbol'] = ''
    dataFrame['protein name'] = ''
    dataFrame['protein symbol'] = ''
    dataFrame['mat_pep'] = ''

    # Reading the gene & protein coordinates of SARS-CoV-2 genome
    with open(args.gene_positions) as fp:

        GENE_PROTEIN_POSITIONS_DICT = json.load(fp)

        for pokay_id in pokay_id_set:
            for entry in GENE_PROTEIN_POSITIONS_DICT.keys():

                # If a gene record matches pokay_id, get 'gene_name' and 'gene_symbol'
                if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="gene" and ("gene" in GENE_PROTEIN_POSITIONS_DICT[entry].keys()) and GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]==pokay_id:
                    # extract protein names and symbols (ontology) from JSON entry
                    gene_name = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_name"]["label"]
                    gene_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_symbol"]["label"]
                    gene_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_orientation"]["label"]
                    strand_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["strand_orientation"]["label"]
                    # add gene names and symbols to dataframe
                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene name"] = gene_name
                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene symbol"] = gene_symbol
                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene"] = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
                    if args.accession=='NC_063383.1': # add gene orientation and strand orientation for MPOX
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene orientation"] = gene_orientation
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "strand orientation"] = strand_orientation

                # If a CDS record matches pokay_id in 'gene' or 'protein_alias', get 'protein_name' and 'protein_symbol'
                if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="CDS" and ("gene" in GENE_PROTEIN_POSITIONS_DICT[entry].keys()):
                    if (GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]==pokay_id):
                        # extract protein names and symbols (ontology) from JSON entry
                        protein_name = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_name"]["label"]
                        protein_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_symbol"]["label"]
                        # add protein names and symbols to dataframe
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein name"] = protein_name
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein symbol"] = protein_symbol
                    elif (pokay_id in GENE_PROTEIN_POSITIONS_DICT[entry]["protein_alias"]):
                        gene = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
                        # extract protein names and symbols (ontology) from JSON entry
                        protein_name = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_name"]["label"]
                        protein_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_symbol"]["label"]
                        # add protein names and symbols to dataframe
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein name"] = protein_name
                        dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein symbol"] = protein_symbol
                        # get 'gene_name' and 'gene_symbol'
                        for entry in GENE_PROTEIN_POSITIONS_DICT.keys():
                            if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="gene" and ("gene" in GENE_PROTEIN_POSITIONS_DICT[entry].keys()) and GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]==gene:
                                # extract protein names and symbols (ontology) from JSON entry
                                gene_name = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_name"]["label"]
                                gene_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_symbol"]["label"]
                                gene_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_orientation"]["label"]
                                strand_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["strand_orientation"]["label"]
                                # add gene names and symbols to dataframe
                                dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene name"] = gene_name
                                dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene symbol"] = gene_symbol
                                dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene"] = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
                                if args.accession=='NC_063383.1': # add gene orientation and strand orientation for MPOX
                                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene orientation"] = gene_orientation
                                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "strand orientation"] = strand_orientation


                # If a mature peptide record matches pokay_id in the 'protein_alias' list, get name of parent, and from there get names and symbols
                if (GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="mature_protein_region_of_CDS") \
                    and (pokay_id in GENE_PROTEIN_POSITIONS_DICT[entry]["protein_alias"]):
                    # extract pokay_id, mat_pep, and parent id from JSON entry
                    mat_pep = entry
                    parent = GENE_PROTEIN_POSITIONS_DICT[entry]["Parent"]
                    # get parent protein name and symbol from corresponding CDS entry, matching "Parent" to "ID"
                    for entry in GENE_PROTEIN_POSITIONS_DICT.keys():
                        if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="CDS" and (GENE_PROTEIN_POSITIONS_DICT[entry]["ID"]==parent):
                            parent_protein_name = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_name"]["label"]
                            parent_protein_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_symbol"]["label"]
                            dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein name"] = parent_protein_name
                            dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "protein symbol"] = parent_protein_symbol
                            parent_gene = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
                            # get parent gene name and symbol from corresponding parent gene entry
                            for entry in GENE_PROTEIN_POSITIONS_DICT.keys():
                                if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="gene" and (GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]==parent_gene):
                                    # extract gene names and symbols (ontology) from JSON entry
                                    parent_gene_name = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_name"]["label"]
                                    parent_gene_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_symbol"]["label"]
                                    # add gene names and symbols to dataframe
                                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene name"] = parent_gene_name
                                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene symbol"] = parent_gene_symbol
                                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "gene"] = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
                    # fill in mat_pep
                    dataFrame.loc[dataFrame["pokay_id"]==pokay_id, "mat_pep"] = mat_pep

#    dataFrame[['pokay_id', 'gene name', 'gene symbol', 'protein name', 'protein symbol', 'mat_pep', 'gene']].drop_duplicates().to_csv('../test_data/pokay_id.tsv', sep='|', index=False)
#    print("saved as '../test_data/pokay_id.tsv'")

    #before merging, unnest the 'original mutation description' column
    dataFrame["original mutation description"] = dataFrame["original mutation description"].str.split(',')
    dataFrame = unnest_multi(dataFrame, ["original mutation description"], reset_index=False)
    dataFrame['index1'] = dataFrame.index

    # if mutation index: add HGVS mutations names, nucleotide positions from it
    if args.mutation_index != 'n/a':
        mutation_index = pd.read_csv(args.mutation_index, sep='\t')
        # merge on "mutation" and "protein symbol" in index ("original mutation description" and "protein symbol" in functional annotation file)
        mutation_index = mutation_index.rename(columns={"mutation": "original mutation description", 'pos':'nucleotide position',
                                                        'hgvs_aa_mutation':'amino acid mutation','hgvs_nt_mutation':'nucleotide mutation',
                                                        'hgvs_alias':'amino acid mutation alias', 'gene_symbol':'gene'})
        #print("columns", mutation_index.columns)
        # remove index columns that don't have a nucleotide mutation entry
        initial_length = mutation_index.shape[0]
        mutation_index = mutation_index[mutation_index['nucleotide mutation'].notna()]
        after_length = mutation_index.shape[0]
        print("Removed " + str(initial_length - after_length) + "/" + str(after_length) + " mutation index rows that are missing a nucleotide mutation")
        # remove doubled columns from dataFrame
        index_cols_to_use = ['nucleotide position', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias']
        dataFrame = dataFrame.drop(columns=index_cols_to_use)
        merged_dataFrame = pd.merge(dataFrame, mutation_index, on=['original mutation description', 'gene'], how='left') #, 'mat_pep'
        #dups = mutation_index[mutation_index.duplicated(subset=['nucleotide position', 'original mutation description'], keep=False)]
        #dups = dups.sort_values(by='nucleotide position')
        #dups.to_csv('madeline_testing/dups.tsv', sep='\t', index=False)

        # keep only specified columns, discarding 'pokay_id' here
        merged_dataFrame = merged_dataFrame[['index1'] + dataFrame_cols]

        # drop mutations not found in the index
        merged_dataFrame = merged_dataFrame[merged_dataFrame["nucleotide position"].notna()]

    # allow mutation index to be optional
    if args.mutation_index == 'n/a':
        merged_dataFrame = dataFrame

    # convert all columns to string type and fillna with empty strings
    for column in merged_dataFrame.columns:
        merged_dataFrame[column] = merged_dataFrame[column].astype("string")
        merged_dataFrame[column] = merged_dataFrame[column].fillna('')

    # convert nucleotide positions to int format (still str)
    merged_dataFrame["nucleotide position"] = merged_dataFrame["nucleotide position"].str.replace('.0', '', regex=False)


    # create agg dictionary
    to_list = ['original mutation description', 'nucleotide position', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias']
    list_dict = dict.fromkeys(to_list, ','.join)
    to_first = list(set(dataFrame_cols) - set(to_list))
    first_dict = dict.fromkeys(to_first, 'first')
    agg_dict = list_dict.copy()
    agg_dict.update(first_dict)

    # perform groupby and aggregation
    merged_dataFrame = merged_dataFrame.groupby(by=['index1'], as_index=False).agg(agg_dict)

    #remove strings of commas in mutation name columns
    for column in ['nucleotide position', 'nucleotide mutation', 'amino acid mutation', 'amino acid mutation alias']:
        merged_dataFrame.loc[merged_dataFrame[column].str.contains(',,', regex=False), column] = ''
        merged_dataFrame.loc[merged_dataFrame[column]==',', column] = ''

    # reorder columns and drop 'index1'
    merged_dataFrame = merged_dataFrame[dataFrame_cols]

    # sort by nucleotide position and reindex
    merged_dataFrame = merged_dataFrame.sort_values(by='nucleotide position')
    merged_dataFrame = merged_dataFrame.reindex()

    write_tsv(dframe=merged_dataFrame)

    # Adding PMIDs requires running a shell script, dois2pmcids.sh, with this txt file of DOIs as the input.
    # After that, run addPMIDs2functionalannotation.py.
    ###TO DO: modify dois2pmcids.sh to take the whole TSV as input and add PMIDs directly to the TSV to streamline this
    if args.save_dois != None:
        dois = merged_dataFrame[merged_dataFrame["DOI"]!=''].drop_duplicates(subset='DOI')
        dois["DOI"] = dois["DOI"].str.replace("doi:", "", regex=False)
        dois["DOI"].to_csv(args.save_dois, header=False, index=False)