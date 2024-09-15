# temporary code for making JSON files out of ROBOT tables
# ROBOT table: 
# example run:
# python robot2json.py --robot "/home/madeline/Downloads/VIRUS-MVP GENEPIO_ROBOT Tables - sc2_gene_protein_data.tsv" --ontology_genes '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_genes.json' --ontology_proteins '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_proteins.json'
# remove orf1a gene name and symbol from TSV as they aren't used in the existing JSON

import pandas as pd
import json
from json import loads, dumps
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a annotated VCF file to a GVF '
                    'file with functional annotation')
    parser.add_argument('--robot', type=str, default=None, required=True,
                        help='Path to the modified ROBOT file of gene and protein names and symbols')
    parser.add_argument('--ontology_genes', type=str, required=True,
                        help='Path to save the ontology gene terms JSON to')
    parser.add_argument('--ontology_proteins', type=str, required=True,
                        help='Path to save the ontology protein terms JSON to')
    return parser.parse_args()

if __name__ == '__main__':

    # define args
    args = parse_args()
    robot_path = args.robot
    genes_file_path = args.ontology_genes
    proteins_file_path = args.ontology_proteins

    # load ROBOT file
    sars = pd.read_csv(robot_path, header=0, sep='\t', usecols=['JSON_match', 'Ontology ID', 'parent class', 'label', 'is about', 'is about label'])

    # process genes
    gene_names = sars[sars['parent class'].str.contains("gene name")]
    gene_names['gene_name'] = gene_names['label']
    gene_names['gene_name_id'] = gene_names['Ontology ID']
    gene_names['Dbxref'] = gene_names['JSON_match'].str.replace("GeneID", "NCBIGene")
    gene_names = gene_names[['Dbxref', 'gene_name', 'gene_name_id', 'is about label']]

    gene_symbols = sars[sars['parent class'].str.contains("gene symbol")]
    gene_symbols['gene_symbol'] = gene_symbols['label']
    gene_symbols['gene_symbol_id'] = gene_symbols['Ontology ID']
    gene_symbols = gene_symbols[['gene_symbol', 'gene_symbol_id', 'is about label']]

    genes = pd.merge(gene_names, gene_symbols, on='is about label', how='left')

    # Define a custom function to create a nested structure
    def custom_nested_structure_genes(row):
        return {row['Dbxref']: {'type': 'gene', 'Dbxref': row['Dbxref'], 'gene_name': {'label': row['gene_name'], 'uri': row['gene_name_id']}, 'gene_symbol': {'label': row['gene_symbol'], 'uri': row['gene_symbol_id']}}}

    # Apply the custom function to each row of the DataFrame
    json_data_custom = genes.apply(custom_nested_structure_genes, axis=1).tolist() # saves a list of dicts
    json_data = {k: v for d in json_data_custom for k, v in d.items()} # transform to a single dict

    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(json_data, indent=4)

    # Write the JSON string to a file
    with open(genes_file_path, "w") as genes_json_file:
        genes_json_file.write(json_str)


    # now do proteins

    protein_names = sars[sars['parent class'].str.contains("protein name")]
    protein_names['protein_name'] = protein_names['label']
    protein_names['protein_name_id'] = protein_names['Ontology ID']
    protein_names['protein_id'] = protein_names['JSON_match']
    protein_names = protein_names[['protein_id', 'protein_name', 'protein_name_id', 'is about label']]

    protein_symbols = sars[sars['parent class'].str.contains("protein symbol")]
    protein_symbols['protein_symbol'] = protein_symbols['label']
    protein_symbols['protein_symbol_id'] = protein_symbols['Ontology ID']
    protein_symbols = protein_symbols[['protein_symbol', 'protein_symbol_id', 'is about label']]

    proteins = pd.merge(protein_names, protein_symbols, on='is about label', how='left')

    # Define a custom function to create a nested structure
    def custom_nested_structure_proteins(row):
        return {row['protein_id']: {'type': 'CDS', 'protein_id': row['protein_id'], 'protein_name': {'label': row['protein_name'], 'uri': row['protein_name_id']}, 'protein_symbol': {'label': row['protein_symbol'], 'uri': row['protein_symbol_id']}}}
    
    # Apply the custom function to each row of the DataFrame
    json_data_custom = proteins.apply(custom_nested_structure_proteins, axis=1).tolist() # list of dicts
    json_data = {k: v for d in json_data_custom for k, v in d.items()} # transform to a single dict

    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(json_data, indent=4)

    # Write the JSON string to a file
    with open(proteins_file_path, "w") as protein_json_file:
        protein_json_file.write(json_str)