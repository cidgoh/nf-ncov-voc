# temporary code for making JSON files out of ROBOT tables
# ROBOT table: 
# example runs:
# python robot2json.py --robot "/home/madeline/Downloads/VIRUS-MVP GENEPIO_ROBOT Tables - sc2_gene_protein_data.tsv" --ontology_genes '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_genes.json' --ontology_proteins '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_proteins.json'
# python robot2json.py --robot ~/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_063383.1/NC_063383.1_ROBOT.tsv --ontology_genes '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_063383.1/NC_063383.1_ontology_genes.json' --ontology_proteins '/home/madeline/Desktop/git_temp/nf-ncov-voc/assets/virus_ontologyTerms/NC_063383.1/NC_063383.1_ontology_proteins.json'
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
    sars = pd.read_csv(robot_path, header=0, sep='\t', usecols=['Dbxref', 'Ontology ID', 'parent class', 'label'])

    # process genes
    gene_names = sars[sars['parent class'].str.contains("gene name")].copy(deep=True)
    gene_names['gene_name'] = gene_names['label']
    gene_names['gene_name_id'] = gene_names['Ontology ID']
    gene_names['Dbxref'] = gene_names['Dbxref'].str.replace("GeneID", "NCBIGene") #'JSON_match'
    gene_names = gene_names[['Dbxref', 'gene_name', 'gene_name_id']]
    gene_names.to_csv("GENE_NAMES_HAVOC.tsv", sep='\t', header=True, index=False)

    gene_symbols = sars[sars['parent class'].str.contains("gene symbol")].copy(deep=True)
    gene_symbols['gene_symbol'] = gene_symbols['label']
    gene_symbols['gene_symbol_id'] = gene_symbols['Ontology ID']
    gene_symbols = gene_symbols[['Dbxref', 'gene_symbol', 'gene_symbol_id']]
    gene_symbols['Dbxref'] = gene_symbols['Dbxref'].str.replace("GeneID", "NCBIGene") #'JSON_match'
    gene_symbols.to_csv("GENE_SYMBOLS_HAVOC.tsv", sep='\t', header=True, index=False)

    genes = pd.merge(gene_names, gene_symbols, on='Dbxref', how='left')

    # add strand orientation
    genes['strand_orientation'] = "sense strand orientation"
    genes['strand_orientation_id'] = "GENEPIO:0700002"
    
    # add gene orientation
    genes['gene_orientation'] = "forward gene orientation"
    genes['gene_orientation_id'] = "GENEPIO:0700004"
    reverse_gene_ids = ["NCBIGene:72551595", "NCBIGene:72551594", "NCBIGene:72551593", "NCBIGene:72551592"]
    #print(genes[genes['gene_symbol'].duplicated(keep=False)]) #show which genes are duplicated
    genes.loc[genes["Dbxref"].isin(reverse_gene_ids), 'gene_orientation'] = "reverse gene orientation"
    genes.loc[genes["Dbxref"].isin(reverse_gene_ids), 'gene_orientation_id'] = "GENEPIO:0700005"

    genes.to_csv("GENES_HAVOC.tsv", sep='\t', header=True, index=False)

    # Define a custom function to create a nested structure
    # ncbigene and the gene name are all the same, this is causing HAVOC
    def custom_nested_structure_genes(row):
        return {row['Dbxref']: {'type': 'gene', 'Dbxref': row['Dbxref'], 
                                'gene_name': {'label': row['gene_name'], 'uri': row['gene_name_id']}, 
                                'gene_symbol': {'label': row['gene_symbol'], 'uri': row['gene_symbol_id']}, 
                                'strand_orientation': {'label': row['strand_orientation'], 'uri': row['strand_orientation_id']},
                                'gene_orientation': {'label': row['gene_orientation'], 'uri': row['gene_orientation_id']}
                                }}

    # Apply the custom function to each row of the DataFrame
    json_data_custom = genes.apply(custom_nested_structure_genes, axis=1).tolist() # saves a list of dicts
    json_data = {k:v for d in json_data_custom for k, v in d.items()} # transform to a single dict

    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(json_data, indent=4)
    #print(json_str)

    # Write the JSON string to a file
    with open(genes_file_path, "w") as genes_json_file:
        genes_json_file.write(json_str)

    print("Gene terms saved as " + genes_file_path)


    # now do proteins

    protein_names = sars[sars['parent class'].str.contains("protein name")].copy(deep=True)
    protein_names['protein_name'] = protein_names['label']
    protein_names['protein_name_id'] = protein_names['Ontology ID']
    protein_names['protein_id'] = protein_names['Dbxref'].str.split(pat=",", expand=True)[0].str.split(":", expand=True)[1] # extract eg. "YP_010377182.1" from "GenBank:YP_010377182.1,GeneID:72551595"
    protein_names = protein_names[['protein_id', 'protein_name', 'protein_name_id']]

    protein_symbols = sars[sars['parent class'].str.contains("protein symbol")].copy(deep=True)
    protein_symbols['protein_symbol'] = protein_symbols['label']
    protein_symbols['protein_symbol_id'] = protein_symbols['Ontology ID']
    protein_symbols['protein_id'] = protein_symbols['Dbxref'].str.split(pat=",", expand=True)[0].str.split(":", expand=True)[1] # extract eg. "YP_010377182.1" from "GenBank:YP_010377182.1,GeneID:72551595"
    protein_symbols = protein_symbols[['protein_id', 'protein_symbol', 'protein_symbol_id']]

    proteins = pd.merge(protein_names, protein_symbols, on='protein_id', how='left')

    # Define a custom function to create a nested structure
    def custom_nested_structure_proteins(row):
        return {row['protein_id']: {'type': 'CDS', 'protein_id': row['protein_id'],
                                    'protein_name': {'label': row['protein_name'], 'uri': row['protein_name_id']}, 
                                    'protein_symbol': {'label': row['protein_symbol'], 'uri': row['protein_symbol_id']}}}
    
    # Apply the custom function to each row of the DataFrame
    json_data_custom = proteins.apply(custom_nested_structure_proteins, axis=1).tolist() # list of dicts
    json_data = {k: v for d in json_data_custom for k, v in d.items()} # transform to a single dict

    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(json_data, indent=4)

    # Write the JSON string to a file
    with open(proteins_file_path, "w") as protein_json_file:
        protein_json_file.write(json_str)

    print("Protein terms saved as " + proteins_file_path)