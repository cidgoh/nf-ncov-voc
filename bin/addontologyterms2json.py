# This script adds the ontology terms in assets/virus_ontologyTerms to the master JSON file. The added keys are 'gene_name', 'gene_symbol', 'protein_name', and 'protein_symbol'.
#Example run:
#python addontologyterms2json.py --gene_positions ~/Downloads/NC_045512.2.json --ontology_genes ../assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_genes.json --ontology_proteins ../assets/virus_ontologyTerms/NC_045512.2/NC_045512.2_ontology_proteins.json --savefile NC_045512.2_ontologytermsadded.json

import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a annotated VCF file to a GVF '
                    'file with functional annotation')
    parser.add_argument('--gene_positions', type=str, default=None, required=True,
                        help='Path to the JSON file of gene positions to add ontology terms to')
    parser.add_argument('--ontology_genes', type=str, required=True,
                        help='Path to the ontology gene terms file')
    parser.add_argument('--ontology_proteins', type=str, required=True,
                        help='Path to the ontology protein terms file')
    parser.add_argument('--savefile', type=str, required=True,
                        help='Filename for the output JSON file')
    return parser.parse_args()

if __name__ == '__main__':

    # define args
    args = parse_args()
    gene_positions = args.gene_positions
    ontology_genes = args.ontology_genes
    ontology_proteins = args.ontology_proteins
    savefile = args.savefile
    
    # Open the gene_positions JSON file
    with open(gene_positions) as fp:
        GENE_PROTEIN_POSITIONS_DICT = json.load(fp)

    # Open the ontology_genes JSON file
    with open(ontology_genes) as fp:
        ONTOLOGY_GENES_DICT = json.load(fp)

    # Open the ontology_proteins JSON file
    with open(ontology_proteins) as fp:
        ONTOLOGY_PROTEINS_DICT = json.load(fp)

    # loop through all CDS and gene regions in gene_positions
    for entry in GENE_PROTEIN_POSITIONS_DICT.keys():

        if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="CDS":
            protein_id = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_id"]
            # find matching CDS entry in ontology_proteins JSON
            for protein_entry in ONTOLOGY_PROTEINS_DICT.keys():
                if ONTOLOGY_PROTEINS_DICT[protein_entry]["protein_id"] == protein_id:
                    protein_name = ONTOLOGY_PROTEINS_DICT[protein_entry]["protein_name"]
                    protein_symbol = ONTOLOGY_PROTEINS_DICT[protein_entry]["protein_symbol"]
            GENE_PROTEIN_POSITIONS_DICT[entry]["protein_name"] = protein_name
            GENE_PROTEIN_POSITIONS_DICT[entry]["protein_symbol"] = protein_symbol

        if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="gene":
            gene_ids = GENE_PROTEIN_POSITIONS_DICT[entry]["Dbxref"].replace("GeneID", "NCBIGene")
            # find matching CDS entry in ontology_proteins JSON
            for gene_entry in ONTOLOGY_GENES_DICT.keys():
                if ONTOLOGY_GENES_DICT[gene_entry]["Dbxref"] == gene_ids:
                    gene_name = ONTOLOGY_GENES_DICT[gene_entry]["gene_name"]
                    gene_symbol = ONTOLOGY_GENES_DICT[gene_entry]["gene_symbol"]
            GENE_PROTEIN_POSITIONS_DICT[entry]["gene_name"] = gene_name
            GENE_PROTEIN_POSITIONS_DICT[entry]["gene_symbol"] = gene_symbol


    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(GENE_PROTEIN_POSITIONS_DICT, indent=4)

    # Write the JSON string to a file
    with open(savefile, "w") as json_file:
        json_file.write(json_str)