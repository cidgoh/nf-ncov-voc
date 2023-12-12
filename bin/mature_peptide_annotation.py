#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script uses a variant called vcf file and annotates the
mutations with mature peptides using the SARS-CoV-2 genome features
from NCBI.

"""

import argparse
from cyvcf2 import VCF, Writer
import json

def parse_args():
    """
    Parses command line arguments for the script.

    Returns:
        argparse.Namespace: An object containing the parsed arguments.
    """
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Adds mature peptide annotation to VCF files ')
    
    # Add command line arguments to the parser
    parser.add_argument('--vcf_file', type=str, default=None,
                        help='Variant calling output file in VCF format')
    parser.add_argument('--output_vcf', type=str, default=None,
                        help='Output VCF file')
    parser.add_argument('--json', type=str, default=None,
                        help='viral JSON file')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    vcf_file = args.vcf_file
    output_vcf = args.output_vcf
    json_file = args.json

    # The above code assigns the values of the parsed arguments to variables.

    if vcf_file is None or output_vcf is None or json_file is None:
        raise ValueError("Please provide all the required arguments")

    # The above code checks if all the required arguments are provided, and raises an error if any of them is missing.

    
    with open(json_file, 'r') as f:
        gene_protein = json.load(f)
    
    data_vcf = VCF(vcf_file)
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep', 'Description': 'Mature Peptide ID',
         'Type': 'String', 'Number': '.'})
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep_desc',
         'Description': 'Mature Peptide Description',
         'Type': 'String', 'Number': '.'})
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep_acc',
         'Description': 'Mature Peptide Accession Number',
         'Type': 'String', 'Number': '.'})
    
    
    w = Writer(output_vcf, data_vcf)
    for record in data_vcf:
        record.INFO["mat_pep"] = "n/a"
        record.INFO["mat_pep_desc"] = "n/a"
        record.INFO["mat_pep_acc"] = "n/a"
        for key in gene_protein:
            if gene_protein[key]["type"] == "mature_protein_region_of_CDS":
                if int(gene_protein[key]["start"]) <= int(record.POS) <= int(gene_protein[key]["end"]):
                    record.INFO["mat_pep"] = str("".join(gene_protein[key]["protein_alias"])).replace(";", ",")
                    record.INFO["mat_pep_desc"] = str("".join(gene_protein[key]["Note"])).replace(";", ",")
                    record.INFO["mat_pep_acc"] = str("".join(gene_protein[key]["ID"])).replace(";", ",")
        w.write_record(record)    
    w.close()
    data_vcf.close()