#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script flags mutations for the problematics sites in the
nf-ncov-voc workflow curated by
EMBL-EBI at https://github.com/W-L/ProblematicSites_SARS-CoV2/
"""

import argparse
import pandas as pd
from cyvcf2 import VCF, Writer


def parse_args():
    parser = argparse.ArgumentParser(
        description='Adds mutation tag to problematic sites based on '
                    'https://github.com/W-L/ProblematicSites_SARS'
                    '-CoV2#human-friendly-version-of-the-vcf-file')
    parser.add_argument('--vcffile', type=str, default=None,
                        help='Variant calling output file in VCF '
                             'format')
    parser.add_argument('--filter_vcf', type=str, default=None,
                        help='Problematic sites in SARS-CoV-2 genomes '
                             'in VCF file')
    parser.add_argument('--output_vcf', type=str, default=None,
                        help='Output VCF file')
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_args()
    
    vcf_file = args.vcffile
    output_vcf = args.output_vcf
    prob_vcf_file = args.filter_vcf

    data_vcf = VCF(vcf_file)
    prob_vcf = VCF(prob_vcf_file)

    # Reading the VCF file and adding 2 more attributes into INFO header
    
    data_vcf.add_info_to_header(
        {'ID': 'ps_filter', 'Description': 'mask/caution',
         'Type': 'String', 'Number': '1'})
    data_vcf.add_info_to_header(
        {'ID': 'ps_exc',
         'Description': 'Reasons for mask/caution',
         'Type': 'String', 'Number': '1'})

    # create a new vcf Writer using the input vcf as a template.

    w = Writer(output_vcf, data_vcf)
    filtered_positions = {}
    for prob_record in prob_vcf:
        #if prob_record.FILTER == "mask" or prob_record.FILTER == "caution":
        filtered_positions[prob_record.POS] = [prob_record.ALT, prob_record.FILTER, prob_record.INFO.get("EXC")]
        
    print(filtered_positions)
    
    # Iterate over the records in the input vcf.
    for record in data_vcf:
        record.INFO["ps_filter"] = "n/a"
        record.INFO["ps_exc"] = "n/a"
        if record.POS in filtered_positions and record.ALT in filtered_positions[record.POS][0]:
            record.INFO["ps_filter"] = filtered_positions[record.POS][1]
            record.INFO["ps_exc"] = filtered_positions[record.POS][2]
        w.write_record(record)

    w.close()
    data_vcf.close()
    prob_vcf.close()
