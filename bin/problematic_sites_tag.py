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

    # Reading the VCF file and adding 2 more attributes into INFO header
    data_vcf = VCF(args.vcffile)
    data_vcf.add_info_to_header(
        {'ID': 'ps_filter', 'Description': 'Mask/Caution',
         'Type': 'String', 'Number': '1'})
    data_vcf.add_info_to_header(
        {'ID': 'ps_exc',
         'Description': 'Reasons for mask/caution',
         'Type': 'String', 'Number': '1'})

    # create a new vcf Writer using the input vcf as a template.
    fname = args.output_vcf

    w = Writer(fname, data_vcf)
    prob_vcf_columns = ['CHROM', 'POS', 'ID', 'REF',
                        'ALT', 'QUAL', 'FILTER',
                        'INFO']
    prob_vcf_df = pd.DataFrame(index=range(0, 478),
                               columns=prob_vcf_columns)

    # Make a dataframe from VCF object for indexed searching
    row = 0
    prob_vcf = VCF(args.filter_vcf)
    for v in prob_vcf:
        prob_vcf_df.iloc[row, [0]] = "MN908947.3"
        prob_vcf_df.iloc[row, [1]] = v.POS
        if not v.ALT:
            prob_vcf_df.iloc[row, [2]] = "."
        else:
            prob_vcf_df.iloc[row, [2]] = v.ID
        prob_vcf_df.iloc[row, [3]] = v.REF
        if not v.ALT:
            prob_vcf_df.iloc[row, [4]] = "."
        else:
            prob_vcf_df.iloc[row, [4]] = ",".join(v.ALT)
        if not v.QUAL:
            prob_vcf_df.iloc[row, [5]] = "."
        else:
            prob_vcf_df.iloc[row, [5]] = v.QUAL
        prob_vcf_df.iloc[row, [6]] = v.FILTER
        prob_vcf_df.iloc[row, [7]] = v.INFO.get('EXC')
        row = row + 1

    # Searching records and adding TAGs into INFO column
    for record in data_vcf:
        record.INFO["ps_filter"] = ""
        record.INFO["ps_exc"] = ""
        if (record.POS in prob_vcf_df["POS"].values) & \
                (record.REF in prob_vcf_df["REF"].values) & \
                (record.ALT in prob_vcf_df["ALT"].values):
            record.INFO["ps_filter"] = prob_vcf_df.loc[
                prob_vcf_df['POS']
                == record.POS, 'FILTER'].item()
            record.INFO["ps_exc"] = prob_vcf_df.loc[prob_vcf_df[
                                                        'POS'] ==
                                                    record.POS,
                                                    'INFO'].item()
        w.write_record(record)

    w.close()
    data_vcf.close()
