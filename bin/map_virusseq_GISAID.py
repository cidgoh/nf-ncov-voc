#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script is a utility script made for mapping VirusSeq Data
portal dataset to GISAID metadata in order to fetch lineage
information that is not updated at data portal. If pangolin is used,
this script is not required in the nf-ncov-voc workflow.

"""

import argparse
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Map VirusSeq data to GISAID Metadata for '
                    'pangolin')
    parser.add_argument('--virusseq', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--gisaid', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--output', type=str, default=None,
                        help='Metadata file (.tsv) format')

    return parser.parse_args()


def write_metadata(dataframe):
    dataframe.to_csv(args.output,
                     sep="\t",
                     quoting=csv.QUOTE_NONE,
                     index=False, header=True)


if __name__ == '__main__':
    args = parse_args()

    virus_seq_df = pd.read_csv(args.virusseq, sep="\t",
                               low_memory=False)
    print("Number of sequences in VirusSeq Data Portal: ",
          len(virus_seq_df))

    print("Sequences in VirusSeq Data Portal with GISAID accessions: ",
          len(virus_seq_df.loc[
                  virus_seq_df['GISAID accession'].notna()]))

    gisaid_df = pd.read_csv(args.gisaid, sep="\t", low_memory=False)
    gisaid_df = gisaid_df.loc[gisaid_df['country'] == 'Canada']
    print("Number of Canadian sequences in GISAID: ", len(gisaid_df))

    df = pd.merge(gisaid_df, virus_seq_df, how='right',
                  left_on='gisaid_epi_isl',
                  right_on='GISAID accession', indicator=True)

    df = df.loc[df['_merge'] == 'both']
    df = df[['strain', 'virus', 'gisaid_epi_isl',
             'genbank_accession', 'pango_lineage', 'date', 'region',
             'country', 'division', 'location', 'submitting_lab',
             'date_submitted', 'fasta header name']]
    df = df.drop(['strain'], axis=1)
    df = df.rename(columns={'fasta header name': 'strain'})
    write_metadata(dataframe=df)
