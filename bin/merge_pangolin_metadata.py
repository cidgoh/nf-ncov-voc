#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script merges Pangolin report (assigned lineages) with the
metadata file which allows data extraction and filtering based on
lineage information in nf-ncov-voc workflow.

"""

import argparse
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merges pangolin output report and metadata file '
                    'using isolate as key')
    parser.add_argument('--metadata', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--pangolin', type=str, default=None,
                        help='Pangolin report (.csv) format')
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

    metadata_df = pd.read_csv(args.metadata, sep="\t")
    pangolin_df = pd.read_csv(args.pangolin)

    metadata_df['fasta header name'] = metadata_df['fasta header name'].str.strip()
    pangolin_df['taxon'] = pangolin_df['taxon'].str.strip()

    merged_df = pd.merge(metadata_df, pangolin_df, left_on='fasta header name',
                         right_on='taxon')
    merged_df = merged_df.rename(columns={"lineage": "pango_lineage",
                                          "fasta header name": "strain"})

    write_metadata(dataframe=merged_df)
