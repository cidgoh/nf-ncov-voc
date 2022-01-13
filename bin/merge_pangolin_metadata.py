#!/usr/bin/env python3

import argparse
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Map VirusSeq data to GISAID Metadata for '
                    'pangolin')
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

    merged_df = pd.merge(metadata_df, pangolin_df, left_on='isolate',
                         right_on='taxon')
    merged_df = merged_df.rename(columns={"lineage": "pango_lineage",
                                          "isolate": "strain"})

    write_metadata(dataframe=merged_df)
