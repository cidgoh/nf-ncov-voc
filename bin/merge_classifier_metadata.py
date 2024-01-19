#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

This script merges classifier report (assigned lineages) with the
metadata file which allows data extraction and filtering based on
lineage information in nf-ncov-voc workflow.

"""

import argparse
import pandas as pd
import csv


def parse_args():
    """
    Parses the command-line arguments for merging pangolin output report and metadata file.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Merges pangolin output report and metadata file '
                    'using isolate as key')
    parser.add_argument('--metadata', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--classifier_report', type=str, default=None,
                        help='classifier output report (.csv) format')
    parser.add_argument('--classifier', type=str, default=None,
                        help='classifier name (pangolin or nextclade)')                        
    parser.add_argument('--output', type=str, default=None,
                        help='Metadata file (.tsv) format')

    return parser.parse_args()


def write_metadata(dataframe):
    """
    Writes the merged dataframe to a metadata file.

    Args:
        dataframe (pandas.DataFrame): The merged dataframe.
    """
    dataframe.to_csv(args.output, sep="\t", compression='gzip', quoting=csv.QUOTE_NONE, index=False, header=True)


if __name__ == '__main__':
    args = parse_args()

    metadata_df = pd.read_csv(args.metadata, sep="\t", compression='gzip')
    classifier_df = pd.read_csv(args.classifier_report, sep="\t")

    if args.classifier == "pangolin":
        # Strip leading and trailing whitespaces from the columns
        metadata_df['isolate'] = metadata_df['isolate'].str.strip()
        classifier_df['taxon'] = classifier_df['taxon'].str.strip()

        # Merge the metadata and classifier dataframes based on 'fasta_header_name' and 'taxon' columns
        merged_df = pd.merge(metadata_df, classifier_df, left_on='isolate', right_on='taxon')
        # Rename the columns
        merged_df = merged_df.rename(columns={"lineage": "pango_lineage","fasta_header_name": "strain"})
    
    elif args.classifier == "nextclade":
        # Strip leading and trailing whitespaces from the columns
        metadata_df['isolate'] = metadata_df['isolate'].str.strip()
        classifier_df['seqName'] = classifier_df['seqName'].str.strip()

        # Merge the metadata and classifier dataframes based on 'fasta_header_name' and 'taxon' columns
        merged_df = pd.merge(metadata_df, classifier_df, left_on='isolate', right_on='seqName')
        # Rename the columns
        merged_df = merged_df.rename(columns={"lineage": "nextclade_lineage", "fasta_header_name": "strain"})
    else:
        print("Please specify the classifier name (pangolin or nextclade)")
        exit(1)
    
    write_metadata(dataframe=merged_df)
