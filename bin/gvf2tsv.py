#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:47:11 2021

@authors: madeline & zohaib

This script converts GVF files to TSVs, for later conversion to a PDF
case report.

"""
import argparse
import pandas as pd
from functions import separate_attributes


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts GVF files to a TSV report')
    parser.add_argument('--gvf_file', type=str, default=None,
                        help='Paths to GVF files to process')
    parser.add_argument('--outtsv', type=str,
                        default="surveillance_report.tsv",
                        help='Filepath for finished .tsv')

    return parser.parse_args()


def gvf2tsv(gvf):
    # read in gvf
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']

    df = pd.read_csv(gvf, sep='\t', names=gvf_columns)
    # remove pragmas and original header
    df = df[~df['#seqid'].str.contains("#")]
    # restart index from 0
    df = df.reset_index(drop=True)

    # expand #attributes column into multiple columns for each attribute,
    # keeping the original #attributes column
    df = separate_attributes(df)
    # change all labels to lowercase
    df.columns = [x.lower() for x in df.columns]
    # drop original #attributes column, and also #source
    df = df.drop(labels=['#source', '#attributes'], axis=1)
    
    # remove '#' from column names
    df.columns = df.columns.str.replace("#", "")

    # drop unwanted columns
    df = df.drop(labels=['seqid', 'type', 'end', 'strand', 'score', 'phase', 'id'], axis=1)

    # rename 'dp' column to 'sequence_depth', make 'viral_lineage'
    # plural
    df = df.rename(columns={'sample_size': 'obs_sample_size',
                            'viral_lineage': 'viral_lineages',
                            'URL': 'citation_url'})

    return df


if __name__ == '__main__':

    args = parse_args()
    outfile = args.outtsv
    gvf_file = args.gvf_file

    # convert gvf file to tsv
    print("")
    print("Processing:")
    print(gvf_file)

    gvf_df = gvf2tsv(gvf=gvf_file)
    gvf_df.to_csv(outfile, sep='\t', index=False)

    print("Processing complete.")
    print("Surveillance report saved as: " + outfile)
