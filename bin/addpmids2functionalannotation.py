#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: madeline

This script adds PMIDs to the Pokay functional annotation
key (tsv) file that can be used by nf-ncov-voc workflow for
annotation of variant called files.

"""

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script adds PMIDs to the functional annotations file')
    parser.add_argument('--functional_annotations', type=str,
                        required=True, help='TSV file of functional annotations')
    parser.add_argument('--pmid_tsv', type=str, required=True,
                        help='TSV file containing PMIDs and DOIs')
    parser.add_argument('--outputfile', type=str, required=True,
                        help='output file (.TSV) format')
    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()

    # load functional annotations
    functional_annotation_df = pd.read_csv(args.functional_annotations, sep='\t', header=0)
    functional_annotation_columns = functional_annotation_df.columns

    # load PMIDs file
    pmids_df = pd.read_csv(args.pmid_tsv, sep=',', header=0)
    pmids_df = pmids_df[pmids_df["PMID"].notna()]
    pmids_df = pmids_df.drop_duplicates(subset='DOI', keep='first')
    pmids_df = pmids_df[['PMID', 'DOI']]
    pmids_df['DOI'] = "doi:" + pmids_df['DOI']

    #merge them
    functional_annotation_df = functional_annotation_df.drop(columns='PMID')
    functional_annotation_df = pd.merge(functional_annotation_df, pmids_df, on='DOI', how='left')
    functional_annotation_df["PMID"] = functional_annotation_df["PMID"].astype('Int64')

    # rearrange columns
    functional_annotation_df = functional_annotation_df[functional_annotation_columns]

    # make sure years are integers
    #functional_annotation_df['publication year'] = functional_annotation_df['publication year'].astype(str).replace('.0', 'CAT', regex=False)
    functional_annotation_df["publication year"] = functional_annotation_df["publication year"].astype('Int64')

    # save to TSV
    functional_annotation_df.to_csv(args.outputfile, sep='\t', header=True, index=False)