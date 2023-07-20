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




#reference,accession,genome,genes/ORF1a/name,genes/ORF1a/coordinates/from,genes/ORF1a/coordinates/to,genes/ORF1b/name,genes/ORF1b/coordinates/from,genes/ORF1b/coordinates/to,genes/S/name,genes/S/coordinates/from,genes/S/coordinates/to,genes/ORF3a/name,genes/ORF3a/coordinates/from,genes/ORF3a/coordinates/to,genes/E/name,genes/E/coordinates/from,genes/E/coordinates/to,genes/M/name,genes/M/coordinates/from,genes/M/coordinates/to,genes/ORF6/name,genes/ORF6/coordinates/from,genes/ORF6/coordinates/to,genes/ORF7a/name,genes/ORF7a/coordinates/from,genes/ORF7a/coordinates/to,genes/ORF7b/name,genes/ORF7b/coordinates/from,genes/ORF7b/coordinates/to,genes/ORF8/name,genes/ORF8/coordinates/from,genes/ORF8/coordinates/to,genes/N/name,genes/N/coordinates/from,genes/N/coordinates/to,genes/ORF10/name,genes/ORF10/coordinates/from,genes/ORF10/coordinates/to,proteins/NSP1/name,proteins/NSP1/gene,proteins/NSP1/description,proteins/NSP1/coordinates/from,proteins/NSP1/coordinates/to,proteins/NSP1/g.coordinates/from,proteins/NSP1/g.coordinates/to,proteins/NSP2/name,proteins/NSP2/gene,proteins/NSP2/coordinates/from,proteins/NSP2/coordinates/to,proteins/NSP2/g.coordinates/from,proteins/NSP2/g.coordinates/to,proteins/NSP3/name,proteins/NSP3/gene,proteins/NSP3/coordinates/from,proteins/NSP3/coordinates/to,proteins/NSP3/g.coordinates/from,proteins/NSP3/g.coordinates/to,proteins/NSP4/name,proteins/NSP4/gene,proteins/NSP4/coordinates/from,proteins/NSP4/coordinates/to,proteins/NSP4/g.coordinates/from,proteins/NSP4/g.coordinates/to,proteins/NSP5/name,proteins/NSP5/gene,proteins/NSP5/description,proteins/NSP5/coordinates/from,proteins/NSP5/coordinates/to,proteins/NSP5/g.coordinates/from,proteins/NSP5/g.coordinates/to,proteins/NSP6/name,proteins/NSP6/gene,proteins/NSP6/coordinates/from,proteins/NSP6/coordinates/to,proteins/NSP6/g.coordinates/from,proteins/NSP6/g.coordinates/to,proteins/NSP7/name,proteins/NSP7/gene,proteins/NSP7/coordinates/from,proteins/NSP7/coordinates/to,proteins/NSP7/g.coordinates/from,proteins/NSP7/g.coordinates/to,proteins/NSP8/name,proteins/NSP8/gene,proteins/NSP8/coordinates/from,proteins/NSP8/coordinates/to,proteins/NSP8/g.coordinates/from,proteins/NSP8/g.coordinates/to,proteins/NSP9/name,proteins/NSP9/gene,proteins/NSP9/coordinates/from,proteins/NSP9/coordinates/to,proteins/NSP9/g.coordinates/from,proteins/NSP9/g.coordinates/to,proteins/NSP10/name,proteins/NSP10/gene,proteins/NSP10/coordinates/from,proteins/NSP10/coordinates/to,proteins/NSP10/g.coordinates/from,proteins/NSP10/g.coordinates/to,proteins/NSP11/name,proteins/NSP11/gene,proteins/NSP11/coordinates/from,proteins/NSP11/coordinates/to,proteins/NSP11/g.coordinates/from,proteins/NSP11/g.coordinates/to,proteins/NSP12/name,proteins/NSP12/gene,proteins/NSP12/description,proteins/NSP12/coordinates/from,proteins/NSP12/coordinates/to,proteins/NSP12/g.coordinates/from,proteins/NSP12/g.coordinates/to,proteins/NSP13/name,proteins/NSP13/gene,proteins/NSP13/description,proteins/NSP13/coordinates/from,proteins/NSP13/coordinates/to,proteins/NSP13/g.coordinates/from,proteins/NSP13/g.coordinates/to,proteins/NSP14/name,proteins/NSP14/gene,proteins/NSP14/description,proteins/NSP14/coordinates/from,proteins/NSP14/coordinates/to,proteins/NSP14/g.coordinates/from,proteins/NSP14/g.coordinates/to,proteins/NSP15/name,proteins/NSP15/gene,proteins/NSP15/description,proteins/NSP15/coordinates/from,proteins/NSP15/coordinates/to,proteins/NSP15/g.coordinates/from,proteins/NSP15/g.coordinates/to,proteins/NSP16/name,proteins/NSP16/gene,proteins/NSP16/description,proteins/NSP16/coordinates/from,proteins/NSP16/coordinates/to,proteins/NSP16/g.coordinates/from,proteins/NSP16/g.coordinates/to,proteins/S/name,proteins/S/gene,proteins/S/description,proteins/S/coordinates/from,proteins/S/coordinates/to,proteins/S/g.coordinates/from,proteins/S/g.coordinates/to,proteins/ORF3a/name,proteins/ORF3a/gene,proteins/ORF3a/coordinates/from,proteins/ORF3a/coordinates/to,proteins/ORF3a/g.coordinates/from,proteins/ORF3a/g.coordinates/to,proteins/E/name,proteins/E/gene,proteins/E/description,proteins/E/coordinates/from,proteins/E/coordinates/to,proteins/E/g.coordinates/from,proteins/E/g.coordinates/to,proteins/M/name,proteins/M/gene,proteins/M/description,proteins/M/coordinates/from,proteins/M/coordinates/to,proteins/M/g.coordinates/from,proteins/M/g.coordinates/to,proteins/ORF6/name,proteins/ORF6/gene,proteins/ORF6/description,proteins/ORF6/coordinates/from,proteins/ORF6/coordinates/to,proteins/ORF6/g.coordinates/from,proteins/ORF6/g.coordinates/to,proteins/ORF7a/name,proteins/ORF7a/gene,proteins/ORF7a/description,proteins/ORF7a/coordinates/from,proteins/ORF7a/coordinates/to,proteins/ORF7a/g.coordinates/from,proteins/ORF7a/g.coordinates/to,proteins/ORF7b/name,proteins/ORF7b/gene,proteins/ORF7b/description,proteins/ORF7b/coordinates/from,proteins/ORF7b/coordinates/to,proteins/ORF7b/g.coordinates/from,proteins/ORF7b/g.coordinates/to,proteins/ORF8/name,proteins/ORF8/gene,proteins/ORF8/description,proteins/ORF8/coordinates/from,proteins/ORF8/coordinates/to,proteins/ORF8/g.coordinates/from,proteins/ORF8/g.coordinates/to,proteins/N/name,proteins/N/gene,proteins/N/description,proteins/N/coordinates/from,proteins/N/coordinates/to,proteins/N/g.coordinates/from,proteins/N/g.coordinates/to,proteins/ORF10/name,proteins/ORF10/gene,proteins/ORF10/description,proteins/ORF10/coordinates/from,proteins/ORF10/coordinates/to,proteins/ORF10/g.coordinates/from,proteins/ORF10/g.coordinates/to


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merges pangolin output report and metadata file '
                    'using isolate as key')
    parser.add_argument('--g', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--pangolin', type=str, default=None,
                        help='Pangolin report (.csv) format')
    parser.add_argument('--output', type=str, default=None,
                        help='Metadata file (.tsv) format')

    return parser.parse_args()


def write_metadata(dataframe):
    dataframe.to_csv(args.output, sep="\t", compression='gzip',
                     quoting=csv.QUOTE_NONE, index=False, header=True)
    


if __name__ == '__main__':
    args = parse_args()

    metadata_df = pd.read_csv(args.metadata, sep="\t", compression='gzip')
    pangolin_df = pd.read_csv(args.pangolin)

    metadata_df['fasta_header_name'] = metadata_df['fasta_header_name'].str.strip()
    pangolin_df['taxon'] = pangolin_df['taxon'].str.strip()

    merged_df = pd.merge(metadata_df, pangolin_df, left_on='fasta_header_name',
                         right_on='taxon')
    merged_df = merged_df.rename(columns={"lineage": "pango_lineage",
                                          "fasta_header_name": "strain"})

    
    write_metadata(dataframe=merged_df)

