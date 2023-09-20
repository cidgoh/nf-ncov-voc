#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script uses assets/ncov_variants_who_variants.tsv as key and
lists the lineages corresponding to VOCs, VOIs and VUMs in input
dataset for extracting metadata.

"""

import argparse
import pandas as pd
import csv
#from functions import *


def parse_args():
    parser = argparse.ArgumentParser(
        description='convert mpox gisaid to covid viralai format')
    parser.add_argument('--metadata', type=str, default=None,
                        help='metadata file')
    parser.add_argument('--outfile', type=str, default=None,
                        help='list of lineages in output file')
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    
    Metadata = pd.read_csv(args.metadata, compression='gzip', sep="\t", low_memory=False)
    #Metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    Metadata.columns = map(str.lower, Metadata.columns)
    Metadata.columns = Metadata.columns.str.replace(' / ', ' ')
    Metadata[['clade', 'lineage']] = Metadata['clade lineage'].str.split(' ', 1, expand=True)
    Metadata.columns = Metadata.columns.str.replace(' ', '_')
    Metadata = Metadata.rename(columns={'collection_date': 'sample_collection_date', 'host': 'host_scientific_name', 'virus_name': 'isolate'})
    Meatdata = Metadata["lineage"].dropna()
    Metadata = Metadata[~Metadata["lineage"].str.contains("probable", na=False)]
    Metadata = Metadata.replace(to_replace='Human', value='Homo Sapiens')
    #Metadata = Metadata[Metadata['host_scientific_name'] == 'Human'] = "Homo sapiens"

    #print(Metadata.columns)
    #print(Metadata.head())
    
    Metadata.to_csv(args.outfile, sep="\t", compression='gzip',
                    quoting=csv.QUOTE_NONE, index=False, header=True)