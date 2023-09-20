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
from functions import *
import datetime




def parse_args():
    parser = argparse.ArgumentParser(
        description='List of Variants of Concern and Interest from')
    parser.add_argument('--variants', type=str, default=None,
                        help='WHO variants OR custom variants')
    parser.add_argument('--metadata', type=str, default=None,
                        help='metadata file')
    parser.add_argument('--virusseq', action='store_true', default=False,
                        help='virusseq updated lineages only')                    
    parser.add_argument('--outfile', type=str, default=None,
                        help='list of lineages in output file')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    Metadata = pd.read_csv(args.metadata, compression='gzip', sep="\t", low_memory=False)

    parsed_lineages=[]
    if args.virusseq:
            Metadata['last_updated'] = pd.to_datetime(Metadata['last_updated']).dt.date
            filtered_df = Metadata[Metadata['last_updated'] > datetime.datetime.now().date() - pd.to_timedelta("7day")]
            parsed_lineages=filtered_df['lineage'].unique()

    else:
        metadata_lineages = Metadata['lineage'].unique()

        if not args.variants == None:
            variants = pd.read_csv(args.variants, sep="\t",
                                low_memory=False)
            lineages = parse_variant_file(dataframe = variants)

            for metadata_lineage in metadata_lineages:
                for who_lin in lineages:
                    if "*" in who_lin:
                        who_lin = who_lin[:-1]
                        if isinstance(metadata_lineage, str) and metadata_lineage.startswith(who_lin):
                            parsed_lineages.append(metadata_lineage)
                    else:
                        if metadata_lineage == who_lin:
                            parsed_lineages.append(metadata_lineage)

        else:
            parsed_lineages=metadata_lineages
            parsed_lineages = [x for x in parsed_lineages if str(x) != 'nan']
            
    with open(args.outfile, 'w') as f:
        for item in sorted(set(parsed_lineages)):
            f.write("%s\n" % item)
