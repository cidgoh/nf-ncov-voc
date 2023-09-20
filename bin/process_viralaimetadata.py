#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description='Download metadata from ViralAi')
    parser.add_argument('--alias', type=str, default=None,
                        help='.tsv file')
    parser.add_argument('--csv', type=str, default=None,
                        help='.csv file')
    parser.add_argument('--outfile', type=str, default=None,
                        help='output file')                        

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    df = pd.read_csv(args.csv, sep=None, compression='gzip', engine='python')
    df['raw_lineage'] = df['lineage']
    
    with open(args.alias, 'r') as j:
        alias_dic = json.loads(j.read())
    
    for k in list(alias_dic.keys()):
        if k.startswith('X'):
            del alias_dic[k]
    
    df['raw_lineage']=(df['raw_lineage'].str.extract(r"([A-Z]+)",
                                                    expand=False)
                                                            .map(alias_dic) + "." + df['raw_lineage'].str.split(".",1).str[1]).fillna(df['raw_lineage'])

    # Sort by sample_collection_date and write it to csv
    df.to_csv(args.outfile, encoding='utf-8', index=False, sep='\t',
              compression='gzip')