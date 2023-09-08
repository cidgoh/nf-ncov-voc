#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script harmonizes metadata file from different sources into a valid format for virus-mvp

"""

import argparse
import pandas as pd
import csv
import yaml



def parse_args():
    parser = argparse.ArgumentParser(
        description='convert metadata files from different sources to general format')
    parser.add_argument('--metadata', type=str, default=None,
                        help='metadata file')
    parser.add_argument('--outfile', type=str, default=None,
                        help='list of lineages in output file')
    parser.add_argument('--config', type=str, default=None,
                        help='config file')
    parser.add_argument('--data', type=str, default=None,
                        help='sata source in yaml file')                     
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    content = None

    with open(args.config) as f:
        try:
            content = yaml.safe_load(f)
        except yaml.YAMLError as e:
            print(e)

    sep=content["Format"]["sep"]
    suffix=content["Format"]["suffix"]
    if(sep=='tab'):
        delim="\t"
    else:
        delim=","
    
    if(suffix=='.gz'):
        compression='gzip'
    else:
        compression=None
    
    Metadata = pd.read_csv(args.metadata, compression=compression, sep=delim, low_memory=False)
    Metadata.columns = map(str.lower, Metadata.columns)
    

    required = ([k for k,v in content["Required"].items() if v == True])
    optional = ([k for k,v in content["Optional"].items() if v == True])
    columns_config = required + optional
    
    data=args.data
    
    columns_meta = []
    for c in columns_config:
        columns_meta.append(content[data][c])
    
    Metadata = Metadata[columns_meta]
    header_dic = {v: k for k, v in content[data].items()}
    
    Metadata.rename(columns=header_dic, inplace=True)
    Metadata.to_csv(args.outfile, sep="\t", compression='gzip',
                    quoting=csv.QUOTE_NONE, index=False, header=True)