#!/usr/bin/env python3

import argparse
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='List of Variants of Concern and Interest from')
    parser.add_argument('--variants', type=str, default=None,
                        help='WHO variants OR custom variants')
    parser.add_argument('--metadata', type=str, default=None,
                        help='metadata file')
    parser.add_argument('--outfile', type=str, default=None,
                        help='list of lineages in output file')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    who_lineages = []
    variants = pd.read_csv(args.variants, sep="\t",
                           low_memory=False)

    for var in variants["pango_lineage"]:
        print(var)
        if "," in var:
            temp = var.split(",")
            who_lineages.extend(temp)
        else:
            who_lineages.append(var)

    Metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    lineages = Metadata['pango_lineage'].unique()

    parsed_lineages=[]
    for lineage in lineages:
        for who_lin in who_lineages:
            if "*" in who_lin:
                who_lin = who_lin[:-1]

                if isinstance(lineage, str) and lineage.startswith(
                        who_lin):
                    parsed_lineages.append(lineage)
            else:
                if lineage == who_lin:
                    parsed_lineages.append(lineage)

    with open(args.outfile, 'w') as f:
        for item in sorted(parsed_lineages):
            f.write("%s\n" % item)
