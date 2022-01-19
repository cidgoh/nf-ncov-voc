#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: zohaib

This script extracts metadata for each VOC, VOI and VUM from the
provided Metadata file based on the assigned lineages. This script
also filters sequences based on the provided criteria.


"""

import argparse
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extracts Variants of Concern and Interest from '
                    'Metadata file')
    parser.add_argument('--table', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--voc', type=str, default=None,
                        help='VOC e.g. B.1.1.7')
    parser.add_argument('--samplingsize', type=int, default=0,
                        help='Sample size, if "0" all sequences '
                             'extracted; Default=0')
    parser.add_argument('--startdate', type=str, default=None,
                        help='Date of submission from (yyyy-mm-dd); '
                             'Default=None')
    parser.add_argument('--enddate', type=str, default=None,
                        help='Date of submission to (yyyy-mm-dd); '
                             'Default=None')
    return parser.parse_args()


def sub_sampling(dataframe, subsampling):
    if (subsampling > 0) and (dataframe.shape[0] > subsampling):
        dataframe = dataframe.sample(n=subsampling, replace=False)
    return dataframe


def write_ids(dataframe):
    ids = dataframe['strain'].tolist()
    with open(args.voc + ".txt", 'w') as \
            filehandle:
        filehandle.writelines("%s\n" % id for id in ids)


def write_metadata(dataframe):
    dataframe.to_csv(args.voc +
                     "_Metadata.tsv", sep="\t",
                     quoting=csv.QUOTE_NONE, index=False, header=True)


def data_filtering(dataframe):
    dataframe = dataframe[dataframe['host (scientific name)'].str.lower() ==
                          'Homo sapiens'.lower()]
    if 'length' in dataframe.columns:
        dataframe = dataframe[dataframe['length'] >= 29000]
    if (not args.startdate == None) and (not args.enddate == None) \
            and ('sample collection date' in dataframe.columns):
        sdate = pd.to_datetime(args.startdate).date()
        edate = pd.to_datetime(args.enddate).date()
        dataframe = dataframe[dataframe[
            'sample collection date'].isin(pd.date_range(sdate, edate))]
    return dataframe


if __name__ == '__main__':
    args = parse_args()

    Metadata = pd.read_csv(args.table, sep="\t", low_memory=False,
                           parse_dates=['sample collection date'])

    if 'sample collection date' in Metadata.columns:
        Metadata['sample collection date'] = pd.to_datetime(Metadata[
                                            'sample collection date'],
                                            format='%Y-%m-%d',
                                            errors='coerce')

    """ Filtering for human associated and consensus sequence of
        at least 29Kb """
    Metadata = Metadata[Metadata['pango_lineage'] == args.voc]
    Metadata = data_filtering(dataframe=Metadata)
    Metadata = sub_sampling(dataframe=Metadata,
                            subsampling=args.samplingsize)
    write_ids(dataframe=Metadata)
    write_metadata(dataframe=Metadata)
