#!/usr/bin/env python3

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
    parser.add_argument('--gisaid', type=bool,
                        choices=[True], help='Bool Value to '
                                             'confirm if the '
                                             'metadata file is '
                                             'from GISAID; '
                                             'Default=True')
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
    with open(dataframe['pango_lineage'].iloc[1] + ".txt", 'w') as \
            filehandle:
        filehandle.writelines("%s\n" % id for id in ids)


def write_metadata(dataframe):
    dataframe.to_csv(dataframe['pango_lineage'].iloc[1] +
                     "_Metadata.tsv", sep="\t",
                     quoting=csv.QUOTE_NONE, index=False, header=True)


def data_filtering(dataframe):
    dataframe = dataframe[(dataframe['host'].str.lower() ==
                           'Human'.lower()) & (
                                  dataframe['length'] >=
                                  29000)]
    if (not args.startdate == None) and (not args.enddate == None):
        sdate = pd.to_datetime(args.startdate).date()
        edate = pd.to_datetime(args.enddate).date()
        dataframe = dataframe[dataframe[
            'date_submitted'].isin(pd.date_range(sdate, edate))]
    return dataframe


if __name__ == '__main__':
    args = parse_args()

    if not args.gisaid:
        Metadata = pd.read_csv(args.table, sep="\t", low_memory=False)
        df = Metadata[Metadata['pango_lineage'] == args.voc]
        df = sub_sampling(dataframe=df,subsampling=args.samplingsize)
        write_ids(dataframe=df)
        write_metadata(dataframe=df)

    else:
        Metadata = pd.read_csv(args.table, sep="\t", low_memory=False,
                               parse_dates=['date_submitted'])

        Metadata['date_submitted'] = pd.to_datetime(Metadata[
                                    'date_submitted'],
                                    format='%Y-%m-%d', errors='coerce')

        """ Filtering for human associated and consensus sequence of
            at least 29Kb """
        df = Metadata[Metadata['pango_lineage'] == args.voc]
        df = data_filtering(dataframe=df)
        df = sub_sampling(dataframe=df, subsampling=args.samplingsize)
        write_ids(dataframe=df)
        write_metadata(dataframe=df)