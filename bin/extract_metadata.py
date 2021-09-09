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
                        help='Variant of Concern e.g. B.1.1.7, '
                             'P.1 etc; ; Default=None')
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


if __name__ == '__main__':
    args = parse_args()
    if not args.gisaid:
        Metadata = pd.read_csv(args.table, sep="\t", low_memory=False)
        df = Metadata[Metadata['pango_lineage'] == args.voc]

        if (args.samplingsize > 0) and (
                df.shape[0] > args.samplingsize):
            df = df.sample(n=args.samplingsize, replace=False)
        if args.samplingsize == 0:
            sampling = ""
        else:
            sampling = str(args.samplingsize)

        ids = df['strain'].tolist()

        with open(args.voc + ".txt", 'w') as filehandle:
            filehandle.writelines("%s\n" % id for id in ids)

        df.to_csv(
            args.voc + "_Metadata.tsv",
            sep="\t",
            quoting=csv.QUOTE_NONE, index=False, header=True)

    else:
        Metadata = pd.read_csv(args.table, sep="\t", low_memory=False,
                               parse_dates=['date_submitted'])
        Metadata['date_submitted'] = pd.to_datetime(
            Metadata['date_submitted'],
            format='%Y-%m-%d',
            errors='coerce')

        """ Filtering for human associated and consensus sequence of
        at least 29Kb """
        Human_associated = Metadata[(Metadata['host'].str.lower() ==
                                     'Human'.lower()) & (
                                            Metadata['length'] >=
                                            29000)]

        sdate = pd.to_datetime(args.startdate).date()
        edate = pd.to_datetime(args.enddate).date()
        Human_associated = Human_associated[
            Human_associated['date_submitted'].isin(
                pd.date_range(sdate, edate))]

        """ Filtering from Main dataset """

        df = Human_associated[
            Human_associated['pango_lineage'] == args.voc]
        if (args.samplingsize > 0) and (
                df.shape[0] > args.samplingsize):
            df = df.sample(n=args.samplingsize, replace=False)
        if args.samplingsize == 0:
            sampling = ""
        else:
            sampling = str(args.samplingsize)
        # df1 = df[['Virus name']].copy()
        # df1 = df.filter(['Virus name'], axis=1)
        # print(df1.size)
        # df1['Virus name']=df1['Virus name'].str.split("/")[2]
        ids = df['strain'].tolist()
        # ids=[i.rsplit('/')[1] for i in ids] with open(args.outdir +
        # args.voc + "_ids_" + sampling + "GISAID.txt", 'w') as
        # filehandle: filehandle.writelines("%s\n" % id for id in ids)
        with open(args.voc + ".txt", 'w') as filehandle:
            filehandle.writelines("%s\n" % id for id in ids)

        df.to_csv(
            args.voc + "_Metadata.tsv",
            sep="\t",
            quoting=csv.QUOTE_NONE, index=False, header=True)
