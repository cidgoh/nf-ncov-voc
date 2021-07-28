#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extracts Variants of Concern from GISAID Metadata '
                    'and Sequences')
    parser.add_argument('--table', type=str, default=None,
                        help='GISAID Metadata file (.tsv) format')
    parser.add_argument('--fasta', type=str, default=None,
                        help='GISAID Sequence file in .fasta or .fa ('
                             'optionally compressed with gzip) '
                             'format')
    #parser.add_argument('--outdir', type=str, default='./',
    #                    help='Output directpry')
    parser.add_argument('--voc', type=str, default=None,
                        help='Variant of Concern e.g. B.1.1.7, P.1 etc')
    parser.add_argument('--samplingsize', type=int, default=0,
                        help='Sample size, if "0" all sequences '
                             'extracted')
    parser.add_argument('--startdate', type=str, default='2020-01-01',
                        help='Date of submission from (yyyy-mm-dd)')
    parser.add_argument('--enddate', type=str, default='2020-01-01',
                        help='Date of submission to (yyyy-mm-dd)')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    Metadata = pd.read_csv(args.table, sep="\t", low_memory=False,
                           parse_dates=['date_submitted'])
    Metadata['date_submitted'] = pd.to_datetime(Metadata['date_submitted'],
                                      format='%Y-%m-%d',
                                      errors='coerce')

    #print(Metadata['pango_lineage'].value_counts())
    """ Filtering for human associated and consensus sequence of
    at least 29Kb """
    Human_associated = Metadata[(Metadata['host'].str.lower() ==
                                 'Human'.lower()) & (Metadata['length'] >=
                                             29000)]

    sdate = pd.to_datetime(args.startdate).date()
    edate = pd.to_datetime(args.enddate).date()
    Human_associated = Human_associated[
        Human_associated['date_submitted'].isin(
            pd.date_range(sdate, edate))]
    #print(Human_associated['Pango lineage'].value_counts())
    """ Filtering from Main dataset """
    # for voc in voc_list:
    df = Human_associated[Human_associated['pango_lineage'] == args.voc]
    if (args.samplingsize > 0) and (df.shape[0] > args.samplingsize):
        df = df.sample(n=args.samplingsize, replace=False)
    if args.samplingsize == 0:
        sampling = ""
    else:
        sampling = str(args.samplingsize)
    #df1 = df[['Virus name']].copy()
    #df1 = df.filter(['Virus name'], axis=1)
    #print(df1.size)
    #df1['Virus name']=df1['Virus name'].str.split("/")[2]
    ids = df['strain'].tolist()
    #ids=[i.rsplit('/')[1] for i in ids]
    #with open(args.outdir + args.voc + "_ids_" + sampling + "GISAID.txt", 'w') as filehandle:
    #    filehandle.writelines("%s\n" % id for id in ids)
    with open(args.voc + ".txt", 'w') as filehandle:
        filehandle.writelines("%s\n" % id for id in ids)

    # df.to_csv(
    #    args.outdir + args.voc + "_ids_" + sampling + "GISAID.txt",
    #    sep="\t",
    #    columns=["Virus name"], quoting=csv.QUOTE_NONE, index=False,
    #    header=False)
    #df.to_csv(
    #    args.outdir + args.voc + "_" + sampling + "GISAID.tsv",
    #    sep="\t",
    #    quoting=csv.QUOTE_NONE, index=False, header=True)
    df.to_csv(
        args.voc + "_Metadata.tsv",
        sep="\t",
        quoting=csv.QUOTE_NONE, index=False, header=True)
