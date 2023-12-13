#!/usr/bin/env python
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
import datetime as dt

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extracts Variants of Concern and Interest from '
                    'Metadata file')
    parser.add_argument('--table', type=str, default=None,
                        help='Metadata file (.tsv) format')
    # parser.add_argument('--criteria', type=str, default="lineage",
    #                     help='Criteria for grouping samples together'
    #                          '- e.g., "lineage or time')
    parser.add_argument('--group', type=str, default=None,
                        help='VOC e.g. B.1.1.7')
    parser.add_argument('--samplingsize', type=int, default=0,
                        help='Sample size, if "0" all sequences '
                             'extracted; Default=0')
    parser.add_argument('--start_date', type=str, default=None,
                        help='Date of submission from (yyyy-mm-dd); '
                             'Default=None')
    parser.add_argument('--end_date', type=str, default=None,
                        help='Date of submission to (yyyy-mm-dd); '
                             'Default=None')
    parser.add_argument('--window', type=int, default=7,
                        help='Number of days between start and end date '
                             'Default=7')
    parser.add_argument('--outtable', type=str, default=None,
                        help='Metadata file (.tsv) format')
    parser.add_argument('--outids', type=str, default=None,
                        help='ids file (.txt) format')                                                                       
                                                       
    return parser.parse_args()


def sub_sampling(dataframe, subsampling):
    if (subsampling > 0) and (dataframe.shape[0] > subsampling):
        dataframe = dataframe.sample(n=subsampling, replace=False)
    return dataframe


def write_ids(dataframe, start_date, end_date, criterion):
    ids = dataframe['isolate'].tolist()
    #with open(args.outids, 'w') as filehandle:
    #         filehandle.writelines("%s\n" % id for id in ids)
    if criterion == "lineage":
        with open(args.outids, 'w') as filehandle:
            filehandle.writelines("%s\n" % id for id in ids)
    elif criterion == "time":
        with open( str(start_date)+ "_" + str(end_date) + ".txt", 'w') as filehandle:
            filehandle.writelines("%s\n" % id for id in ids)
    else:
        with open(args.outids, 'w') as filehandle:
            filehandle.writelines("%s\n" % id for id in ids)

def data_filtering(dataframe):
    #if 'geo_loc_name_state_province_territory' in dataframe.columns and not location == None:
    #    dataframe = dataframe[dataframe['geo_loc_name_state_province_territory'].str.lower() ==
    #                            location.lower()]
    #if 'host_scientific_name' in dataframe.columns:
    #    dataframe = dataframe[dataframe['host_scientific_name'].str.lower() ==
    #                      'Homo sapiens'.lower()]
    #if 'length' in dataframe.columns:
    #    dataframe = dataframe[dataframe['length'] >= 29000]
    if (not args.start_date == None) and (not args.end_date == None) \
            and ('sample_collection_date' in dataframe.columns):
        sdate = pd.to_datetime(args.start_date).date()
        edate = pd.to_datetime(args.end_date).date()
        dataframe = dataframe[dataframe[
            'sample_collection_date'].isin(pd.date_range(sdate, edate))]
    return dataframe


def filter_metadata(dataframe):
    if 'host_scientific_name' in dataframe.columns:
        dataframe = dataframe[dataframe['host_scientific_name'].str.lower() ==
                              'Homo sapiens'.lower()]
    if 'length' in dataframe.columns:
        dataframe = dataframe[dataframe['length'] >= 29000]        
    return dataframe

    

def write_metadata(dataframe, start_date, end_date, criterion):
    # dataframe.to_csv(args.outtable, sep="\t", compression='gzip',
    #                   quoting=csv.QUOTE_NONE, index=False, header=True)
    if criterion == "lineage":
        dataframe.to_csv(args.outtable, sep="\t", compression='gzip',
                     quoting=csv.QUOTE_NONE, index=False, header=True)
    elif criterion == "time":
        dataframe.to_csv( str(start_date)+ "_" + str(end_date) + 
                     "_Metadata.tsv.gz", sep="\t", compression='gzip',
                     quoting=csv.QUOTE_NONE, index=False, header=True)
    else:
        dataframe.to_csv(args.outtable, sep="\t", compression='gzip',
                     quoting=csv.QUOTE_NONE, index=False, header=True)




if __name__ == '__main__':
    args = parse_args()
    if not args.group == None:
        group = args.group.split(":")[1]
        criteria = args.group.split(":")[0]
        
    Metadata = pd.read_csv(args.table, sep="\t", low_memory=False, compression='gzip',
                           parse_dates=['sample_collection_date'])

    if 'sample_collection_date' in Metadata.columns:
        Metadata['sample_collection_date'] = pd.to_datetime(Metadata[
                                            'sample_collection_date'],
                                            format='%Y-%m-%d',
                                            errors='coerce')

        sdate = pd.to_datetime(args.start_date, format='%Y-%m-%d')
        edate = pd.to_datetime(args.end_date, format='%Y-%m-%d')
        window=args.window


    """ Filtering for human associated and consensus sequence of
        at least 29Kb """ 
   

    if criteria == "lineage":    
        Metadata = Metadata[Metadata['lineage'] == group]
        
        Metadata = data_filtering(dataframe=Metadata)
        Metadata = sub_sampling(dataframe=Metadata, subsampling=args.samplingsize)
        write_ids(dataframe=Metadata, start_date=sdate, end_date=edate, criterion=criteria)
        write_metadata(dataframe=Metadata, start_date=sdate, end_date=edate, criterion=criteria)
    
    elif criteria == "time":
        while sdate <= edate:
            query_date = sdate + pd.DateOffset(days=window - 1)
            #print(sdate, query_date)
            sub_meta = Metadata.query('sample_collection_date >= @sdate and sample_collection_date <= @query_date')
            #print(len(sub_meta))
            sub_meta = sub_sampling(dataframe=sub_meta, subsampling=args.samplingsize)
            
            write_ids(dataframe=sub_meta, start_date=str(sdate)[:10], end_date=str(query_date)[:10], criterion=criteria)
            write_metadata(dataframe=sub_meta, start_date=str(sdate)[:10], end_date=str(query_date)[:10], criterion=criteria)
            
            sdate += pd.DateOffset(days=window)


    else:
        Metadata = filter_metadata(dataframe=Metadata)
        Metadata = Metadata[Metadata[criteria] == group.replace("_", " ")]
        
        Metadata = data_filtering(dataframe=Metadata)
        Metadata = sub_sampling(dataframe=Metadata, subsampling=args.samplingsize)
        write_ids(dataframe=Metadata, start_date=sdate, end_date=edate, criterion=criteria)
        write_metadata(dataframe=Metadata, start_date=sdate, end_date=edate, criterion=criteria)
        