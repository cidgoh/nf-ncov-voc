#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author: Zohaib Anwar
"""

# Importing necessary libraries
import argparse
import pandas as pd
import csv
import datetime

# Defining a function to parse the command line arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description='List of Variants of Concern and Interest from')
    parser.add_argument('--variants', type=str, default=None,
                        help='WHO variants OR custom variants')
    parser.add_argument('--metadata', type=str, default=None,
                        help='metadata file')
    parser.add_argument('--criteria', type=str, default=None,
                        help='criteria for grouping')                      
    parser.add_argument('--virusseq', action='store_true', default=False,
                        help='virusseq updated lineages only')                    
    parser.add_argument('--outfile', type=str, default=None,
                        help='list of lineages in output file')
    parser.add_argument('--start_date', type=str, default=None,
                        help='start date for grouping'),
    parser.add_argument('--end_date', type=str, default=None,
                        help='end date for grouping')
    parser.add_argument('--logfile', type=str, default=None,
                        help='log of updated dataset')

    return parser.parse_args()

# Defining a function to parse the variants file
def parse_variant_file(dataframe):
    who_lineages = []
    for var in dataframe["pango_lineage"]:
        if "," in var:
            for temp in var.split(","):
                if not "[" in var:
                    who_lineages.append(temp)
                    who_lineages.append(temp+".*")
                    
                else:
                    parent=temp[0]
                    child=temp[2:-1].split("|")
                    for c in child:
                        who_lineages.append(parent + str(c))
                        who_lineages.append(parent+str(c)+".*")
        else:
            who_lineages.append(var)
    return who_lineages

# Defining a function to group the data by lineage


def group_by_lineage(dataframe, log, virusseq, end_date):
    """
    This function groups the input dataframe by lineage and returns a list of parsed lineages.
    If the virusseq argument is provided, it filters the dataframe by last_updated date and returns new lineages and sequences.
    If the variants argument is provided, it parses the variant file and returns lineages that match with metadata lineages.
    If neither virusseq nor variants argument is provided, it returns all metadata lineages.

    Args:
    - dataframe: pandas dataframe containing metadata information
    - log: log file path
    - virusseq: boolean flag indicating whether to filter dataframe by last_updated date
    - args: command line arguments

    Returns:
    - parsed_lineages: list of parsed lineages
    """
    parsed_lineages=[]
    
    # filter dataframe by last_updated date if virusseq argument is provided
    if virusseq:
        dataframe['last_updated'] = pd.to_datetime(dataframe['last_updated']).dt.date
        filtered_df = dataframe[dataframe['last_updated'] > datetime.datetime.strptime(args.end_date, '%Y-%m-%d').date() - pd.to_timedelta("7day")]
        Metadata_lineage = dataframe['alias'].unique()
        parsed_lineages=filtered_df['alias'].unique()
        
        new_seqs=len(filtered_df['isolate'])
        new_lineages = list(set(parsed_lineages) - set(Metadata_lineage))    
    else:
        metadata_lineages = dataframe['alias'].unique()

        # parse variant file and return lineages that match with metadata lineages
        if args.variants is not None:
            variants = pd.read_csv(args.variants, sep="\t", low_memory=False)
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
            
            new_lineages = None

        # return all metadata lineages if neither virusseq nor variants argument is provided
        else:
            parsed_lineages=metadata_lineages
            parsed_lineages = [x for x in parsed_lineages if str(x) != 'nan']
            new_lineages = None
    
    # write to log file
    if new_lineages is not None:
        with open(log, 'w') as f:
            f.write("Total sequences:\t%s\n" % num_seqs)
            f.write("New sequences:\t%s\n" % new_seqs)
            f.write("New lineages:\t%s\n" % new_lineages)
    else:
        with open(log, 'w') as f:
            f.write("Total sequences:\t%s\n" % num_seqs)
            f.write("New sequences:\t%s\n" % num_seqs)
            f.write("New lineages:\t%s\n" % new_lineages)            
    
    return parsed_lineages


def group_by_column(dataframe, column_name, variable=None):
    """
    Groups the given dataframe by the unique values in the specified column.
    
    Args:
    - dataframe: pandas.DataFrame object to group by column
    - column_name: str, name of the column to group by
    - variable: str, optional, value to filter the groups by
    
    Returns:
    - list of str, unique values in the specified column if variable is None
    - list of str, [variable] if variable is in the unique values of the specified column
    """
    
    column_unique = dataframe[column_name].unique()
    groups = []
    if variable is None:
        groups = column_unique
    else:
        if variable in column_unique:
            groups.append(variable)
    
    return groups

import pandas as pd

def group_by_time(dataframe, start_date, end_date):
    """
    Groups the given dataframe by sample collection date within the given start and end dates.
    
    Args:
    dataframe (pandas.DataFrame): The dataframe to group by date.
    start_date (str): The start date in YYYY-MM-DD format.
    end_date (str): The end date in YYYY-MM-DD format.
    
    Returns:
    list: A list of date ranges in the format "start_date:end_date".
    """
    dataframe['sample_collection_date'] = pd.to_datetime(dataframe['sample_collection_date']).dt.date
    groups = []
    date = str(start_date) + "_" + str(end_date)
    groups.append(date)
    
    return groups


if __name__ == '__main__':
    args = parse_args()
    Metadata = pd.read_csv(args.metadata, compression='gzip', sep="\t", low_memory=False)
    num_seqs=len(Metadata)
    
    grouping_criteria = args.criteria

    if args.criteria == "lineage":
        parsed_groups  = group_by_lineage(dataframe = Metadata, log = args.logfile, virusseq = args.virusseq, end_date = args.end_date)
        group="lineage"

    elif args.criteria == "time":
        parsed_groups = group_by_time(dataframe = Metadata, start_date = args.start_date, end_date = args.end_date)
        group="time"

    else:
        parsed_groups = group_by_column(dataframe = Metadata, column_name = args.criteria, variable = args.variable)
        group=args.criteria
        
    with open(args.outfile, 'w') as f:
        for item in sorted(set(parsed_groups)):
            f.write("%s:%s\n" % (group,item))
        