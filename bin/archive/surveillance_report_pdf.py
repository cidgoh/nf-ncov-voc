#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@authors: madeline & zohaib

This script produces summary of surveillance report in a tex format.
It uses tsv file generated in gvf2tsv as an input and
functional indicator file to produce tex file which can be
used for producing PDF file.

This script is based on Jared Simpson's code in ncov-tools
https://github.com/jts/ncov-tools/blob/master/workflow/scripts/generate_report.py

"""

import argparse
import pandas as pd
import os
import argparse
from datetime import datetime
import glob
import csv
import sys
import re
from collections import OrderedDict
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarizes raw surveillance file (TSV) into '
                    'indicator centric and mutation centric readable '
                    'PDF report')
    parser.add_argument('--tsv', type=str, default=None,
                        help='Path to surveillance report TSV file')
    parser.add_argument('--functions_table', type=str, default=None,
                        help='TSV file containing Pokay '
                             'Functional categories:Indicators '
                             'category mappings')
    parser.add_argument('--metadata', type=str, default='n/a',
                        help='Metadata file for contextual data')
    parser.add_argument('--frequency_threshold', type=float,
                        default=0.25,
                        help='Alternate frequency threshold cutoff '
                             'for inclusion in report')
    parser.add_argument('--virusseq', type=bool, default=False,
                        help='VirusSeq dataset')

    return parser.parse_args()


def summarize_functions(tsv, functions_df_template):
    # load functions_df template
    df = pd.read_csv(functions_df_template, sep='\t', header=0)

    # populate the Mutations column row by row
    for row in df['Sub-categories from POKAY']:
        row_mutations_set = set()
        # get list of Pokay categories
        category_list = row.split(',')
        for category in category_list:
            # remove trailing spaces to enable matching
            category = category.rstrip()
            # find mutation names that match that category and add
            # them to the set
            cat_mutations = tsv[tsv['function_category'] == category][
                'name']
            cat_mutations_set = set(cat_mutations)
            # add category mutations set to row mutations set
            row_mutations_set.update(cat_mutations_set)
        # save row mutations set in the 'Mutations' column, sorted
        # alphabetically
        row_list = sorted(list(row_mutations_set))
        row_str = ', '.join(str(e) for e in row_list)
        mask = df['Sub-categories from POKAY'] == row
        df.loc[mask, 'Mutations'] = row_str

    return df


def calculate_frequency(dataframe, newcolname, divisor):
    if dataframe['ao'][dataframe['ao'].astype(str).str.contains(
            ",")].empty:  # if there are no commas anywhere in the 'ao' column, calculate AF straight out
        dataframe[newcolname] = dataframe['ao'].astype(int) / dataframe[
            divisor].astype(int)
    else:  # if there is a comma, add the numbers together to calculate alternate frequency
        # tsv_df['added_ao'] = tsv_df['ao'].astype(str).apply(lambda x: sum(map(int, x.split(','))))
        split_series = dataframe['ao'].str.split(pat=',').apply(
            pd.Series)
        # rename series columns to 'ao_0', 'ao_1', etc.
        split_series.columns = ['ao_' + str(name) for name in
                                split_series.columns.values]
        # ensure all counts are numeric
        for column in split_series.columns:
            split_series[column] = pd.to_numeric(split_series[column],
                                                 errors='coerce')
        # append series to tsv_df
        dataframe = pd.concat([dataframe, split_series], axis=1)
        # calculate frequency for each column
        colnames = [i for i in dataframe.columns.values.tolist() if
                    'ao_' in i]
        for column in colnames:
            dataframe[column] = dataframe[column] / dataframe[
                divisor].astype(int)
        # combine columns into one column
        dataframe[newcolname] = dataframe[colnames].apply(
            lambda row: ','.join(row.values.astype(str)), axis=1)
        # drop split columns
        dataframe = dataframe.drop(labels=colnames, axis=1)
        dataframe[newcolname] = dataframe[newcolname].str.replace(
            ',nan', '')

    return dataframe


def add_source_hyperref(dataframe):
    dataframe["Citation"] = "\href{" + dataframe["Citation URL"].astype(
        str) + "}{" + dataframe["Citation"] + "}"
    return dataframe


def summarize_mutations(tsv, functions_dataframe):
    named_mutations = functions_dataframe['Mutations'].values.tolist()
    named_mutations = ', '.join(str(e) for e in named_mutations).split(
        ', ')
    named_mutations = set(named_mutations)
    if metadata != 'n/a':
        tsv_df_cols = ['name', 'function_category',
                       'function_description', 'viral_lineages',
                       'citation', 'ao', 'dp', 'reference_seq',
                       'variant_seq', 'citation_url']
    else:
        tsv_df_cols = ['name', 'function_category',
                       'function_description',
                       'citation', 'ao', 'dp', 'reference_seq',
                       'variant_seq', 'citation_url']

    # create empty dataframe
    df = pd.DataFrame(columns=tsv_df_cols)

    for mutation in named_mutations:
        # get rows of the tsv for that mutation
        tsv_rows = tsv[tsv['name'] == mutation]
        # keep certain columns of the tsv rows
        tsv_rows = tsv_rows[tsv_df_cols]
        # concatenate dfs
        df = pd.concat((df, tsv_rows))

    # remove clade-defining values from strains column
    # renaming 'viral_clade_defining' to 'viral_lineages'
    if metadata != 'n/a':
        df['viral_lineages'] = df[
            'viral_lineages'].str.replace(r"=.*?;", ",", regex=True)
        # remove trailing commas
        df['viral_lineages'] = df[
            'viral_lineages'].str.rstrip(' ').str.rstrip(',')

        # rename mutations_df columns
        # Removing 'Frequency (Variant)' for now
        # Renaming 'dp' to 'sequence_depth'
        final_mutations_df_cols = ['Mutations', 'Sub-category',
                                   'Function', 'Lineages', 'Citation',
                                   'Alternate Allele Obs',
                                   'Sequence Depth', 'Reference Allele',
                                   'Alternate Allele', 'Citation URL']
    else:
        final_mutations_df_cols = ['Mutations', 'Sub-category',
                                   'Function', 'Citation',
                                   'Alternate Allele Obs',
                                   'Sequence Depth', 'Reference Allele',
                                   'Alternate Allele', 'Citation URL']
    renaming_dict = dict(zip(tsv_df_cols, final_mutations_df_cols))
    df = df.rename(columns=renaming_dict)

    # add 'Frequency (Functional)' column
    # adding if condition
    # if there are no commas
    # anywhere in the 'ao' column, calculate AF straight out
    #if df['Alternate Allele Obs'][df['Alternate Allele Obs'].astype(
    #        str).str.contains("n/a")].empty:
    #    df['Alternate Frequency'] = np.nan
    if df['Alternate Allele Obs'][df['Alternate Allele Obs'].astype(
            str).str.contains(",")].empty:
        df['Alternate Frequency'] = round(
            df['Alternate Allele Obs'].astype(
                int) / df['Sequence Depth'].astype(int), 2)
    else:
        # if there is a comma, add the numbers together to calculate
        # alternate frequency tsv_df['added_ao'] = tsv_df[
        # 'ao'].astype(str).apply(lambda x: sum(map(int, x.split(',
        # '))))
        df['Agg Alternate Allele'] = df['Alternate Allele Obs'].apply(
            lambda x: sum(map(int, x.split(','))))
        # rename series columns to 'ao_0', 'ao_1', etc.
        df['Alternate Frequency'] = round(df['Agg Alternate ' \
                                             'Allele'].astype(
            int) / df['Sequence Depth'].astype(int), 2)

    df = add_source_hyperref(dataframe=df)
    if not df['Alternate Frequency'].isnull().values.any():
        mask = df['Alternate Frequency'] >= args.frequency_threshold
        df = df[mask]
    if metadata != 'n/a':
        mutations_df_cols = ['Mutations', 'Sub-category',
                             'Function', 'Lineages', 'Citation',
                             'Sequence Depth', 'Reference Allele',
                             'Alternate Allele',
                             'Alternate Frequency']
    else:
        mutations_df_cols = ['Mutations', 'Sub-category',
                             'Function', 'Citation',
                             'Sequence Depth', 'Reference Allele',
                             'Alternate Allele', 'Alternate Frequency']
    df = df[mutations_df_cols]
    return df


#
# utility class to assist in converting
# a tsv into a latex table - it can
# rename header columns (name_map)
# perform arbitrary transforms (row_func)
# and filter out columns/rows
#


class TableFormatter:
    def __init__(self):
        self.name_map = dict()
        self.row_func = dict()
        self.column_filter = dict()
        self.row_accept = None
        self.table_spec = ""
        self.size = "normalsize"


# latex starting boilerplate to set up the document class, packages, etc
def write_preamble():
    p = r'''
\documentclass{article}
\usepackage[margin=0.5in, right=1.125in, footskip=35pt]{geometry}
\usepackage{fancyheadings}

\usepackage{lastpage}
\pagestyle{fancy}
\usepackage{hyperref}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{array}
\usepackage{float}
\usepackage{microtype}
\pagestyle{fancy}
\geometry{textwidth=6.0in}


\lhead{\small
        \vspace{-2.0\baselineskip}\href{
        https://github.com/cidgoh/nf-ncov-voc}{nf-ncov-voc}}
\rhead{\small
        \vspace{-2.0\baselineskip}\thepage\ of \pageref{LastPage}}
\cfoot{\small
        \vspace{-2.0\baselineskip}\href{mailto:mzanwar@sfu.ca}{
        Contact Us}}
\rfoot{\small
        \vspace{-2.0\baselineskip}CIDGOH\textsuperscript{
        \textcopyright}}
\renewcommand{\footrulewidth}{1pt}% default is 0pt


\usepackage{longtable}

\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\emergencystretch 3em


\begin{document}



'''
    print(p)


# latex ending boilerplate
def write_postamble():
    p = r'''
\end{document}'''
    print(p)


# count the number of rows in a tsv file
def count_tsv(filename):
    c = 0
    with(open(filename)) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            c += 1
    return c


def escape_latex(s):
    s = s.replace("_", "\_")
    s = s.replace("&", "\&")
    s = s.replace("%", "\%")
    s = s.replace("#", "\#")
    s = s.replace("$", "\$")
    if not "href" in s:
        s = s.replace("{", "\{")
        s = s.replace("}", "\}")
    s = s.replace("~", "\textasciitilde")
    s = s.replace("^", "\textasciicircum}")
    return s


# high-level function to transform a tsv into a latex table
# performing remapping of column names and arbitrary transformation
# of values within each column (e.g. escaping characters latex
# doesn't like) using table formatter
def df_to_table(df, table_formatter):
    reader = df.to_dict(into=OrderedDict, orient='records')

    rows = list()
    header = list()

    for row in reader:

        # remove columns
        for k in table_formatter.column_filter:
            del row[k]

        # skip rows that fail the row filter, if any
        if table_formatter.row_accept is not None and not \
                table_formatter.row_accept(
                    row):
            continue

        if len(header) == 0:

            # remap column names
            for k in row.keys():
                if k in table_formatter.name_map:
                    header.append(table_formatter.name_map[k])
                else:
                    header.append(escape_latex(k))

        # transform with row func
        for k in row:
            if k in table_formatter.row_func:
                row[k] = table_formatter.row_func[k](row[k])
        rows.append(row.values())

    # write latex to stdout
    write_table(table_formatter.table_spec, header, rows,
                table_formatter.size)


# latex for displaying a table
def write_table(spec, header, rows, size):
    # print(r"\begin{center}")
    print(r"\%s" % size)

    print(r"\begin{longtable}{%s}" % spec)

    print(r"\hline")
    print(" & ".join(header) + r" \\ \hline")
    print(r"\endhead")
    for r in rows:
        print(" & ".join(
            [escape_latex(str(v)) for v in r]) + r" \\ \hline")

    print(r"\end{longtable}")
    print(r"\normalsize")
    # print(r"\end{center}")


# write the large per-sample QC table
def write_func_summary(df):
    tf = TableFormatter()
    tf.size = "scriptsize"
    tf.table_spec = "{|p{4.0cm}|p{7.0cm}|p{5.0cm}|}"

    print(r"\section*{Indicator}")
    print(
        r"This table contains key indicators "
        r"identified")
    df_to_table(df, tf)


def write_mutation_summary(df):
    tf = TableFormatter()
    tf.size = "scriptsize"
    if metadata != 'n/a':
        tf.table_spec = "{|p{1.2cm}|p{2.5cm}|p{3.3cm}|p{1.8cm}|p{1.5cm}|p{1.0cm}|p{1.6cm}|p{1.3cm}|p{1.3cm}|}"
    else:
        tf.table_spec = "{|p{1.2cm}|p{2.5cm}|p{3.3cm}|p{1.8cm}|p{1.0cm}|p{1.3cm}|p{1.3cm}|p{1.3cm}|}"

    print(r"\section*{Mutation Significance}")
    print(
        r"This table contains key functional impacts of mutations "
        r"identified")
    df_to_table(df, tf)


if __name__ == '__main__':

    args = parse_args()
    functions_template = args.functions_table
    report_tsv = args.tsv
    metadata = args.metadata
    tsv_df = pd.read_csv(report_tsv, sep='\t', header=0)

    if metadata != 'n/a':
        metadata_df = pd.read_csv(metadata, sep="\t", low_memory=False, compression='gzip',
                                  parse_dates=[
                                      'sample_collection_date'])

        metadata_df['sample_collection_date'] = pd.to_datetime(
            metadata_df['sample_collection_date'], format='%Y-%m-%d',
            errors='coerce')
        lineages = []
        for lineage in tsv_df['viral_lineages']:
            lineage = lineage.split(", ")
            lineages.extend(lineage)

        base = os.path.basename(args.tsv)
        variant = base.split('_')[0]
        # make functions_df
        functions_df = summarize_functions(tsv=tsv_df,
                                           functions_df_template=
                                           functions_template)

        # make mutations_df
        mutations_df = summarize_mutations(tsv=tsv_df,
                                           functions_dataframe=functions_df)
        write_preamble()

        print(r"\section*{Surveillance report}")
        print(
            r"Surveillance generated by nf-ncov-voc for %s variant" % (
                variant))

        print(r"\subsection*{Date }")
        print(
            r"This report is generated on %s using %s number of genomes collected between %s and %s "
            % (
                datetime.today().strftime('%Y-%m-%d'),
                len(metadata_df.index),
                pd.to_datetime(metadata_df['sample_collection_date'].min()).date(),
                pd.to_datetime(metadata_df['sample_collection_date'].max()).date()))

        print(r"\section*{Pango Lineages}")
        print(r"{Pango Lineages in this report }%s "
              % sorted(set(lineages)))
        write_func_summary(df=functions_df)
        write_mutation_summary(df=mutations_df)

        if args.virusseq is True:
            print(r"\newpage")
            print("The results here are in whole or "
                  "part based upon data hosted at the "
                  "Canadian VirusSeq Data Portal: "
                  " \href{https://virusseq-dataportal.ca/}{"
                  "https://virusseq-dataportal.ca/}."
                  "We wish to acknowledge the following "
                  "organisations/laboratories for "
                  "contributing data to the Portal: "
                  "Canadian Public Health Laboratory "
                  "Network (CPHLN), CanCOGGeN VirusSeq "
                  "and the list of labs available at "
                  "\href{https://virusseq-dataportal.ca/acknowledgements"
                  "}{https://virusseq-dataportal.ca/acknowledgements})")

        write_postamble()

        # mutations_df.to_csv('mutations.tsv', sep='\t', index=False)
    else:
        # base = os.path.basename(args.tsv)
        # variant = base.split('_')[0]
        # make functions_df
        functions_df = summarize_functions(tsv=tsv_df,
                                           functions_df_template=
                                           functions_template)
        # make mutations_df
        mutations_df = summarize_mutations(tsv=tsv_df,
                                           functions_dataframe=functions_df)
        write_preamble()
        print(r"\section*{Surveillance report}")
        write_func_summary(df=functions_df)
        write_mutation_summary(df=mutations_df)
        write_postamble()
