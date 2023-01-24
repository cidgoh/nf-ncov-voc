#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: madeline

This script converts VCF files that have been annotated into GVF
files, including the functional annotation. Required user
input is a VCF file.


"""

import argparse
import pandas as pd
import numpy as np
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a annotated VCF file to a GVF '
                    'file with functional annotation')
    parser.add_argument('--vcffile', type=str, default=None,
                        help='Path to a snpEFF-annotated VCF file')
    parser.add_argument('--functional_annotations', type=str,
                        default=None, help='TSV file of functional '
                                           'annotations')
    parser.add_argument('--size_stats', type=str, default='n/a',
                        help='Statistics file for for size extraction')
    parser.add_argument('--clades', type=str, default=None,
                        help='TSV file of WHO strain names and '
                             'VOC/VOI status')
    parser.add_argument('--clades_threshold', type=float,
                        default=0.75,
                        help='Alternate frequency cutoff for '
                             'clade-defining mutations')
    parser.add_argument('--gene_positions', type=str,
                        default=None,
                        help='gene positions in JSON format')
    parser.add_argument('--names_to_split', type=str,
                        default=None,
                        help='.tsv of multi-aa mutation names to '
                             'split up into individual aa names')
    parser.add_argument('--strain', type=str,
                        default='n/a',
                        help='Lineage; user mode is if strain="n/a"')
    parser.add_argument('--outvcf', type=str,
                        help='Filename for the output GVF file')
    parser.add_argument("--names", help="Save mutation names without "
                                        "functional annotations to "
                                        "TSV files for "
                                        "troubleshooting purposes",
                        action="store_true")
    return parser.parse_args()


def map_pos_to_gene_protein(pos, aa_names, GENE_PROTEIN_POSITIONS_DICT):
    """This function is inspired/lifted from Ivan's code.
    Map a series of nucleotide positions to SARS-CoV-2 genes.
    See https://www.ncbi.nlm.nih.gov/nuccore/MN908947.
    :param pos: Nucleotide position pandas series from VCF
    :param aa_names: aa_names pandas series
    :param GENE_PROTEIN_POSITIONS_DICT: Dictionary of gene positions from cov_lineages
    :type pos: int
    :return: series containing SARS-CoV-2 chromosome region names at each
    nucleotide position in ``pos``
    """
    # make a dataframe of the same length as
    # pos to put gene names in (+ other things)
    df = pos.astype(str).to_frame()

    # loop through genes dict to get gene names
    df["gene_names"] = df["POS"]
    for gene in GENE_PROTEIN_POSITIONS_DICT["genes"]:
        # get nucleotide coordinates for this gene
        start = GENE_PROTEIN_POSITIONS_DICT["genes"][gene]["coordinates"]["from"]
        end = GENE_PROTEIN_POSITIONS_DICT["genes"][gene]["coordinates"]["to"]
        # for all the mutations that are found in this region,
        # assign this gene name
        gene_mask = pos.astype(int).between(start, end, inclusive="both")
        if gene == "Stem-loop":  # no stem_loop entry in SARS-CoV-2.json
            df["gene_names"][gene_mask] = gene + ",3\' UTR"
        else:
            df["gene_names"][gene_mask] = gene
    # label all mutations that didn't belong to any gene as "intergenic"
    df["gene_names"][df["gene_names"].str.isnumeric()] = "intergenic"

    # loop through proteins dict to get protein names
    df["protein_names"] = df["POS"]
    for protein in GENE_PROTEIN_POSITIONS_DICT["proteins"]:
        start = GENE_PROTEIN_POSITIONS_DICT["proteins"][protein][
            "g.coordinates"]["from"]
        end = GENE_PROTEIN_POSITIONS_DICT["proteins"][protein]["g.coordinates"]["to"]
        protein_name = GENE_PROTEIN_POSITIONS_DICT["proteins"][protein]["name"]
        # get protein names for all mutations that are within range
        protein_mask = pos.astype(int).between(start, end, inclusive="both")
        df["protein_names"][protein_mask] = protein_name
    # label all mutations that didn't belong to any protein as "n/a"
    df["protein_names"][df["protein_names"].str.isnumeric()] = "n/a"

    return(df["gene_names"], df["protein_names"])


def clade_defining_threshold(threshold, df, sample_size):
    """Specifies the clade_defining attribute as True if AF >
    threshold, False if AF <= threshold, and n/a if the VCF is for a
    single genome """

    if sample_size == 1:
        df["#attributes"] = df["#attributes"].astype(str) + \
                            "clade_defining=n/a;"
    else:
        df.loc[df.AF > threshold, "#attributes"] = df.loc[
                                                       df.AF > threshold, "#attributes"].astype(
            str) + "clade_defining=True;"
        df.loc[df.AF <= threshold, "#attributes"] = df.loc[
                                                        df.AF <= threshold, "#attributes"].astype(
            str) + "clade_defining=False;"
    return df


gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
               '#score', '#strand', '#phase', '#attributes']
vcf_colnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                'FILTER', 'INFO', 'FORMAT', 'unknown']


def vcftogvf(var_data, strain, GENE_PROTEIN_POSITIONS_DICT, names_to_split, sample_size):
    df = pd.read_csv(var_data, sep='\t', names=vcf_colnames)
    # remove pragmas
    df = df[~df['#CHROM'].str.contains("#")]
    # restart index from 0
    df = df.reset_index(drop=True)

    new_df = pd.DataFrame(index=range(0, len(df)), columns=gvf_columns)

    # fill in first 7 GVF columns, excluding 'type'
    new_df['#seqid'] = df['#CHROM']
    new_df['#source'] = '.'
    new_df['#start'] = df['POS']
    # this needs fixing
    new_df['#end'] = (df['POS'].astype(int) + df['ALT'].str.len() -
                      1).astype(str)
    new_df['#score'] = '.'
    new_df['#strand'] = '+'
    new_df['#phase'] = '.'

    # parse INFO column
    # sort out problematic sites tag formats
    df['INFO'] = df['INFO'].str.replace('ps_filter;', 'ps_filter=;')
    df['INFO'] = df['INFO'].str.replace('ps_exc;', 'ps_exc=;')
    df['INFO'] = df['INFO'].str.replace('=n/a', '')

    # parse EFF entry in INFO
    # series: extract everything between parentheses as elements of a
    # list
    eff_info = df['INFO'].str.findall('\((.*?)\)')
    # take first element of list
    eff_info = eff_info.apply(pd.Series)[0]
    # split at pipe, form dataframe
    eff_info = eff_info.str.split(pat='|').apply(pd.Series)
    # hgvs names
    hgvs = eff_info[3].str.rsplit(pat='c.').apply(pd.Series)
    hgvs_protein = hgvs[0].str[:-1]

    hgvs_nucleotide = 'g.' + hgvs[1] # change 'c.' to 'g.' for nucleotide names

    new_df['nt_name'] = hgvs_nucleotide
    new_df['aa_name'] = hgvs_protein

    # change nucleotide names of the form "g.C*4378A" to g.C4378AN;
    # change vcf_gene to "intergenic" here
    asterisk_mask = hgvs_nucleotide.str.contains('\*')
    hgvs_nucleotide[asterisk_mask] = 'g.' + df['REF'] + df['POS'] + \
                                     df['ALT']
    eff_info[5][asterisk_mask] = "intergenic"

    # use nucleotide name where protein name doesn't exist (for
    # 'Name' attribute)
    Names = hgvs[0].str[:-1]
    # fill in empty protein name spaces with nucleotide names ("c."...)
    Names[~Names.str.contains("p.")] = hgvs_nucleotide

    # take "p." off the protein names
    Names = Names.str.replace("p.", "", regex=True)

    new_df["Names"] = Names

    new_df['vcf_gene'] = eff_info[5]
    new_df['mutation_type'] = eff_info[1]

    new_df["multi_name"] = ''
    new_df["multiaa_comb_mutation"] = ''

    # gene and protein name extraction
    gene_names, protein_names = map_pos_to_gene_protein(
        df['POS'].astype(int), new_df['aa_name'], GENE_PROTEIN_POSITIONS_DICT)
    new_df['#attributes'] = 'chrom_region=' + gene_names + ';'
    new_df['#attributes'] = new_df['#attributes'] + 'protein=' + \
        protein_names + ';'

    # make 'INFO' column easier to extract attributes from
    # split at ;, form dataframe
    info = df['INFO'].str.split(pat=';').apply(pd.Series)
    for column in info.columns:
        split = info[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0]
        if isinstance(title, str):
            title = title.lower()
            content = split[1]
            # ignore "tag=" in column content
            info[column] = content
            # make attribute tag as column label
            info.rename(columns={column: title}, inplace=True)

    # fill in 'type' column
    new_df['#type'] = info['type']

    # add 'INFO' attributes by name
    for column in ['ps_filter', 'ps_exc', 'mat_pep_id',
                   'mat_pep_desc', 'mat_pep_acc']:
        # drop nans if they exist
        info[column] = info[column].fillna('')
        new_df['#attributes'] = new_df['#attributes'].astype(str) + \
            column + '=' + info[column].astype(str) + ';'

    # add ro, ao, dp
    unknown = df['unknown'].str.split(pat=':').apply(pd.Series)

    #if sample_size == 1:
    #    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
    #                            'ro=n/a;ao=n/a;dp=1;'
    #else:
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'ro=' + unknown[3].astype(str) + ';'
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'ao=' + unknown[5].astype(str) + ';'
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'dp=' + info['dp'].astype(str) + ';'

    # add sample_size attribute
    # print(sample_size)
    new_df['#attributes'] = new_df['#attributes'] + "sample_size=" + \
                            str(sample_size) + ';'

    # add alternate frequency (AF) column for clade-defining cutoff (
    # af=ao/dp)
    # if there are no commas anywhere in the 'ao' column, calculate
    # AF straight out
    if unknown[5][unknown[5].str.contains(",")].empty:
        new_df['AF'] = unknown[5].astype(int) / info['dp'].astype(int)
    # if there is a comma, add the numbers together to calculate
    # alternate frequency
    else:
        new_df['added_ao'] = unknown[5].apply(lambda x: sum(map(int,
                                                                x.split(
                                                                    ','))))
        new_df['AF'] = new_df['added_ao'].astype(int) / info[
            'dp'].astype(int)

    # add columns copied straight from Zohaib's file
    for column in ['REF', 'ALT']:
        key = column.lower()
        if key == 'ref':
            key = 'Reference_seq'
        elif key == 'alt':
            key = 'Variant_seq'
        new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                                key + '=' + df[column].astype(str) + ';'

    # split multi-aa names from the vcf into single-aa names (multi-row)
    # load names_to_split spreadsheet
    multiaanames = pd.read_csv(names_to_split, sep='\t', header=0)
    # multi-aa names that are in the gvf (list form)
    names_to_separate = np.intersect1d(multiaanames, new_df[
                                                         'Names'].str[
                                                     2:])
    # now split the rows apart in-place
    for multname in names_to_separate:
        # relevant part of 'split_into' column
        splits = multiaanames[multiaanames['name'] == multname][
            'split_into'].copy()
        # list of names to split into
        splits_list_1 = splits.str.split(pat=',').tolist()[0]
        splits_list = [s.replace("'", '') for s in splits_list_1]
        # new_df index containing multi-aa name
        split_index = new_df.index.get_loc(new_df.index[new_df[
                                                            'Names'].str[
                                                        2:] == multname][
                                               0])
        # copy of rows to alter
        seprows = new_df.loc[[split_index]].copy()
        # delete original combined mutation rows
        new_df = new_df.drop(split_index)

        i = 0
        for sepname in splits_list:
            # single-aa name
            seprows['Names'] = "p." + sepname
            # original multi-aa name
            seprows["multi_name"] = multname
            # other single-aa names corresponding to this multi-aa
            # mutation
            seprows["multiaa_comb_mutation"] = splits.tolist()[
                0].replace("'" + sepname + "'", '').replace(",,",
                                                            ',').replace(
                "','", "', '").strip(',')
            new_df = pd.concat([new_df.loc[:split_index + i], seprows,
                                new_df.loc[
                                split_index + i:]]).reset_index(
                drop=True)
            i += 1

    # add attributes
    new_df['#attributes'] = 'Name=' + new_df["Names"] + ';' + new_df[
        '#attributes'].astype(str)
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'nt_name=' + new_df['nt_name'] + ';'
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'aa_name=' + new_df['aa_name'] + ';'
    # gene names
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'vcf_gene=' + new_df['vcf_gene'] + ';'
    # mutation type
    new_df['#attributes'] = new_df['#attributes'].astype(str) + \
                            'mutation_type=' + new_df[
                                'mutation_type'] + ';'

    # add strain name, multi-aa notes, sample_size
    new_df['#attributes'] = new_df[
                                '#attributes'] + 'viral_lineage=' + strain + ';'
    new_df['#attributes'] = new_df['#attributes'] + "multi_aa_name=" + \
                            new_df["multi_name"] + ';'

    new_df['#attributes'] = new_df[
                                '#attributes'] + "multiaa_comb_mutation=" + \
                            new_df["multiaa_comb_mutation"] + ';'
    new_df['#attributes'] = new_df[
                                '#attributes'] + "alternate_frequency=" + \
                            new_df['AF'].astype(str) + ';'

    # only keep the columns needed for a gvf file, plus
    # multiaa_comb_mutation to add to comb_mutation later
    new_df = new_df[gvf_columns + ['multiaa_comb_mutation', 'AF']]
    # new_df.to_csv('new_df.tsv', sep='\t', index=False, header=False)
    return new_df


# takes 4 arguments: the output df of vcftogvf.py, the functional
# annotation file, the clade defining mutations tsv, the strain name,
# and the names_to_split tsv.


def add_functions(gvf, annotation_file, clade_file, strain):
    attributes = gvf["#attributes"].str.split(pat=';').apply(pd.Series)

    # remember this includes nucleotide names where there are no
    # protein names
    gvf["mutation"] = attributes[0].str.split(pat='=').apply(pd.Series)[1]

    # merge annotated vcf and functional annotation files by
    # 'mutation' column in the gvf
    # load functional annotations spreadsheet
    df = pd.read_csv(annotation_file, sep='\t', header=0)

    for column in df.columns:
        df[column] = df[column].str.lstrip()
        # add functional annotations
    merged_df = pd.merge(df, gvf, on=['mutation'], how='right')

    print(merged_df)
    # collect all mutation groups (including reference mutation) in a
    # column, sorted alphabetically
    # this is more roundabout than it needs to be; streamline with
    # grouby() later
    merged_df["mutation_group"] = merged_df["comb_mutation"].astype(
        str) + ", '" + merged_df["mutation"].astype(str) + "', " + \
                                  merged_df[
                                      'multiaa_comb_mutation'].astype(
                                      str)
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.replace("nan, ", "")
    merged_df["mutation_group"] = merged_df[
        "mutation_group"].str.rstrip(' ').str.rstrip(',')

    # separate the mutation_group column into its own df with one
    # mutation per column
    mutation_groups = merged_df["mutation_group"].str.split(
        pat=',').apply(pd.Series)
    mutation_groups = mutation_groups.apply(
        lambda s: s.str.replace("'", ""))
    mutation_groups = mutation_groups.apply(
        lambda s: s.str.replace(" ", ""))
    # now each mutation has a column instead
    mutation_groups = mutation_groups.transpose()
    # sort each column alphabetically
    sorted_df = mutation_groups

    for column in mutation_groups.columns:
        sorted_df[column] = mutation_groups.sort_values(by=column,
                                                        ignore_index=True)[
            column]
    sorted_df = sorted_df.transpose()

    # since they're sorted, put everything back into a single cell,
    # don't care about dropna
    df3 = sorted_df.apply(lambda x: ','.join(x.astype(str)), axis=1)
    unique_groups = df3.drop_duplicates()
    unique_groups_multicol = sorted_df.drop_duplicates()
    # for sanity checking
    merged_df["mutation_group_labeller"] = df3

    # make a unique id for mutation groups that have all members
    # represented in the vcf
    # for groups with missing members, delete those functional
    # annotations
    merged_df["id"] = 'NaN'
    id_num = 0
    for row in range(unique_groups.shape[0]):
        group_mutation_set = set(unique_groups_multicol.iloc[row])
        # remove nan and 'nan' from set
        group_mutation_set = {x for x in group_mutation_set if (x == x
                                                                and x != 'nan')}
        gvf_all_mutations = set(gvf['mutation'].unique())
        indices = merged_df[merged_df.mutation_group_labeller ==
                            unique_groups.iloc[row]].index.tolist()
        # if all mutations in the group are in the vcf file, include
        # those rows and give them an id
        if group_mutation_set.issubset(gvf_all_mutations):
            merged_df.loc[merged_df.mutation_group_labeller ==
                          unique_groups.iloc[row], "id"] = "ID_" + \
                                                           str(id_num)
            id_num += 1
        else:
            # if not, drop group rows, leaving the remaining indices
            # unchanged
            merged_df = merged_df.drop(indices)

            # change semicolons in function descriptions to colons
    merged_df['function_description'] = merged_df[
        'function_description'].str.replace(';', ':')
    # change heteozygosity column to True/False
    merged_df['heterozygosity'] = merged_df['heterozygosity'] == \
                                  'heterozygous'
    # remove trailing spaces from citation
    merged_df['citation'] = merged_df['citation'].str.strip()
    # add key-value pairs to attributes column
    for column in ['function_category', 'source', 'citation',
                   'comb_mutation', 'function_description',
                   'heterozygosity']:
        key = column.lower()
        # replace NaNs with empty string
        merged_df[column] = merged_df[column].fillna('')
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + key + '=' + merged_df[column].astype(str) + ';'
    
    # get clade_defining status, and then info from clades file
    # load clade-defining mutations file
    if clade_file != 'None':
        clades = pd.read_csv(clade_file, sep='\t', header=0)
        clades = clades.replace(np.nan, '', regex=True)

        # find the relevant pango_lineage line in the clade file that
        # matches args.strain (call this line "var_to_match")

        cladefile_strain = 'None'
        available_strains = []
        for var in clades['pango_lineage'].tolist():
            if "," in var:
                for temp in var.split(","):
                    if "[" not in var:
                        available_strains.append(temp)
                        if strain.startswith(temp):
                            var_to_match = var
                    else:
                        parent = temp[0]
                        child = temp[2:-3].split("|")
                        for c in child:
                            available_strains.append(parent + str(c))
                            available_strains.append(parent + str(c) + ".*")
                            if strain.startswith(parent + str(c)):
                                var_to_match = var
            else:
                available_strains.append(var)
                if strain.startswith(var):
                    var_to_match = var
                    
        #print("var_to_match", var_to_match)


        for strain in available_strains:
            #for pango_strain in strain.replace("*", "").split(','):
            if args.strain.startswith(strain): # this will ignore asterisks: "BQ.1.1" returns "BQ" not "BQ.*" as the cladefile_strain
                cladefile_strain = strain
                #print("cladefile_strain", cladefile_strain)

    # if strain in available_strains:
    if cladefile_strain != 'None':
        # find the index of the relevant row
        var_index = clades.index[clades['pango_lineage'] == var_to_match].tolist()[0]
        # extract status, WHO strain name, etc. from clades file
        who_variant = clades.loc[var_index, 'variant']
        variant_type = clades.loc[var_index, 'variant_type']
        voi_designation_date = clades.loc[var_index, 'voi_designation_date']
        voc_designation_date = clades.loc[var_index, 'voc_designation_date']
        vum_designation_date = clades.loc[var_index, 'vum_designation_date']
        status = clades.loc[var_index, 'status']

        # get True/False/n/a designation for clade-defining status
        merged_df = clade_defining_threshold(args.clades_threshold,
                                             merged_df, sample_size)

        # add remaining attributes from clades file
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + "variant=" + who_variant + ';' + "variant_type=" + \
                                   variant_type + ';' + "voi_designation_date=" + \
                                   voi_designation_date + ';' + \
                                   "voc_designation_date=" + \
                                   voc_designation_date + ';' + \
                                   "vum_designation_date=" + \
                                   vum_designation_date + ';' + \
                                   "status=" + status + ';'
    else:
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + "clade_defining=n/a;" + "who_variant=n/a;" + \
                                   "variant_type=n/a;" + \
                                   "voi_designation_date=n/a;" + \
                                   "voc_designation_date=n/a;" + \
                                   "vum_designation_date=n/a;" + \
                                    "status=n/a;"

    # add ID to attributes
    merged_df["#attributes"] = 'ID=' + merged_df['id'].astype(
        str) + ';' + merged_df["#attributes"].astype(str)

    if args.names:
        # get list of names in tsv but not in functional annotations,
        # and vice versa, saved as a .tsv
        tsv_names = gvf["mutation"].unique()
        functional_annotation_names = df["mutation"].unique()
        print(str(np.setdiff1d(tsv_names,
                               functional_annotation_names).shape[0])
              + "/" + str(tsv_names.shape[0]) + " mutation names were "
                                                "not found in "
                                                "functional_annotations")
        leftover_names = pd.DataFrame({'in_tsv_only': np.setdiff1d(
            tsv_names, functional_annotation_names)})
        leftover_names["strain"] = strain
        '''
        clade_names = clades["mutation"].unique()
        leftover_clade_names = pd.DataFrame({'unmatched_clade_names':np.setdiff1d(clade_names, tsv_names)})
        leftover_clade_names["strain"] = strain
        '''
        return merged_df[gvf_columns], leftover_names, gvf[
            "mutation"].tolist()  # , leftover_clade_names

    else:
        return merged_df[gvf_columns]


def find_sample_size(table, lineage):

    if table != 'n/a':
        strain_tsv_df = pd.read_csv(table, delim_whitespace=True,
                                    usecols=['file', 'num_seqs'])

        # not user mode
        if lineage != 'n/a':
            num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(
                lineage + ".qc.")]['num_seqs'].values
            sample_size = num_seqs[0]

        # user-uploaded vcf
        else:
            filename_to_match = args.vcffile.split(".sorted")[0] \
                # looks like "strain.qc"
            num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(
                filename_to_match)]['num_seqs'].values
            sample_size = num_seqs[0]

    # user-uploaded fasta
    elif table == 'n/a' and lineage == 'n/a':
        sample_size = 'n/a'

    return sample_size



if __name__ == '__main__':

    args = parse_args()
    with open(args.gene_positions) as fp:
        GENE_PROTEIN_POSITIONS_DICT = json.load(fp)
    annotation_file = args.functional_annotations
    clade_file = args.clades

    # make empty list in which to store mutation names from all
    # strains in the folder together
    all_strains_mutations = []
    # empty dataframe to hold unmatched names
    leftover_df = pd.DataFrame()
    # unmatched_clade_names = pd.DataFrame() #empty dataframe to hold
    # unmatched clade-defining mutation names
    pragmas = pd.DataFrame([['##gff-version 3'], ['##gvf-version '
                                                  '1.10'], [
                                '##species NCBI_Taxonomy_URI=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049']])  # pragmas are in column 0

    file = args.vcffile
    sample_size = find_sample_size(args.size_stats, args.strain)
    print("Processing: " + file)

    # create gvf from annotated vcf (ignoring pragmas for now)
    gvf = vcftogvf(file, args.strain, GENE_PROTEIN_POSITIONS_DICT,
                   args.names_to_split, sample_size)
    # add functional annotations
    if args.names:
        '''
        annotated_gvf, leftover_names, mutations, \
        leftover_clade_names = add_functions(gvf,
                                             annotation_file,
                                             clade_file, args.strain)
        '''
        annotated_gvf, leftover_names, mutations = add_functions(gvf,
                                                                 annotation_file,
                                                                 clade_file,
                                                                 args.strain)
    else:
        annotated_gvf = add_functions(gvf, annotation_file,
                                      clade_file, args.strain)
    # add pragmas to df, then save to .gvf
    # columns are now 0, 1, ...
    annotated_gvf = pd.DataFrame(np.vstack([annotated_gvf.columns,
                                            annotated_gvf]))
    final_gvf = pragmas.append(annotated_gvf)
    filepath = args.outvcf  # outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)

    # get name troubleshooting reports
    if args.names:
        all_strains_mutations.append(mutations)
        leftover_df = leftover_df.append(leftover_names)
        # unmatched_clade_names = unmatched_clade_names.append(
        # leftover_clade_names)
        # save unmatched names (in tsv but not in
        # functional_annotations) across all strains to a .tsv file
        leftover_names_filepath = "leftover_names.tsv"
        leftover_df.to_csv(leftover_names_filepath, sep='\t',
                           index=False)
        print("")
        print("Mutation names not found in functional annotations "
              "file saved to " + leftover_names_filepath)
        '''
        #save unmatched clade-defining mutation names to a .tsv file
        leftover_clade_names_filepath = "leftover_clade_defining_names.tsv"
        unmatched_clade_names.to_csv(leftover_clade_names_filepath,sep='\t', index=False)
        print("Clade-defining mutation names not found in the annotated VCF saved to " + leftover_clade_names_filepath)
        '''
        # print number of unique mutations across all strains
        flattened = [val for sublist in all_strains_mutations for val in
                     sublist]
        arr = np.array(flattened)
        print("# unique mutations in VCF file: ",
              np.unique(arr).shape[0])

    print("")
    print("Processing complete.")
