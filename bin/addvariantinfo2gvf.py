#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 13:41:42 2023

@author: madeline
This script annotates GVF files with variant information.

The attributes completed by this script are: 
["variant", "variant_type", "voi_designation_date",
"voc_designation_date", "vum_designation_date", "status"]
"""

import argparse
import pandas as pd
import numpy as np
from functions import separate_attributes, rejoin_attributes



def parse_args():
    parser = argparse.ArgumentParser(
        description='Adds variant information to a GVF file')
    parser.add_argument('--ingvf', type=str, default=None,
                        help='Path to a GVF file')
    parser.add_argument('--outgvf', type=str,
                        help='Filename for the output GVF file')
    parser.add_argument('--strain', type=str,
                        default='n/a',
                        help='Lineage; user mode is if strain="n/a"')
    parser.add_argument('--clades', type=str, default='n/a',
                        help='TSV file of WHO strain names and '
                             'VOC/VOI status')
    return parser.parse_args()



def add_variant_information(clade_file, gvf, strain):    
    # get info from clades file
    # load clade-defining mutations file
    ### MZA: need to clean up this and add this into separate function "variant_info" 

    # expand #attributes into columns to fill in separately
    gvf = separate_attributes(gvf)
    
    variant_attributes = ["variant", "variant_type", "voi_designation_date",
                   "voc_designation_date", "vum_designation_date",
                   "status"]
    
    if clade_file=='n/a':
        gvf[[variant_attributes]] = "n/a"

    elif clade_file != 'n/a':
        clades = pd.read_csv(clade_file, sep='\t', header=0)
        clades = clades.replace(np.nan, '', regex=True) #this is needed to append empty strings to attributes (otherwise datatype mismatch error)

        # find the relevant pango_lineage line in the clade file that
        # matches args.strain (call this line "var_to_match")
        ### replace this part with func "parse_variant_file" in functions.py
        #available_strains = parse_variant_file(clades)

        var_to_match = 'None'
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

        # if strain in available_strains:
        if var_to_match !='None':

            print("var_to_match", var_to_match)
            # find the index of the relevant row
            var_index = clades.index[clades['pango_lineage'] == var_to_match].tolist()[0]
            print("var_index", var_index)
            print(clades.loc[var_index])
            # extract status, WHO strain name, etc. from clades file
            who_variant = clades.loc[var_index, 'variant']
            variant_type = clades.loc[var_index, 'variant_type']
            voi_designation_date = clades.loc[var_index, 'voi_designation_date']
            voc_designation_date = clades.loc[var_index, 'voc_designation_date']
            vum_designation_date = clades.loc[var_index, 'vum_designation_date']
            status = clades.loc[var_index, 'status']

            # add attributes from variant info file
            gvf["variant"] = who_variant
            gvf["variant_type"] = variant_type
            gvf["voi_designation_date"] = voi_designation_date
            gvf["voc_designation_date"] = voc_designation_date
            gvf["vum_designation_date"] = vum_designation_date
            gvf["status"] = status

        else:
            gvf[[variant_attributes]] = "n/a"
            
            
        # merge attributes back into a single column
        gvf = rejoin_attributes(gvf, empty_attributes)

                                    
                                    
    return(gvf)



 
if __name__ == '__main__':

    args = parse_args()
    
    gvf_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']
    
    pragmas = pd.DataFrame([['##gff-version 3'], ['##gvf-version '
                                                 '1.10'], [
                                '##species NCBI_Taxonomy_URI=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2697049']])  # pragmas are in column 0

    empty_attributes = 'ID=;Name=;chrom_region=;protein=;ps_filter=;ps_exc=; \
        mat_pep_id=;mat_pep_desc=;mat_pep_acc=; ro=;ao=;dp=;sample_size=; \
        Reference_seq=;Variant_seq=;nt_name=;aa_name=;vcf_gene=; \
        mutation_type=; viral_lineage=;multi_aa_name=;multiaa_comb_mutation=; \
        alternate_frequency=;function_category=;source=; citation=; \
        comb_mutation=;function_description=;heterozygosity=;clade_defining=; \
        variant=;variant_type=;voi_designation_date=;voc_designation_date=; \
        vum_designation_date=;status=;'
    empty_attributes = empty_attributes.replace(" ", "")                                              

    # read in gvf file
    gvf = pd.read_csv(args.ingvf, sep='\t', names=gvf_columns, index_col=False)

    # remove pragmas and original header row
    gvf = gvf[~gvf['#seqid'].astype(str).str.contains("#")]

    # separate pragmas
    #pragmas = gvf[gvf['#seqid'].str.contains("##")]
        
    # add variant info
    variant_annotated_gvf = add_variant_information(
        args.clades, gvf, args.strain)
    
    # add pragmas to df, then save to .gvf
    # columns are now 0, 1, ...
    final_gvf = pd.DataFrame(np.vstack([variant_annotated_gvf.columns,
                                            variant_annotated_gvf]))
    final_gvf = pragmas.append(final_gvf)
    filepath = args.outgvf  # outdir + strain + ".annotated.gvf"
    print("Saved as: ", filepath)
    print("")
    final_gvf.to_csv(filepath, sep='\t', index=False, header=False)

