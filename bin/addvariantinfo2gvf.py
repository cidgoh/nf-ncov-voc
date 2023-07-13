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
from functions import separate_attributes, rejoin_attributes, get_variant_info
from functions import empty_attributes, gvf_columns, vcf_columns


def add_variant_information(clade_file, gvf, strain):    
    # get variant info from clades file

    # expand #attributes into columns to fill in separately
    gvf = separate_attributes(gvf)
    
    variant_attributes = ["variant", "variant_type", "voi_designation_date",
                   "voc_designation_date", "vum_designation_date",
                   "status"]
    
    if clade_file=='n/a':
        gvf[[variant_attributes]] = "n/a"
    
    elif clade_file != 'n/a':
        # load variant info file
        clades = pd.read_csv(clade_file, sep='\t', header=0)
        clades = clades.fillna('')
        
        # retrieve relevant variant information from clades file
        x = get_variant_info(strain, clades)
        
        # if the strain is listed in the file,
        # add variant attributes to the GVF
        if x.strain_in_cladefile==True:
            gvf["variant"] = x.who_variant
            gvf["variant_type"] = x.variant_type
            gvf["voi_designation_date"] = x.voi_designation_date
            gvf["voc_designation_date"] = x.voc_designation_date
            gvf["vum_designation_date"] = x.vum_designation_date
            gvf["status"] = x.status
        else:
            gvf[[variant_attributes]] = "n/a"
    
    # merge attributes back into a single column
    gvf = rejoin_attributes(gvf, empty_attributes)

                                    
                                    
    return(gvf)


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

 
if __name__ == '__main__':

    args = parse_args()                                          

    # read in gvf file
    gvf = pd.read_csv(args.ingvf, sep='\t', names=gvf_columns, index_col=False)

    # remove pragmas and original header row
    pragmas = gvf[gvf['#seqid'].astype(str).str.contains("##")]
    pragmas.columns = range(9)
    pragmas = pragmas.fillna('')
    gvf = gvf[~gvf['#seqid'].astype(str).str.contains("#")]
        
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

