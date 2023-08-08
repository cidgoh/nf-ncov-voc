#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 2023

@author: madeline

This script converts a GFF to a dictionary of gene names and positions, and adds that dictionary
to a JSON file like gene_positions.JSON for SARS-CoV-2.

Later on, this script will be extended to also parse protein positions from the GFF.

example start_dict file:
{
  "reference": "Nigera-2018",
  "accession": "NC_063383.1",
  "species": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10244",
  "genome": "ACTG"
 }
"""

import pandas as pd
import numpy as np
from json import load, loads, dumps
from functions import separate_attributes
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts a annotated VCF file to a GVF '
                    'file with functional annotation')
    parser.add_argument('--gff_file', type=str, default=None,
                        help='Path to the GFF with genome annotation')
    parser.add_argument('--start_dict', type=str, default=None,
                        help='Path to the GFF with genome annotation')
    parser.add_argument('--savefile', type=str,
                        default=None, help='JSON filename to save results to')
    return parser.parse_args()


if __name__ == '__main__':
    
    args = parse_args()
    gff_file = args.gff_file #"NC_063383.1.gff"
    out_json = args.savefile #"gene_positions_mpox.json"
    
    gene_colors = [
    "rgb(217, 173, 61)",
     "rgb(80, 151, 186)",
      "rgb(230, 112, 48)",
      "rgb(142, 188, 102)",
      "rgb(229, 150, 55)",
      "rgb(170, 189, 82)",
      "rgb(223, 67, 39)",
      "rgb(196, 185, 69)",
      "rgb(117, 182, 129)",
      "rgb(96, 170, 158)" 
    ]
    
    # read in start_dict
    with open(args.start_dict) as fp:
        start_dict = load(fp)
    
    # read in gff
    gff_columns = ['#seqid', '#source', '#type', '#start', '#end',
                   '#score', '#strand', '#phase', '#attributes']
    gff = pd.read_csv(gff_file, sep='\t', names=gff_columns)
    gff = gff[~gff['#seqid'].astype(str).str.contains("##")]
    
    # keep only rows where #type==gene
    gff = gff[gff['#type'].astype(str)=="gene"]
    gff['#start'] = gff['#start'].astype('Int64')
    gff['#end'] = gff['#end'].astype('Int64')
    
    # separate attributes into individual columns
    gff = separate_attributes(gff)
    
    # get genes info: name and coordinates (gene name is from "Name" attribute)
    genes = gff[['Name', '#start', '#end']]
    genes = gff.reset_index(drop=True)
    genes = genes.drop_duplicates('Name', keep="first")
    
    # rename columns
    genes = genes.rename(columns={"#start":"start", "#end":"end", "Name":"name"})
    
    # add second gene name column
    genes['gene'] = genes['name']
    
    # add gene colors, repeating in order
    genes = genes.join(pd.DataFrame(gene_colors * int(len(genes)/len(gene_colors)+1),
        columns=['color']))
    
    '''
    # add intergenic gene region-specific colors
    # if doing this, should get UTR coordinates from the gff
    gene_specific_colors = {"3' UTR": "black", "5' UTR": "black", "INTERGENIC": "grey"}
    specific_colors_df = pd.DataFrame(np.zeros((3, 5)), columns=['gene', 'name', 'start', 'end', 'color'])
    specific_colors_df['name'] = gene_specific_colors.keys()
    specific_colors_df['gene'] = specific_colors_df['name']
    specific_colors_df['color'] = gene_specific_colors.values()
    specific_colors_df[['start', 'end']] = "NA"
    genes = pd.concat([genes, specific_colors_df], axis = 0)
    '''
    
    # convert gene info from gff to dictionary
    
    gene_info = (genes.groupby(['gene', 'name', 'color'], sort=False) 
          .apply(lambda x: x[['start','end']].to_dict('records')[0]) 
          .reset_index() 
          .rename(columns={0:'coordinates'})
          .groupby(['gene', 'color'], sort=False) 
          .apply(lambda x: x[['name','color','coordinates']].to_dict('records')[0]) 
          .reset_index() 
          .rename(columns={0:'information'})) 
            
    gene_info_dict = dict(zip(gene_info.gene, gene_info.information))
    
    # add gene_info_dict to start_dict
    start_dict["genes"] = gene_info_dict
    
    # convert dictionary to JSON
    json_object = dumps(start_dict, indent=2)  
    #print(json_object)
    
    # write JSON object to file
    with open(out_json, "w") as outfile:
        outfile.write(json_object)
