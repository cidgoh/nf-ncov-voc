#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import argparse
import math


def add_color(dict, colors):
    # Add the color key-value pair to feature_dict based on feature type
    if dict["type"] == "CDS":
        used_color = colors.pop(0)
        dict["color"] = used_color
        colors.append(used_color)
    elif dict["type"] in ["stem_loop", "five_prime_UTR", "three_prime_UTR"]:
        dict["color"] = "rgb(0, 0, 0)"
    elif dict["type"] == "mature_protein_region_of_CDS":
        dict["color"] = "rgb(128,0,128)"
    
    return dict 


def add_species_pragma(gff_list):
    # Add the pragma key-value pair to feature_dict based on feature type
    for line in gff_list:
            if line.startswith("##species"):
                species = line.split(" ")[1].strip()
    return species


def add_accession_pragma(gff_list):
    # Add the pragma key-value pair to feature_dict based on feature type
    for line in gff_list:
            if line.startswith("##sequence-region"):
                accession = line.split(" ")[1].strip()
    return accession


def add_alias(dict, alias):
    # Add an alias key-value pair to feature_dict based on feature type
    if not alias is None:
            if dict["type"] in ["CDS", "mature_protein_region_of_CDS"]:
                if dict["product"] in alias:
                    dict["protein_alias"] = alias[dict["product"]]
    
    return dict 


def add_protein_coordinates(dict):
    # Add protein coordinates to the dictionary of features
    if dict["type"] == "mature_protein_region_of_CDS":
        aa_start = int(dict["ID"].split(":")[1].split("..")[0])
        aa_end = int(dict["ID"].split(":")[1].split("..")[1])
        dict["aa_start"]  = aa_start
        dict["aa_end"]  = aa_end
    elif dict["type"] == "CDS":
        aa_start = 1
        aa_end = math.floor((dict["end"] - dict["start"])/3)
        dict["aa_start"]  = aa_start
        dict["aa_end"]  = aa_end
    
    return dict


def gff_to_json(gff_file_path, json_file_path, colors_list, alias_dic):
    # Open the GFF file
    with open(gff_file_path, "r") as gff_file:
        # Create a dictionary to hold the features
        features = {}
        id_list = []
        intergenic_dic = {}
        intergenic_dic["type"] = "INTERGENIC"
        intergenic_dic["color"] = "rgb(128,128,128)"
        
        # Iterate over each line in the GFF file
        for line in gff_file:
            # Skip comment lines
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    accession = line.split(" ")[1].strip()
                elif line.startswith("##species"):
                    species = line.split(" ")[1].strip()
                continue
            
        
            # Split the line into fields
            fields = line.strip().split("\t")

            # Create a dictionary to hold the feature information
            feature_dict = {}
            #print(fields)
            feature_dict["type"] = fields[2]
            feature_dict["start"] = int(fields[3])
            feature_dict["end"] = int(fields[4])

            # Parse the attributes field into a dictionary
            attributes = {}
            for attribute in fields[8].split(";"):
                if attribute.strip() == "":
                    continue
                key, value = attribute.split("=")
                attributes[key] = value
            feature_dict.update(attributes)

            # List of dictionaries to be added to the feature dictionary
            if feature_dict["ID"] in id_list:
                feature_dict["product"] = feature_dict["product"]+"-i"
            id_list.append(feature_dict["ID"])
            
            # Post-processing the dictionary of each feature
            
            feature_dict = add_protein_coordinates(dict=feature_dict)
            feature_dict = add_color(dict=feature_dict, colors=colors_list)
            feature_dict = add_alias(dict=feature_dict, alias=alias_dic)
            
            if feature_dict["type"] == "region":
                feature_dict["species"] = species
                feature_dict["accession"] = accession

            # Add the feature dictionary to the dictionary of features
            
            print(feature_dict)

            if feature_dict["type"] == "CDS":
                if "protein_alias" in feature_dict:
                    features[feature_dict["protein_alias"]] = feature_dict
                else:
                    if feature_dict["product"] == "ORF1a polyprotein":
                        continue
                    else:
                        features[feature_dict["product"]] = feature_dict
            elif feature_dict["type"] == "mature_protein_region_of_CDS":
                if "protein_alias" in feature_dict and not feature_dict["protein_alias"] in features: 
                    features[feature_dict["protein_alias"]] = feature_dict
                elif "protein_alias" not in feature_dict and feature_dict["product"] in features:
                    features[feature_dict["product"]] = feature_dict
            elif feature_dict["type"] == "gene":
                features[feature_dict["ID"]] = feature_dict
            else:
                features[feature_dict["gbkey"]] = feature_dict

        # Add intergenic regions to the dictionary of features
        features[intergenic_dic["type"]] = intergenic_dic        
    # Convert the dictionary of features to a JSON string
    json_str = json.dumps(features, indent=4)

    # Write the JSON string to a file
    with open(json_file_path, "w") as json_file:
        json_file.write(json_str)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts GFF to a dictionary of gene names and positions, \
        and adds these to a JSON file"
    )
    parser.add_argument("--gff_file", type=str, default=None, help="Path to the GFF with genome annotation")
    parser.add_argument("--json_file", type=str, default=None, help="JSON filename to save results to")
    parser.add_argument("--color_file", type=str, default=None, help="JSON file containing color codes for genes")
    parser.add_argument("--alias_file", type=str, default=None, help="JSON file containing alias codes for proteins")
    return parser.parse_args()


if __name__ == "__main__":

    # Parse the command-line arguments
    args = parse_args()
    gff_file_path = args.gff_file
    json_file_path = args.json_file

    if args.color_file:
        with open(args.color_file) as f:
            gene_colors = json.load(f)
    else:
        # Define the colors to use for each gene
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
            "rgb(96, 170, 158)",
        ]

    if args.alias_file:
        with open(args.alias_file) as f:
            alias_dic = json.load(f)
    else:
        alias_dic = None

    # Check if the input file exists
    if not os.path.exists(gff_file_path):
        print(f"Error: Input file '{gff_file_path}' does not exist.")

    # Call the gff_to_json function to convert the GFF file to JSON
    gff_to_json(gff_file_path, json_file_path, gene_colors, alias_dic)

    # Check if the output file exists
    if not os.path.exists(json_file_path):
        print(f"Error: Output file '{json_file_path}' was not created.")
