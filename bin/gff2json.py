#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import argparse
import math


def add_color(dict, colors, pathogen="sars-cov-2"):
    """Add the color key-value pair to feature_dict based on feature type"""
    if dict["type"] == "CDS":
        used_color = colors.pop(0)
        dict["color"] = used_color
        colors.append(used_color)
    elif dict["type"] in ["stem_loop", "five_prime_UTR", "three_prime_UTR"]:
        dict["color"] = "rgb(0, 0, 0)"
    elif dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
        # For influenza, use different colors for different mature proteins
        if pathogen == "influenza":
            if "product" in dict:
                if "HA1" in dict["product"]:
                    dict["color"] = "rgb(70, 130, 180)"  # Steel blue for HA1
                elif "HA2" in dict["product"]:
                    dict["color"] = "rgb(220, 20, 60)"   # Crimson for HA2
                else:
                    dict["color"] = "rgb(128, 0, 128)"   # Purple for others
            else:
                dict["color"] = "rgb(128, 0, 128)"
        else:
            dict["color"] = "rgb(128, 0, 128)"
    
    return dict 


def add_species_pragma(gff_list):
    """Add the pragma key-value pair to feature_dict based on feature type"""
    for line in gff_list:
        if line.startswith("##species"):
            species = line.split(" ")[1].strip()
            return species
    return None


def add_accession_pragma(gff_list):
    """Add the pragma key-value pair to feature_dict based on feature type"""
    for line in gff_list:
        if line.startswith("##sequence-region"):
            accession = line.split(" ")[1].strip()
            return accession
    return None


def add_alias(dict, alias):
    """Add an alias key-value pair to feature_dict based on feature type"""
    if alias is not None:
        if dict["type"] in ["CDS", "mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
            if "product" in dict and dict["product"] in alias:
                dict["protein_alias"] = alias[dict["product"]]
    
    return dict 


def add_protein_coordinates(dict, pathogen="sars-cov-2"):
    """Add protein coordinates to the dictionary of features"""
    if dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
        # Parse ID format based on pathogen
        if ":" in dict["ID"]:
            coord_part = dict["ID"].split(":")[1]
            if ".." in coord_part:
                aa_start = int(coord_part.split("..")[0])
                aa_end = int(coord_part.split("..")[1])
                dict["aa_start"] = aa_start
                dict["aa_end"] = aa_end
    elif dict["type"] == "CDS":
        aa_start = 1
        aa_end = math.floor((dict["end"] - dict["start"]) / 3)
        dict["aa_start"] = aa_start
        dict["aa_end"] = aa_end
    
    return dict


def gff_to_json(gff_file_path, json_file_path, colors_list, alias_dic, pathogen="sars-cov-2"):
    """Convert GFF file to JSON format"""
    # Open the GFF file
    with open(gff_file_path, "r") as gff_file:
        # Create a dictionary to hold the features
        features = {}
        id_list = []
        intergenic_dic = {}
        intergenic_dic["type"] = "INTERGENIC"
        intergenic_dic["color"] = "rgb(128,128,128)"
        
        species = None
        accession = None
        
        # Iterate over each line in the GFF file
        for line in gff_file:
            # Skip comment lines but extract metadata
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    accession = line.split(" ")[1].strip()
                elif line.startswith("##species"):
                    species = line.split(" ")[1].strip()
                continue
            
            # Split the line into fields
            
            fields = line.strip().split("\t")

            # Skip lines that don't have enough fields (empty or malformed)
            if len(fields) < 9:
                continue

            # Create a dictionary to hold the feature information
            feature_dict = {}
            feature_dict["type"] = fields[2]
            feature_dict["start"] = int(fields[3])
            feature_dict["end"] = int(fields[4])

            # Parse the attributes field into a dictionary
            attributes = {}
            for attribute in fields[8].split(";"):
                if attribute.strip() == "":
                    continue
                if "=" not in attribute:
                    continue
                key, value = attribute.split("=", 1)
                attributes[key] = value
            feature_dict.update(attributes)

            # Handle duplicate IDs
            if feature_dict["ID"] in id_list:
                feature_dict["product"] = feature_dict.get("product", "unknown") + "-i"
            id_list.append(feature_dict["ID"])
            
            # Post-processing the dictionary of each feature
            feature_dict = add_protein_coordinates(dict=feature_dict, pathogen=pathogen)
            feature_dict = add_color(dict=feature_dict, colors=colors_list, pathogen=pathogen)
            feature_dict = add_alias(dict=feature_dict, alias=alias_dic)
            
            if feature_dict["type"] == "region":
                if species:
                    feature_dict["species"] = species
                if accession:
                    feature_dict["accession"] = accession

            # Add the feature dictionary to the dictionary of features
            print(feature_dict)

            # Determine the key to use for this feature
            feature_key = None
            
            if feature_dict["type"] == "CDS":
                # SARS-CoV-2 specific handling
                if pathogen == "sars-cov-2":
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"]
                    else:
                        if "product" in feature_dict and feature_dict["product"] == "ORF1a polyprotein":
                            continue  # Skip ORF1a polyprotein for SARS-CoV-2
                        elif "product" in feature_dict:
                            feature_key = feature_dict["product"]
                        else:
                            feature_key = feature_dict["ID"]
                else:
                    # Influenza and other pathogens
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"]
                    elif "product" in feature_dict:
                        feature_key = feature_dict["product"]
                    elif "gene" in feature_dict:
                        feature_key = feature_dict["gene"]
                    else:
                        feature_key = feature_dict["ID"]
                    
            elif feature_dict["type"] in ["mature_protein_region_of_CDS", "signal_peptide_region_of_CDS"]:
                # SARS-CoV-2 specific handling
                if pathogen == "sars-cov-2":
                    if "protein_alias" in feature_dict and feature_dict["protein_alias"] not in features:
                        feature_key = feature_dict["protein_alias"]
                    elif "protein_alias" not in feature_dict and "product" in feature_dict and feature_dict["product"] in features:
                        feature_key = feature_dict["product"]
                    elif "product" in feature_dict:
                        feature_key = feature_dict["product"]
                    else:
                        feature_key = feature_dict["ID"]
                else:
                    # Influenza and other pathogens - always keep these regions
                    if "protein_alias" in feature_dict:
                        feature_key = feature_dict["protein_alias"]
                    elif "product" in feature_dict:
                        feature_key = feature_dict["product"]
                    else:
                        feature_key = feature_dict["ID"]
                    
            elif feature_dict["type"] == "gene":
                feature_key = feature_dict["ID"]
                
            elif "gbkey" in feature_dict:
                feature_key = feature_dict["gbkey"]
            else:
                feature_key = feature_dict["ID"]
            
            # Skip if feature_key is None (e.g., skipped ORF1a)
            if feature_key is None:
                continue
            
            # Add to features dictionary
            if pathogen == "sars-cov-2":
                # SARS-CoV-2: overwrite duplicates (original behavior)
                if feature_key not in features:
                    features[feature_key] = feature_dict
            else:
                # Influenza and others: keep all regions with unique keys
                if feature_key not in features:
                    features[feature_key] = feature_dict
                else:
                    # If key exists, append a suffix to keep both
                    suffix = 1
                    new_key = f"{feature_key}_{suffix}"
                    while new_key in features:
                        suffix += 1
                        new_key = f"{feature_key}_{suffix}"
                    features[new_key] = feature_dict

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
        and adds these to a JSON file. Supports multiple pathogens including SARS-CoV-2 and Influenza."
    )
    parser.add_argument("--gff_file", type=str, required=True, help="Path to the GFF with genome annotation")
    parser.add_argument("--json_file", type=str, required=True, help="JSON filename to save results to")
    parser.add_argument("--pathogen", type=str, default="sars-cov-2", 
                        choices=["sars-cov-2", "influenza"],
                        help="Pathogen type for proper parsing (default: sars-cov-2)")
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
        exit(1)

    # Call the gff_to_json function to convert the GFF file to JSON
    gff_to_json(gff_file_path, json_file_path, gene_colors, alias_dic, args.pathogen)

    # Check if the output file exists
    if os.path.exists(json_file_path):
        print(f"Success: Output file '{json_file_path}' was created.")
    else:
        print(f"Error: Output file '{json_file_path}' was not created.")