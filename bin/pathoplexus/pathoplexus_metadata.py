#!/usr/bin/env python3

import requests
import pandas as pd
import argparse
import json

# API Endpoint
URL = "https://lapis.pathoplexus.org/mpox/sample/details"

# Headers
HEADERS = {
    "accept": "application/json"
}

def fetch_data(config_file):
    """Fetch metadata from LAPIS API
    params = {
        "limit": limit,
        "dataFormat": "JSON",
        "downloadAsFile": "false",
        "dataUseTerms": "OPEN"
    }
    """
    def load_params_from_json(file_path):
        """Load parameters from a JSON file"""
        with open(file_path, 'r') as file:
            return json.load(file)

    params = load_params_from_json(config_file)
    response = requests.get(URL, params=params, headers=HEADERS)
    
    if response.status_code == 200:
        return response.json()  # Returns JSON data
    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None

def save_to_tsv(data, filename):
    """Save JSON data to a TSV file"""
    if not data or "data" not in data:
        print("No valid data to save.")
        return

    # Convert JSON to DataFrame
    df = pd.DataFrame(data["data"])
    df = df.drop_duplicates(subset=['accession'])
    # Save as TSV
    df.to_csv(filename, sep="\t", index=False)
    print(f"Data saved to {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch MPOX metadata from Pathoplexus")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-c", "--config", required=True, help="parameters JSON file")
    args = parser.parse_args()

    metadata = fetch_data(config_file = args.config)
    save_to_tsv(metadata, args.output)