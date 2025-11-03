#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Spencer Revised Binary data from Python to JSON format.

This script extracts the species data from spencer_revised_binary.py and converts it
to a structured JSON format for easier use in other applications.
"""

import json
import sys
from pathlib import Path

current_dir = Path(__file__).parent

# -------------------- Data extraction -------------------- #
def extract_species_data():
    """Extract binary species data from spencer_revised_binary.py."""
    sys.path.append(str(current_dir))
    try:
        from spencer_revised_binary import species as binary_species
        return binary_species
    except ImportError as e:
        print(f"Error importing spencer_revised_binary: {e}")
        return None

def extract_ternary_data():
    """Extract ternary species data from spencer_revised_ternary.py."""
    sys.path.append(str(current_dir))
    try:
        from spencer_revised_ternary import species as ternary_species
        return ternary_species
    except ImportError as e:
        print(f"Error importing spencer_revised_ternary: {e}")
        return None

# -------------------- Conversion functions -------------------- #
def convert_binary_to_json_format(species_data):
    """
    Convert binary species data to nested JSON format.
    All binary parameters (b0,b1,b2,c_phi,theta,lambda) are included per species pair.
    """
    if not species_data:
        return None

    json_data = {}

    for key, value in species_data.items():
        param, species1, species2 = key

        # Alphabetically sort species for consistent keys
        species_pair = ",".join(sorted([species1, species2]))

        # Initialize species pair if not exist
        if species_pair not in json_data:
            json_data[species_pair] = {
                "ref": "Liu et al., 2024",
                "T_range": [-90, 25]  # ℃
            }

        # Add parameter data
        json_data[species_pair][param] = value

    return json_data

def convert_ternary_to_json_format(ternary_data):
    """
    Convert ternary species data to nested JSON format.
    Each triplet stores psi or zeta parameters.
    """
    if not ternary_data:
        return None

    json_data = {}

    for key, value in ternary_data.items():
        param, species1, species2, species3 = key

        # Alphabetically sort species for consistent keys
        species_triplet = ",".join(sorted([species1, species2, species3]))

        if species_triplet not in json_data:
            json_data[species_triplet] = {
                "ref": "Liu et al., 2024",
                "T_range": [-90, 25], # ℃
                "eq_num": 0
            }
            
        json_data[species_triplet][param] = value
        

    return json_data

# -------------------- JSON saving -------------------- #
def save_to_json(data, output_file):
    """Save data to JSON with proper formatting."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, sort_keys=False)

# -------------------- Main workflow -------------------- #
def main():
    print("Converting Spencer Revised data to JSON...")

    # Extract binary and ternary data
    binary_data = extract_species_data()
    if not binary_data:
        print("Failed to extract binary species data. Exiting.")
        return 1
    print(f"Extracted {len(binary_data)} binary species entries")

    ternary_data = extract_ternary_data()
    if not ternary_data:
        print("Failed to extract ternary species data. Exiting.")
        return 1
    print(f"Extracted {len(ternary_data)} ternary species entries")

    # Convert to JSON format
    binary_json = convert_binary_to_json_format(binary_data)
    ternary_json = convert_ternary_to_json_format(ternary_data)

    if not binary_json or not ternary_json:
        print("Failed to convert data to JSON format. Exiting.")
        return 1

    # Combine into single JSON structure
    combined_json = {
        "binary": binary_json,
        "ternary": ternary_json
    }

    # Save to file
    output_file = current_dir / "pypitzer_parameter.json"
    save_to_json(combined_json, output_file)
    print(f"Successfully converted data to {output_file}")

    # Count entries for reporting
    binary_count = len(combined_json["binary"])
    ternary_count = len(combined_json["ternary"])
    theta_count = sum(1 for bp in combined_json["binary"].values() if "theta" in bp)
    lambda_count = sum(1 for bp in combined_json["binary"].values() if "lambda" in bp)
    psi_count = sum(1 for triplet in combined_json["ternary"].values() if "psi" in triplet)
    zeta_count = sum(1 for triplet in combined_json["ternary"].values() if "zeta" in triplet)

    print(f"JSON file contains:")
    print(f"Binary parameters: {binary_count} species pairs")
    print(f"  - Theta parameters: {theta_count} pairs")
    print(f"  - Lambda parameters: {lambda_count} pairs")
    print(f"Ternary parameters: {ternary_count} triplets")
    print(f"  - Psi parameters: {psi_count} triplets")
    print(f"  - Zeta parameters: {zeta_count} triplets")

    return 0

# -------------------- Entry point -------------------- #
if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
