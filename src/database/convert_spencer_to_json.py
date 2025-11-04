#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Spencer Revised Binary data from Python to JSON format.

This script extracts the species data from spencer_revised_binary.py and converts it
to a structured JSON format for easier use in other applications.
All *float* values are written in scientific notation with 9 significant digits (as strings),
while integers (like eq_num) are left untouched.
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
    """Convert binary species data to nested JSON format."""
    if not species_data:
        return None

    json_data = {}

    for key, value in species_data.items():
        param, species1, species2 = key
        species_pair = ",".join(sorted([species1, species2]))

        if species_pair not in json_data:
            json_data[species_pair] = {
                "ref": "Liu et al., 2024",
                "T_range": [-90, 25],
                "eq_num": 0
            }

        json_data[species_pair][param] = value

    return json_data


def convert_ternary_to_json_format(ternary_data):
    """Convert ternary species data to nested JSON format."""
    if not ternary_data:
        return None

    json_data = {}

    for key, value in ternary_data.items():
        param, species1, species2, species3 = key
        species_triplet = ",".join(sorted([species1, species2, species3]))

        if species_triplet not in json_data:
            json_data[species_triplet] = {
                "ref": "Liu et al., 2024",
                "T_range": [-90, 25],
                "eq_num": 0
            }

        json_data[species_triplet][param] = value

    return json_data

# -------------------- Float formatting helper -------------------- #
def format_floats(obj):
    """
    Recursively format all floats to 9-digit scientific notation as strings,
    but keep integers unchanged.
    """
    if isinstance(obj, float):
        return f"{obj:.9e}"
    elif isinstance(obj, int):  # <-- Leave integers intact
        return obj
    elif isinstance(obj, dict):
        return {k: format_floats(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [format_floats(v) for v in obj]
    else:
        return obj

# -------------------- JSON saving -------------------- #
def save_to_json(data, output_file):
    """
    Save data to JSON while preserving full float precision (~16 digits).
    Integers remain unchanged.
    """
    def convert_floats(obj):
        if isinstance(obj, float):
            # Keep full precision, serialize as number
            return float(f"{obj:.16g}")  # preserves up to 16 significant digits
        elif isinstance(obj, int):
            return obj
        elif isinstance(obj, dict):
            return {k: convert_floats(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_floats(v) for v in obj]
        else:
            return obj

    data_with_full_precision = convert_floats(data)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data_with_full_precision, f, indent=2, ensure_ascii=False, sort_keys=False)

# -------------------- Main workflow -------------------- #
def main():
    print("Converting Spencer Revised data to JSON...")

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

    binary_json = convert_binary_to_json_format(binary_data)
    ternary_json = convert_ternary_to_json_format(ternary_data)

    if not binary_json or not ternary_json:
        print("Failed to convert data to JSON format. Exiting.")
        return 1

    combined_json = {
        "binary": binary_json,
        "ternary": ternary_json
    }

    output_file = current_dir / "pypitzer_parameter.json"
    save_to_json(combined_json, output_file)
    print(f"Successfully converted data to {output_file}")

    # Quick summary
    binary_count = len(combined_json["binary"])
    ternary_count = len(combined_json["ternary"])
    theta_count = sum(1 for bp in combined_json["binary"].values() if "theta" in bp)
    lambda_count = sum(1 for bp in combined_json["binary"].values() if "lambda" in bp)
    psi_count = sum(1 for triplet in combined_json["ternary"].values() if "psi" in triplet)
    zeta_count = sum(1 for triplet in combined_json["ternary"].values() if "zeta" in triplet)

    print("JSON file contains:")
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
