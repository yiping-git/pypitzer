#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Spencer Revised Binary data from Python to JSON format.

This script extracts the species data from spencer_revised_binary.py and converts it
to a structured JSON format for easier use in other applications.
"""

import json
import sys
import os
from pathlib import Path

current_dir = Path(__file__).parent

def extract_species_data():
    """Extract species data from the spencer_revised_binary.py file."""
    # Import the species data from the module
    sys.path.append(str(current_dir))
    
    try:
        from spencer_revised_binary import species as binary_species
        return binary_species
    except ImportError as e:
        print(f"Error importing spencer_revised_binary: {e}")
        return None

def extract_ternary_data():
    """Extract ternary species data from the spencer_revised_ternary.py file."""
    # Import the species data from the module
    sys.path.append(str(current_dir))
    
    try:
        from spencer_revised_ternary import species as ternary_species
        return ternary_species
    except ImportError as e:
        print(f"Error importing spencer_revised_ternary: {e}")
        return None

def convert_binary_to_json_format(species_data):
    """Convert the binary species data to a nested JSON format grouped by species pairs."""
    if not species_data:
        return None
    
    # Initialize the nested structure
    json_data = {
        "binary": {},
        "theta": {},
        "lambda": {}
    }
    
    for key, value in species_data.items():
        param, species1, species2 = key
        
        # Create species pair key
        species_pair = f"{species1},{species2}"
        
        # Determine which category this parameter belongs to
        if param in ['b0', 'b1', 'b2', 'c_phi']:
            category = "binary"
        elif param == 'theta':
            category = "theta"
        elif param == 'lambda':
            category = "lambda"
        else:
            # Default to binary for unknown parameters
            category = "binary"
        
        # Initialize species pair if it doesn't exist
        if species_pair not in json_data[category]:
            json_data[category][species_pair] = {}
        
        # Add the parameter data
        json_data[category][species_pair][param] = value
        json_data[category][species_pair]['ref'] = 'Liu et al., 2024'
        json_data[category][species_pair]['T_range'] = [-90, 25] # ℃

    return json_data

def convert_ternary_to_json_format(ternary_data):
    """Convert the ternary species data to a nested JSON format grouped by species triplets."""
    if not ternary_data:
        return None
    
    # Initialize the nested structure - all species triplets in one object
    json_data = {}
    
    for key, value in ternary_data.items():
        param, species1, species2, species3 = key
        
        # Create species triplet key
        species_triplet = f"{species1},{species2},{species3}"
        
        # Initialize species triplet if it doesn't exist
        if species_triplet not in json_data:
            json_data[species_triplet] = {
                'ref': 'Liu et al., 2024',
                'T_range': [-90, 25]  # ℃
            }
        
        # Add the parameter data
        json_data[species_triplet][param] = value

    return json_data

def save_to_json(data, output_file):
    """Save the data to a JSON file with proper formatting."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, sort_keys=True)

def main():
    """Main function to convert Spencer data to JSON."""
    print("Converting Spencer Revised data to JSON...")
    
    # Extract binary data from Python file
    binary_data = extract_species_data()
    
    if not binary_data:
        print("Failed to extract binary species data. Exiting.")
        return 1
    
    print(f"Extracted {len(binary_data)} binary species entries")
    
    # Extract ternary data from Python file
    ternary_data = extract_ternary_data()
    
    if not ternary_data:
        print("Failed to extract ternary species data. Exiting.")
        return 1
    
    print(f"Extracted {len(ternary_data)} ternary species entries")
    
    # Convert binary data to JSON format
    binary_json = convert_binary_to_json_format(binary_data)
    
    if not binary_json:
        print("Failed to convert binary data to JSON format. Exiting.")
        return 1
    
    # Convert ternary data to JSON format
    ternary_json = convert_ternary_to_json_format(ternary_data)
    
    if not ternary_json:
        print("Failed to convert ternary data to JSON format. Exiting.")
        return 1
    
    # Create output directory if it doesn't exist
    output_dir = current_dir 

    # Combine binary and ternary data into one structure
    combined_json = {
        "binary": binary_json.get("binary", {}),
        "theta": binary_json.get("theta", {}),
        "lambda": binary_json.get("lambda", {}),
        "ternary": ternary_json
    }

    # Save combined data to JSON file
    combined_output_file = output_dir / "spencer_revised_parameters.json"
    save_to_json(combined_json, combined_output_file)
    
    print(f"Successfully converted data to {combined_output_file}")
    
    # Count entries in each category for interaction parameters
    binary_count = len(combined_json.get("binary", {}))
    theta_count = len(combined_json.get("theta", {}))
    lambda_count = len(combined_json.get("lambda", {}))
    ternary_count = len(combined_json.get("ternary", {}))
    psi_count = sum(1 for triplet in combined_json.get("ternary", {}).values() if 'psi' in triplet)
    zeta_count = sum(1 for triplet in combined_json.get("ternary", {}).values() if 'zeta' in triplet)
    
    print(f"JSON file contains:")
    print(f"Binary parameters:")
    print(f"  - Binary parameters: {binary_count} species pairs")
    print(f"  - Theta parameters: {theta_count} species pairs")
    print(f"  - Lambda parameters: {lambda_count} species pairs")
    print(f"Ternary parameters:")
    print(f"  - Total species triplets: {ternary_count}")
    print(f"  - Psi parameters: {psi_count} triplets")
    print(f"  - Zeta parameters: {zeta_count} triplets")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
