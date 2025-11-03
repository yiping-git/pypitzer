#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Spencer Revised Chemical Potential data from Python to JSON format.

This script extracts the chemical potential data from spencer_revised_chemical_potential.py 
and converts it to a structured JSON format with state extraction.
"""

import json
import sys
import os
from pathlib import Path

current_dir = Path(__file__).parent

def extract_chemical_potential_data():
    """Extract chemical potential data from the spencer_revised_chemical_potential.py file."""
    # Import the species data from the module
    sys.path.append(str(current_dir))
    
    try:
        from spencer_revised_chemical_potential import solids as chemical_potential_species
        return chemical_potential_species
    except ImportError as e:
        print(f"Error importing spencer_revised_chemical_potential: {e}")
        return None

def convert_chemical_potential_to_json_format(chemical_potential_data):
    """Convert the chemical potential data to JSON format with state extraction."""
    if not chemical_potential_data:
        return None
    
    json_data = {}
    
    for species_name, coefficients in chemical_potential_data.items():
        # Extract state from species name
        if '(' in species_name and ')' in species_name:
            # Extract the state from parentheses
            state_part = species_name[species_name.find('(')+1:species_name.find(')')]
            # clean_name = species_name[:species_name.find('(')].strip()
            
            # Map state abbreviations to full names
            state_mapping = {
                'AQ': 'aq',
                'l': 'l', 
                'S': 's',
                'g': 'g'
            }
            state = state_mapping.get(state_part.upper(), 's')  # default to solid
        else:
            # No state specified, default to solid
            clean_name = species_name
            state = 's'
        
        # Create the species entry with flattened coefficients
        json_data[species_name] = {
            'state': state,
            'ref': 'Liu et al., 2024',
            'T_range': [-90, 25],  # ℃
            **coefficients  # Flatten the coefficients directly into the object
        }
    
    return json_data

def save_to_json(data, output_file):
    """Save the data to a JSON file with proper formatting."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False, sort_keys=True)

def main():
    """Main function to convert Spencer chemical potential data to JSON."""
    print("Converting Spencer Revised Chemical Potential data to JSON...")
    
    # Extract chemical potential data from Python file
    chemical_potential_data = extract_chemical_potential_data()
    
    if not chemical_potential_data:
        print("Failed to extract chemical potential species data. Exiting.")
        return 1
    
    print(f"Extracted {len(chemical_potential_data)} chemical potential species entries")
    
    # Convert chemical potential data to JSON format
    chemical_potential_json = convert_chemical_potential_to_json_format(chemical_potential_data)
    
    if not chemical_potential_json:
        print("Failed to convert chemical potential data to JSON format. Exiting.")
        return 1
    
    # Create output directory if it doesn't exist
    output_dir = current_dir 

    # Save chemical potential data to JSON file
    chemical_potential_output_file = output_dir / "pypitzer_potential.json"
    save_to_json(chemical_potential_json, chemical_potential_output_file)
    
    print(f"Successfully converted chemical potential data to {chemical_potential_output_file}")
    
    # Count entries in chemical potential data
    chemical_potential_count = len(chemical_potential_json)
    aq_count = sum(1 for species in chemical_potential_json.values() if species.get('state') == 'aq')
    l_count = sum(1 for species in chemical_potential_json.values() if species.get('state') == 'l')
    s_count = sum(1 for species in chemical_potential_json.values() if species.get('state') == 's')
    g_count = sum(1 for species in chemical_potential_json.values() if species.get('state') == 'g')
    
    print(f"Chemical potential JSON file contains:")
    print(f"  - Total species: {chemical_potential_count}")
    print(f"  - Aqueous (aq): {aq_count} species")
    print(f"  - Liquid (l): {l_count} species")
    print(f"  - Solid (s): {s_count} species")
    print(f"  - Gas (g): {g_count} species")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
