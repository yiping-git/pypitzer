import json
import re
from pathlib import Path

current_dir = Path(__file__).parent
json_file = current_dir / "pypitzer_solids.json"

with open(json_file, "r", encoding="utf-8") as f:
    reactions_data = json.load(f)


def clean_state_marks(text: str) -> str:
    """Remove state marks like (aq), (s), (l), (g)."""
    return re.sub(r'\([a-zA-Z]+\)', '', text)


def parse_species(species_str: str):
    """
    Parse species like '2K+' or 'Fe+2' or 'H2O' to get value and type.
    Removes state marks beforehand.
    """
    s = species_str.strip()
    if not s:
        raise ValueError("Empty species string")

    m = re.match(r'^(\d+)\s*(.+)$', s)
    if m:
        coef = int(m.group(1))
        name = m.group(2).strip()
    else:
        coef = 1
        name = s

    if '+' in name and '-' not in name:
        type_ = 'cation'
    elif '-' in name and '+' not in name:
        type_ = 'anion'
    elif '+' in name and '-' in name:
        if name.endswith('+'):
            type_ = 'cation'
        elif name.endswith('-'):
            type_ = 'anion'
        else:
            type_ = 'neutral'
    else:
        type_ = 'neutral'

    return name, {'value': coef, 'type': type_}


def get_solid_stoichiometry(solid_name: str, reactions_dict: dict):
    """
    Return only the RHS species dictionary (no outer solid key).
    """
    if solid_name not in reactions_dict:
        raise ValueError(f"{solid_name} not found in reactions dictionary")

    reaction_str = reactions_dict[solid_name]
    lhs, rhs = reaction_str.split('=', 1)

    rhs_no_states = clean_state_marks(rhs)
    species_tokens = [tok.strip() for tok in re.split(r'\s+\+\s+', rhs_no_states.strip()) if tok.strip()]

    stoich = {}
    for tok in species_tokens:
        name, info = parse_species(tok)
        stoich[name] = info

    return stoich


# Example usage
if __name__ == '__main__':
    solid_input = "2KCl·FeCl2·2H2O(s)"
    stoich_data = get_solid_stoichiometry(solid_input, reactions_data)
    print(json.dumps(stoich_data, indent=2))
