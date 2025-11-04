import json
import re
from pathlib import Path

current_dir = Path(__file__).parent
json_file = current_dir / "pypitzer_reaction.json"

with open(json_file, "r", encoding="utf-8") as f:
    reaction_database = json.load(f)


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


def get_solid_stoichiometry(reaction_data):
    """
    Return only the RHS species dictionary (no outer solid key).
    """
    reaction_str = reaction_data['reaction']
    if reaction_str:
        lhs, rhs = reaction_str.split('=', 1)

        rhs_no_states = clean_state_marks(rhs)
        species_tokens = [tok.strip() for tok in re.split(r'\s+\+\s+', rhs_no_states.strip()) if tok.strip()]

        stoich = {}
        for tok in species_tokens:
            name, info = parse_species(tok)
            stoich[name] = info
        return stoich
    return None

class Reaction:
    def __init__(self, reaction_key):
        self.reaction_key = reaction_key

    def reaction_data(self):
        if self.reaction_key not in reaction_database:
            raise ValueError(f"{self.reaction_key} not found in reactions dictionary")
        return reaction_database[self.reaction_key]

    @property
    def stoich(self):
        return get_solid_stoichiometry(self.reaction_data())
    
    @property
    def analytic(self):
        return self.reaction_data()['analytic']
    
    @property
    def eq_num(self):
        return self.reaction_data()['eq_num']

    


# Example usage
if __name__ == '__main__':
    reaction = Reaction("CaCl2·6H2O(s)")

    print(type(reaction.stoich))
    print(reaction.analytic)
