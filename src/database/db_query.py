import re
import sqlite3
from pathlib import Path
import hashlib

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


def get_solid_stoichiometry(reaction_str):
    """
    Return only the RHS species dictionary (no outer solid key).
    """
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

class NoDataError(Exception):
    """Raised when a database query returns no results for the given temperature or key."""
    pass

class PypitzerDB:
    def __init__(self, db_path=None):
        if db_path is None:
            current_dir = Path(__file__).parent
            db_path = current_dir / "pypitzer.sqlite"
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row

    # ------------------- existing methods -------------------
    @staticmethod
    def preprocess_pair(pair: str) -> str:
        ions = [ion.strip() for ion in pair.split(",")]
        ions.sort()
        return ",".join(ions)

    @staticmethod
    def str_to_floats(s: str):
        if not s:
            return []
        return [float(x) for x in s.split(",")]

    def get_binary(self, pair: str, T: float):
        pair = self.preprocess_pair(pair)
        cur = self.conn.cursor()
        cur.execute("""
            SELECT *
            FROM binary_params
            WHERE pair = ?
            AND T_min < ?
            AND T_max >= ?
        """, (pair, T, T))
        row = cur.fetchone()
        if not row:
            raise NoDataError(f"No data found for pair {pair} at T={T} °C")

        row_dict = dict(row)
        for key in ["b0", "b1", "b2", "c_phi", "lambda", "theta"]:
            row_dict[key] = self.str_to_floats(row_dict.get(key))
        return row_dict


    def get_ternary(self, triplet: str, T: float):
        triplet = self.preprocess_pair(triplet)
        cur = self.conn.cursor()
        cur.execute("""
            SELECT *
            FROM ternary_params
            WHERE triplet = ?
            AND T_min < ?
            AND T_max >= ?
        """, (triplet, T, T))
        row = cur.fetchone()
        if not row:
            raise NoDataError(f"No data found for triplet {triplet} at T={T} °C")

        row_dict = dict(row)
        for key in ["psi", "zeta"]:
            row_dict[key] = self.str_to_floats(row_dict.get(key))
        return row_dict

    # ------------------- new reaction method -------------------
    @staticmethod
    def preprocess_reaction(reaction: str) -> str:
        """Clean reaction string: remove extra spaces, preserve ions, sort reactants/products."""
        lhs, rhs = reaction.split("=")
        lhs_list = [s.strip() for s in lhs.split(" + ")]
        rhs_list = [s.strip() for s in rhs.split(" + ")]
        lhs_sorted = " + ".join(sorted(lhs_list))
        rhs_sorted = " + ".join(sorted(rhs_list))
        return f"{lhs_sorted} = {rhs_sorted}"

    @staticmethod
    def reaction_hash(reaction: str) -> str:
        """Compute SHA256 hash of canonical reaction string."""
        return hashlib.sha256(reaction.encode()).hexdigest()

    def get_reaction(self, reaction: str, T: float):
        reaction_canon = self.preprocess_reaction(reaction)
        r_hash = self.reaction_hash(reaction_canon)
        cur = self.conn.cursor()
        cur.execute("""
            SELECT *
            FROM reaction
            WHERE reaction_id = ?
            AND T_min < ?
            AND T_max >= ?
        """, (r_hash, T, T))
        row = cur.fetchone()
        if not row:
            raise NoDataError(f"No reaction data found for {reaction} at T={T} °C")

        row_dict = dict(row)
        row_dict["analytic"] = self.str_to_floats(row_dict.get("analytic")) # type: ignore
        row_dict["stoich"] = get_solid_stoichiometry(row_dict.get("reaction"))
        return row_dict

    def close(self):
        self.conn.close()


# ------------------- usage example -------------------
if __name__ == "__main__":
    db = PypitzerDB()

    # binary = db.get_binary("Na+,Cl-")
    # ternary = db.get_ternary("Na+,Cl-,K+")
    reaction = db.get_reaction("NaCl-2H2O(s) = Na+(aq) + Cl-(aq) + 2H2O(l)", 25)

    # print("Binary:", binary)
    # print("Ternary:", ternary)
    print("Reaction:", reaction)

    db.close()
