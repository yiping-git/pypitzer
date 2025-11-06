import sqlite3
from pathlib import Path
import hashlib

current_dir = Path(__file__).parent
json_path = current_dir / "pypitzer_reaction.json"

db_path = current_dir / "pypitzer.sqlite"

# --- canonicalize reaction string ---
def preprocess_reaction(reaction: str) -> str:
    lhs, rhs = reaction.split("=")
    lhs_list = [s.strip() for s in lhs.split(" + ")]
    rhs_list = [s.strip() for s in rhs.split(" + ")]
    lhs_sorted = " + ".join(sorted(lhs_list))
    rhs_sorted = " + ".join(sorted(rhs_list))
    return f"{lhs_sorted} = {rhs_sorted}"

def reaction_hash(reaction: str) -> str:
    return hashlib.sha256(reaction.encode()).hexdigest()
    
def insert_reaction(conn: sqlite3.Connection,
                    reaction: str,
                    analytic: list[float],
                    eq_num: int,
                    ref: str,
                    T_min: float,
                    T_max: float):
    """
    Insert a single reaction into the 'reaction' table.

    Parameters
    ----------
    conn : sqlite3.Connection
        Open connection to the SQLite database.
    reaction : str
        Reaction string, e.g., "H2O = H+ + OH-".
    analytic : list of float
        Coefficients for the analytic lnK expression.
    eq_num : int
        Equation number / model identifier.
    ref : str
        Reference for the data.
    T_min : float
        Minimum temperature (inclusive).
    T_max : float
        Maximum temperature (exclusive).
    """
    

    reaction_canon = preprocess_reaction(reaction)
    r_hash = reaction_hash(reaction_canon)
    
    # Convert analytic list to comma-separated string
    analytic_str = ",".join(str(x) for x in analytic)
    
    cur = conn.cursor()
    try:
        cur.execute("""
            INSERT INTO reaction (reaction_id, reaction, analytic, eq_num, ref, T_min, T_max)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (r_hash, reaction_canon, analytic_str, eq_num, ref, T_min, T_max))
        conn.commit()
    except sqlite3.IntegrityError as e:
        raise ValueError(f"Failed to insert reaction: {e}")



if __name__ == "__main__":
    conn = sqlite3.connect(db_path)

    data =[
            {
            "reaction": "LiCl(aq) = Li+(aq) + Cl-(aq)",
            "analytic": [-3.5869,0.007566,-1254.4,0,0,],
            "T_min": 25,
            "T_max": 300,
            },
            {
            "reaction": "LiCl-H2O(s) = Li+(aq) + Cl-(aq) + H2O(l)",
            "analytic": [36.2702,0.005655998,2.696341,-13.2627354,0,],
            "T_min": 25,
            "T_max": 300,
            },
            {
            "reaction": "LiCl-2H2O(s) = Li+(aq) + Cl-(aq) + 2H2O(l)",
            "analytic":[6.1532,-0.0064179,0,0,0,],
            "T_min": 25,
            "T_max": 300,
            },
            {
            "reaction": "KCl(s) = K+(aq) + Cl-(aq)",
            "analytic": [6.496259873, -0.012323368, -1537.997915, 1.507740146, -42632.32102],
            "T_min": 25,
            "T_max": 300,
            },
            {
            "reaction": "NaCl(s) = Na+(aq) + Cl-(aq)",
            "analytic": [-752.24954, -0.11904958, 41385.703, 274.17933, -2480.9109e+03],
            "T_min": 25,
            "T_max": 300,
            },
            {
            "reaction": "LiOH-H2O(s) = Li+(aq) + OH-(aq) + H2O(l)",
            "analytic": [-3.8702256, -0.00914328, 0.4129607, 3.0072101, 0],
            "T_min": 25,
            "T_max": 200,
            },
            {
            "reaction": "LiOH(s) = Li+(aq) + OH-(aq)",
            "analytic": [-256.08895, -0.06047962, 9290.232715, 98.8913449, 142.35284],
            "T_min": 50,
            "T_max": 200,
            },
    ]

    for d in data:
        insert_reaction(
            conn=conn,
            reaction= d['reaction'],
            analytic= d['analytic'],
            eq_num= 5,
            ref= "Lassin et al., 2015",
            T_min= d['T_min'],
            T_max= d['T_max'],
        )
