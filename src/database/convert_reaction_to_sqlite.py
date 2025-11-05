import sqlite3
import json
from pathlib import Path
import hashlib

current_dir = Path(__file__).parent
json_path = current_dir / "pypitzer_reaction.json"

db_path = current_dir / "pypitzer.sqlite"

def preprocess_reaction(reaction: str) -> str:
    """
    Clean reaction string:
    - remove extra spaces
    - sort reactants and products alphabetically
    - do NOT split on '+' inside chemical formulas
    """
    lhs, rhs = reaction.split("=")
    
    # split by ' + ' with spaces
    lhs_list = [s.strip() for s in lhs.split(" + ")]
    rhs_list = [s.strip() for s in rhs.split(" + ")]
    
    lhs_sorted = " + ".join(sorted(lhs_list))
    rhs_sorted = " + ".join(sorted(rhs_list))
    
    return f"{lhs_sorted} = {rhs_sorted}"

def reaction_hash(reaction: str) -> str:
    """
    Compute SHA256 hash of the canonicalized reaction string.
    """
    return hashlib.sha256(reaction.encode()).hexdigest()

def json_to_reaction_db(json_path, db_path):
    # Load JSON
    data = json.load(open(json_path))

    # Connect to SQLite
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Create table
    cur.execute("""
    CREATE TABLE IF NOT EXISTS reaction (
        reaction_id TEXT PRIMARY KEY,
        reaction TEXT,
        analytic TEXT,
        eq_num INTEGER,
        ref TEXT,
        T_min REAL,
        T_max REAL
    )
    """)

    # Insert each reaction
    for key, entry in data.items():
        reaction_str = entry["reaction"]
        reaction_str = preprocess_reaction(reaction_str)
        print(reaction_str)
        r_hash = reaction_hash(reaction_str)

        # convert analytic dict to comma-separated string
        analytic_str = ",".join(str(entry["analytic"].get(f"a{i}", 0)) for i in range(6))

        cur.execute("""
            INSERT OR IGNORE INTO reaction
            (reaction_id, reaction, analytic, eq_num, ref, T_min, T_max)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (
            r_hash,
            reaction_str,
            analytic_str,
            entry.get("eq_num", 0),
            entry.get("ref", ""),
            entry.get("T_range", [None, None])[0],
            entry.get("T_range", [None, None])[1]
        ))

    conn.commit()
    conn.close()

if __name__ =="__main__":
    json_to_reaction_db(json_path, db_path)