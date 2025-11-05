import sqlite3
import json
from pathlib import Path

current_dir = Path(__file__).parent
input_path = current_dir / "pypitzer_parameter.json"

data = json.load(open(input_path))

output_path = current_dir / "pypitzer.sqlite"

conn = sqlite3.connect(output_path)
cur = conn.cursor()

def coeffs(entry, key):
    d = entry.get(key, {})
    return ",".join(str(d.get(f"a{i}", 0)) for i in range(6))

cur.executescript("""
CREATE TABLE IF NOT EXISTS binary_params (
    pair TEXT PRIMARY KEY,
    ref TEXT,
    T_min REAL,
    T_max REAL,
    eq_num INTEGER,
    b0 TEXT,
    b1 TEXT,
    b2 TEXT,
    c_phi TEXT,
    lambda TEXT,
    theta TEXT
);

CREATE TABLE IF NOT EXISTS ternary_params (
    triplet TEXT PRIMARY KEY,
    ref TEXT,
    T_min REAL,
    T_max REAL,
    eq_num INTEGER,
    psi TEXT,
    zeta TEXT
);
""")

# Insert binary
for pair, entry in data["binary"].items():
    cur.execute("""
        INSERT OR REPLACE INTO binary_params
        (pair, ref, T_min, T_max, eq_num, b0, b1, b2, c_phi, lambda, theta)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        pair,
        entry["ref"],
        entry["T_range"][0],
        entry["T_range"][1],
        entry["eq_num"],
        coeffs(entry, "b0"),
        coeffs(entry, "b1"),
        coeffs(entry, "b2"),
        coeffs(entry, "c_phi"),
        coeffs(entry, "lambda"),
        coeffs(entry, "theta"),
    ))

# Insert ternary
for triplet, entry in data["ternary"].items():
    cur.execute("""
        INSERT OR REPLACE INTO ternary_params
        (triplet, ref, T_min, T_max, eq_num, psi, zeta)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, (
        triplet,
        entry["ref"],
        entry["T_range"][0],
        entry["T_range"][1],
        entry["eq_num"],
        coeffs(entry, "psi"),
        coeffs(entry, "zeta"),
    ))

conn.commit()
conn.close()
