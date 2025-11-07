import sqlite3
import json
from pathlib import Path
current_dir = Path(__file__).parent

def preprocess_pair(pair: str) -> str:
    """
    Preprocess a pair string:
    1. Remove spaces
    2. Order ions alphabetically
    """
    ions = [ion.strip() for ion in pair.split(",")]  # remove spaces
    ions.sort()  # alphabetical order
    return ",".join(ions)

def preprocess_data(parmeters: str) -> str:
    """
    Preprocess a pair string:
    1. Remove spaces
    2. Order ions alphabetically
    """
    cleaned_parameters = parmeters.strip() # remove spaces
    return cleaned_parameters 

def preprocess_triplet(triplet: str) -> str:
    """
    Preprocess a ternary triplet string:
    1. Remove spaces
    2. Order ions alphabetically
    """
    ions = [ion.strip() for ion in triplet.split(",")]
    ions.sort()
    return ",".join(ions)


def insert_ternary_param(conn, triplet, ref, T_range, eq_num, psi=None, zeta=None):
    """
    Insert a single ternary parameter into the ternary_params table.

    Parameters:
        conn: sqlite3.Connection
        triplet: str, e.g., "Cl-,K+,Na+"
        ref: str, reference
        T_range: list or tuple of [T_min, T_max]
        eq_num: int, equation index
        psi, zeta: str or None, comma-separated 6 coefficients
            If None, defaults to "0,0,0,0,0,0"
    """
    default_coeff = "0,0,0,0,0,0"
    triplet = preprocess_triplet(triplet)

    cur = conn.cursor()
    cur.execute("""
        INSERT OR IGNORE INTO ternary_params
        (triplet, ref, T_min, T_max, eq_num, psi, zeta)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, (
        triplet,
        ref,
        T_range[0],
        T_range[1],
        eq_num,
        preprocess_data(psi) if psi is not None else default_coeff,
        preprocess_data(zeta) if zeta is not None else default_coeff,
    ))
    conn.commit()
    print(f"Successfully inserted/ignored ternary parameter for '{triplet}'")


if __name__ == "__main__":
    db_path = current_dir / "pypitzer.sqlite"
    conn = sqlite3.connect(db_path)




    data = {
        "CO2,Fe+2,Cl-":     "-1.34260256e3,-7.72286e-1,3.91603e-4,0,2.772680974e4,2.5362319406e2,0",
        "CO2,Fe+2,SO4-2":   "-7.37424392e3,-4.608331,2.489207e-3,0,1.431626076e5,1.412302898e5,0",
        "O2,Na+,Cl-":       "0,0,0,0,-2.739,0,0",
        "O2,Na+,OH-":       "-1.25e-2,0,0,0,0,0,0",
        "O2,Na+,NO3-":      "-1.20e-2,0,0,0,0,0,0",
        "O2,Na+,HCO3-":     "0.0,0,0,0,0,0,0",
        "O2,Na+,CO3-2":     "-1.81e-2,0,0,0,0,0,0",
        "O2,Na+,SO4-2":     "-4.60e-2,0,0,0,0,0,0",
        "O2,K+,Cl-":        "-2.11e-2,0,0,0,0,0,0",
        "O2,K+,OH-":        "2.342e-3,0,0,0,0,0,-8.3615e2",
        "O2,K+,NO3-":       "-2.81e-2,0,0,0,0,0,0",
        "O2,K+,SO4-2":      "0.0,0,0,0,0,0,0",
        "O2,Mg+2,Cl-":      "-5.65e-3,0,0,0,0,0,0",
        "O2,Mg+2,SO4-2":    "0,0,0,0,0,0,0",
        "O2,Ca+2,Cl-":      "-1.69e-2,0,0,0,0,0,0",
        "O2,Ca+2,NO3-":     "0,0,0,0,0,0,0",
        "O2,Fe+2,Cl-":      "-5.65e-3,0,0,0,0,0,0",
        "O2,Fe+2,SO4-2":    "0,0,0,0,0,0,0",
        "O2,H+,Cl-":        "-7.7e-3,0,0,0,0,0,0",
    }

    for k, v in data.items():
        insert_ternary_param(
            conn,
            triplet=k, 
            ref="Marion., 2003",
            T_range=[-90, 25],
            eq_num=1,
            psi=None,
            zeta=v
        )

    conn.close()