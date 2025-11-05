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

    insert_ternary_param(
        conn,
        triplet="FeOH+, Fe+2, Cl-",  # test
        ref="Liu et al., 2024",
        T_range=[-90, 25],
        eq_num=1,
        psi="2.8e-2,0,0,0,0,0,0",
        zeta=None
    )

    conn.close()
