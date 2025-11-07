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

def insert_binary_param(conn, pair, ref, T_range, eq_num,
                        b0=None, b1=None, b2=None, c_phi=None, lambda_=None, theta=None):
    """
    Insert a single binary parameter into the binary_params table.
    
    Parameters:
        conn: sqlite3.Connection
        pair: str, e.g., "Cl-,Na+"
        ref: str, reference
        T_range: list or tuple of [T_min, T_max]
        eq_num: int, equation index
        b0, b1, b2, c_phi, lambda_, theta: str or None, comma-separated 6 coefficients
            If None, defaults to "0,0,0,0,0,0"
    """
    # default for missing coefficient blocks
    default_coeff = "0,0,0,0,0,0"
    pair = preprocess_pair(pair)  # normalize key


    cur = conn.cursor()
    cur.execute("""
        INSERT OR IGNORE INTO binary_params
        (pair, ref, T_min, T_max, eq_num, b0, b1, b2, c_phi, lambda, theta)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        pair,
        ref,
        T_range[0],
        T_range[1],
        eq_num,
        preprocess_data(b0) if b0 is not None else default_coeff,
        preprocess_data(b1) if b1 is not None else default_coeff,
        preprocess_data(b2) if b2 is not None else default_coeff,
        preprocess_data(c_phi) if c_phi is not None else default_coeff,
        preprocess_data(lambda_) if lambda_ is not None else default_coeff,
        preprocess_data(theta) if theta is not None else default_coeff,
    ))
    conn.commit()
    print(f"Inserted/ignored binary parameter for '{pair}'")







if __name__ == "__main__":
    db_path = current_dir / "pypitzer.sqlite"
    conn = sqlite3.connect(db_path)

    
    data = {
    ('Li+,LiCl'):"0,0,0,0,0,0,0,0",
    ('Cl-,LiCl'):"0,0,0,0,0,0,0,0",
    ('K+,LiCl'): "-2.879993582,-1.996147e-3,0,0,0,0.590495613,0.016767927,0",
    ('Na+,LiCl'):"-0.112759482,0.0394693e-3,0,0,0,0,0,0",
}

    for k, v in data.items():
        insert_binary_param(
            conn, 
            pair=k, # test
            ref="Lassin et al., 2015", 
            T_range=[25,300], 
            eq_num=2,
            b0=None, 
            b1=None, 
            b2=None, 
            c_phi=None, 
            lambda_=v, 
            theta=None
        )

