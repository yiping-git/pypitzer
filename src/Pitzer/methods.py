# -*- coding: utf-8 -*-
# Author: Yiping Liu
# Description: This script calculates the average of a list of numbers.
# Version: 1.0
# Last Modified: Nov 11, 2025
import json
import re
import numpy as np
import pandas as pd
import itertools

# from src.database.spencer_revised_chemical_potential import spencer_chemical_potential_db
# from src.database.marion_chemical_potential import marion_chemical_potential_db
# from src.database.lassin_chemical_potential import lassin_chemical_potential_db
# from src.database.solid_data import solids

from src.public.j_x import compute_j_jp
from src.public.low_level import get_charge_number
from pathlib import Path

here = Path(__file__).resolve().parent.parent

db_json_path = here / "database/pypitzer_parameter.json"
with open(db_json_path,"r",encoding="utf-8") as f:
    parameter_db = json.load(f)

# Query methods from json db.
def binary_query(ion_pair):
    binary_pair = ",".join(sorted(list(ion_pair)))
    parameters = parameter_db['binary'][binary_pair]
    return parameters

def ternary_query(ion_pair):
    ternary_pair = ",".join(sorted(list(ion_pair)))
    parameters = parameter_db['ternary'][ternary_pair]
    return parameters


def chemical_potential_lassin(a1, a2, a3, a4, a5, t):
    log_k = a1 + a2 * t + a3 / t + a4 * np.log10(t) + a5 / (t ** 2)
    return log_k * 2.303


def parameter_cal_lassin(a1, a2, a3, a4, a5, a6, a7, a8, t):
    p = a1 + a2 * t + a3 * t ** 2 + a4 * t ** 3 + a5 / t + a6 * np.log(t) + a7 / (t - 263) + a8 / (680 - t)
    return p

def parameter_cal_spencer(a1, a2, a3, a4, a5, a6, t):
    parameter = a1 + a2 * t + a3 * t ** 2 + a4 * t ** 3 + a5 / t + a6 * np.log(t)
    return parameter

def parameter_fun_spencer(a, T):
    return a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] / T + a[5] * np.log(T)


def parameter_cal_moller(a1, a2, a3, a4, a5, a6, a7, a8, t):
    p = a1 + a2 * t + a3 / t + a4 * np.log(t) + a5 / (t - 263) + a6 * t ** 2 + a7 / (680 - t) + a8 / (t - 227)
    return p

def parameter_cal_holmes(a1, a2, a3, a4, a5, a6, t):
    t_r = 298.15
    parameter = a1 + a2 * (1 / t - 1 / t_r) + a3 * np.log(t / t_r) + a4 * (t - t_r) + a5 * (t ** 2 - t_r ** 2) + a6 * np.log(
        t - 260)
    return parameter

def parameter_cal_marion(a1, a2, a3, a4, a5, a6, a7, t):
    # reference: [4]
    parameter = a1 + a2 * t + a3 * t ** 2 + a4 * t ** 3 + a5 / t + a6 * np.log(t) + a7 / (t ** 2)
    return parameter



def a_phi_moller(t):
    a1 = 3.36901532e-01
    a2 = -6.32100430e-04
    a3 = 9.14252359e00
    a4 = -1.35143986e-02
    a5 = 2.26089488e-03
    a6 = 1.92118597e-06
    a7 = 4.52586464e01
    a8 = 0
    return a1 + a2 * t + a3 / t + a4 * np.log(t) + a5 / (t - 263) + a6 * t ** 2 + a7 / (680 - t) + a8 / (t - 227)


def a_phi_spencer(T):
    """
    Calculate the A(φ) (Debye-Hukel constant) with temperature-dependent equation.
    :param T: temperature of solution, [K]
    :return: A(φ)
    """
    a1 = 8.66836498e1
    a2 = 8.48795942e-2
    a3 = -8.88785150e-5
    a4 = 4.88096393e-8
    a5 = -1.32731477e3
    a6 = -1.76460172e1
    a_phi = parameter_cal_spencer(a1, a2, a3, a4, a5, a6, T)
    return a_phi


def g_func(a):
    result = 2 * (
            1 - (1 + a) * np.exp(-a)
    ) / a ** 2
    return result


def g_func_prime(a):
    g_prime = -2 * (
            1 - (1 + a + a ** 2 / 2) * np.exp(-a)
    ) / a ** 2
    return g_prime


def get_beta(parameters, pair, i):
    z1 = get_charge_number(pair[0])
    z2 = get_charge_number(pair[1])

    alpha = 2  # kg^(1/2)⋅mol^(-1/2)
    alpha_1 = 1.4  # kg^(1/2)⋅mol^(-1/2)
    alpha_2 = 12  # kg^(1/2)⋅mol^(-1/2)
    b0, b1, b2 = parameters

    if abs(z1) == 2 and abs(z2) == 2:
        beta_mx = b0 + b1 * g_func(
            alpha_1 * (i ** (1 / 2))
        ) + b2 * g_func(
            alpha_2 * (i ** (1 / 2))
        )
    else:
        beta_mx = b0 + b1 * g_func(
            alpha * i ** (1 / 2)
        )
    return beta_mx


def get_beta_prime(parameters, pair, i):
    charge_number1 = get_charge_number(pair[0])
    charge_number2 = get_charge_number(pair[1])

    alpha = 2  # kg^(1/2)⋅mol^(-1/2)
    alpha_1 = 1.4  # kg^(1/2)⋅mol^(-1/2)
    alpha_2 = 12  # kg^(1/2)⋅mol^(-1/2)
    b0, b1, b2 = parameters
    if abs(charge_number1) == 2 and abs(charge_number2) == 2:
        beta_prime = (
                             b1 * g_func_prime(
                         alpha_1 * (i ** (1 / 2))
                     ) + b2 * g_func_prime(
                         alpha_2 * (i ** (1 / 2))
                     )
                     ) / i
    else:
        beta_prime = b1 * g_func_prime(alpha * i ** (1 / 2)) / i
    return beta_prime


def get_beta_phi(parameters, pair, i):
    z1 = get_charge_number(pair[0])
    z2 = get_charge_number(pair[1])

    alpha = 2  # kg^(1/2)⋅mol^(-1/2)
    alpha_1 = 1.4  # kg^(1/2)⋅mol^(-1/2)
    alpha_2 = 12  # kg^(1/2)⋅mol^(-1/2)

    b0, b1, b2 = parameters

    """for 2-2 type salts"""
    if abs(z1) == 2 and abs(z2) == 2:
        beta_phi = b0 + b1 * np.exp(-alpha_1 * (i ** (1 / 2))) + b2 * np.exp(-alpha_2 * (i ** (1 / 2)))
    else:
        """for 1-1, 1-2, 2-1, 3-1, 4-1 type salts"""
        beta_phi = b0 + b1 * np.exp(-alpha * (i ** (1 / 2)))
    return beta_phi

def get_c(c_phi, z_m, z_x):
    c_mx = c_phi / (2 * (abs(z_m * z_x)) ** 0.5)
    return c_mx

def get_c_gamma(c_phi):
    c_gamma = 3 * c_phi / 2
    return c_gamma


def get_f(a_phi, i):
    """
    :param a_phi: A_phi (Debye-Hukel constant)
    :param i: ionic strength
    :return: expression of "f" function
    """
    b = 1.2
    f = - (4 * i * a_phi / b) * np.log(1 + b * i ** (1 / 2))
    return f


def get_f_gamma(a_phi, i):
    """
    :param a_phi:A_phi (Debye-Hukel constant)
    :param i: ionic strength
    :return: the "f^gamma" function in Pitzer's model
    """
    b = 1.2
    f_gamma = - a_phi * (
            i ** (1 / 2) / (1 + b * i ** (1 / 2)) + (2 / b) * np.log(
        1 + b * i ** (1 / 2))
    )
    return f_gamma


"""
*** End of the Pitzer parameters dealing *** 
"""


def get_x_mn(z_m: int, z_n: int, a_phi: float, i: float):
    """
    calculate the 'x' value of ions 'm' and 'n'.
    :param z_m: charge number of ion 'm'
    :param z_n: charge number of ion 'n'
    :param a_phi: Avogedral's number of this solution
    :param i: Ionic strength of this solution
    :return: 'x_mn' for calculating the 'J' value
    :reference: [2] p9
    """
    x_mn = 6 * z_m * z_n * a_phi * i ** 0.5
    return x_mn


def get_e_theta(z_m, z_n, a_phi, i):
    """
    :param z_m: charge number of species m
    :param z_n: charge number of species n
    :param a_phi:
    :param i: ionic strength
    :return: e_theta and e_theta_prime
    :reference: [1] p123
    """

    x_mn = get_x_mn(z_m, z_n, a_phi, i)
    x_mm = get_x_mn(z_m, z_m, a_phi, i)
    x_nn = get_x_mn(z_n, z_n, a_phi, i)

    mn = compute_j_jp(x_mn)
    mm = compute_j_jp(x_mm)
    nn = compute_j_jp(x_nn)
    j_mn = mn['j_x']
    j_mn_prime = mn['j_xp']
    j_mm = mm['j_x']
    j_mm_prime = mm['j_xp']
    j_nn = nn['j_x']
    j_nn_prime = nn['j_xp']

    e_theta = (z_m * z_n / (4 * i)) * (j_mn - 0.5 * j_mm - 0.5 * j_nn)
    e_theta_prime = -(e_theta / i) + (z_m * z_n / (8 * i ** 2)) * (
            x_mn * j_mn_prime - 0.5 * x_mm * j_mm_prime - 0.5 * x_nn * j_nn_prime)
    return {
        "e_theta": e_theta,
        "e_theta_prime": e_theta_prime,
    }


def calculate_ionic_strength(molalities):
    data = molalities
    ions = data.keys()
    sum_value = 0
    for ion in ions:
        charge_number = get_charge_number(ion)
        sum_value += data[ion] * (charge_number ** 2)
    return sum_value / 2


def calculate_molality(x, species):
    x1, x2 = x
    molalities = {}
    for key, value in species.items():
        if key != 'Cl-':
            molalities[key] = value * x1
    molalities['Cl-'] = x2
    return molalities



def calculate_charge_balance(x, molalities):
    balance = 0
    for species in molalities.keys():
        balance += get_charge_number(species) * molalities[species]
    return balance


def get_chemical_potential(parameters, T, eq_num):
    if eq_num == 0:
        return parameter_fun_spencer(parameters, T)
    else:
        return 0


def group_components(components):
    """
    Find groups from components of ions and neutral species
    :param components: consists of cations, anions and neutral species.
    :return: groups
    """
    cations = [c for c in components if '+' in c]
    anions = [a for a in components if '-' in a]
    neutrals = [n for n in components if '+' not in n and '-' not in n]

    cation_anion_pairs = list(itertools.product(cations, anions))
    cation_pairs = list(itertools.combinations(cations, 2)) if len(cations) >= 2 else []
    anion_pairs = list(itertools.combinations(anions, 2)) if len(anions) >= 2 else []
    neutral_pairs = list(itertools.combinations(neutrals, 2))

    neutral_ion_pairs = list(itertools.product(neutrals, cations + anions))

    neutral_cation_anion_pairs = [(a, *b) for a in neutrals for b in cation_anion_pairs]

    return {
        'cations': cations,
        'anions': anions,
        'neutrals': neutrals,
        'cation_anion_pairs': cation_anion_pairs,
        'cation_pairs': cation_pairs,
        'anion_pairs': anion_pairs,
        'neutral_pairs': neutral_pairs,
        'neutral_ion_pairs': neutral_ion_pairs,
        'neutral_cation_anion_pairs': neutral_cation_anion_pairs
    }



def get_parameter_pypitzer(parameters, T, eq_num):
    if eq_num == 0:
        return parameter_fun_spencer(parameters, T)
    return parameter_fun_spencer(parameters, T)


def ternary_parameter_cal(pair, T):
    parameters = ternary_query(pair)

    para_psi  = list(parameters['psi'].values())
    para_zeta = list(parameters['zeta'].values()) if 'zeta' in parameters else None
    eq_num = parameters['eq_num']

    psi = get_parameter_pypitzer(para_psi, T, eq_num)
    zeta = get_parameter_pypitzer(para_zeta, T, eq_num) if para_zeta else 0
    return (psi, zeta)


def get_beta_012(ion_pair, T):
    parameters = binary_query(ion_pair)

    para_b0 = list(parameters['b0'].values())
    para_b1 = list(parameters['b1'].values())
    para_b2 = list(parameters['b2'].values()) if "b2" in parameters else None
    eq_num = parameters['eq_num']

    b0 = get_parameter_pypitzer(para_b0, T, eq_num) 
    b1 = get_parameter_pypitzer(para_b1, T, eq_num)
    b2 = get_parameter_pypitzer(para_b2, T, eq_num) if para_b2 else 0
    return (b0, b1, b2)


def beta_calculate(ion_pair, ionic_strength, T):
    beta_012 = get_beta_012(ion_pair, T)
    b_phi = get_beta(beta_012, ion_pair, ionic_strength)
    return b_phi


def beta_phi_calculate(ion_pair, ionic_strength, T):
    beta_012 = get_beta_012(ion_pair, T)
    b_phi = get_beta_phi(beta_012, ion_pair, ionic_strength)
    return b_phi

def beta_prime_calculate(ion_pair, ionic_strength, T):
    beta_012 = get_beta_012(ion_pair, T)
    b_prime = get_beta_prime(beta_012, ion_pair, ionic_strength)
    return b_prime

def c_calculate(ion_pair, T):
    parameters = binary_query(ion_pair)
    para_cphi = list(parameters['c_phi'].values())

    eq_num = parameters['eq_num']
    c0 = get_parameter_pypitzer(para_cphi, T, eq_num)

    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])

    c = get_c(c0, z1, z2)
    return c

def lambda_cal(pair,T):
    parameters = binary_query(pair)

    para_lambda = list(parameters['lambda'].values())
    eq_num = parameters['eq_num']
    return get_parameter_pypitzer(para_lambda, T, eq_num)


def cc_phi_calculate(ion_pair, a_phi, ionic_strength, T):
    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])
    parameters = binary_query(ion_pair)

    para_theta = list(parameters['theta'].values())
    eq_num = parameters['eq_num']

    theta = get_parameter_pypitzer(para_theta, T, eq_num)

    if z1 != z2:
        e_thetas = get_e_theta(z1, z2, a_phi, ionic_strength)
        e_theta = e_thetas['e_theta']
        e_theta_prime = e_thetas['e_theta_prime']
    else:
        e_theta = 0
        e_theta_prime = 0

    phi = theta + e_theta
    phi_prime = e_theta_prime

    return {
        'phi': phi,
        'phi_prime': phi_prime
    }


def aa_phi_calculate(ion_pair, a_phi, ionic_strength, T):
    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])

    parameters = binary_query(ion_pair)

    para_theta = list(parameters['theta'].values())
    eq_num = parameters['eq_num']

    theta = get_parameter_pypitzer(para_theta, T, eq_num)
    if z1 != z2:
        e_thetas = get_e_theta(z1, z2, a_phi, ionic_strength)
        e_theta = e_thetas['e_theta']
        e_theta_prime = e_thetas['e_theta_prime']
    else:
        e_theta = 0
        e_theta_prime = 0

    phi = theta + e_theta
    phi_prime = e_theta_prime

    return {
        "phi": phi,
        "phi_prime": phi_prime
    }


"""
***************************
 DATA OF SOLID PHASE
***************************
"""
current_dir = Path(__file__).parent
json_file = here / "database/pypitzer_reaction.json"

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
    stoich = {}
    if reaction_str:
        lhs, rhs = reaction_str.split('=', 1)

        rhs_no_states = clean_state_marks(rhs)
        species_tokens = [tok.strip() for tok in re.split(r'\s+\+\s+', rhs_no_states.strip()) if tok.strip()]

        
        for tok in species_tokens:
            name, info = parse_species(tok)
            stoich[name] = info
    return stoich


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



"""
Reference
[1] Pitzer K S. Activity coefficients in electrolyte solutions[M]. CRC press, 2018.
[2] Pitzer K S. Thermodynamics of electrolytes. V. Effects of higher-order electrostatic terms[J]. Journal of Solution 
Chemistry, 1975, 4(3): 249-265.
[3] Bradley, Daniel J., and Kenneth S. Pitzer. "Thermodynamics of electrolytes. 12. Dielectric properties of water and 
Debye-Hueckel parameters to 350. degree. C and 1 kbar." Journal of physical chemistry 83.12 (1979): 1599-1603.
[4] Marion GM, Catling DC, Kargel JS. Modeling aqueous ferrous iron chemistry at low temperatures with application to 
Mars. Geochimica et cosmochimica Acta. 2003 Nov 15;67(22):4251-66.
"""
