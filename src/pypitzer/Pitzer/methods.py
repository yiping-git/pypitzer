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

from pypitzer.public.j_x import compute_j_jp
from pypitzer.public.low_level import get_charge_number


def parameter_fun_spencer(a, T):
    # eq_num = 0
    return a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] / T + a[5] * np.log(T)

def parameter_fun_marion(a, T):
    # eq_num = 1
    return  a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] / T + a[5] * np.log(T) + a[6] / (T ** 2)

def parameter_fun_eq2(a, T):
    # eq_num = 2
    return a[1] + a[2] * T + a[3] * T ** 2 + a[4] * T ** 3 + a[5] / T + a[6] * np.log(T) + a[7 ]/ (T - 263) + a[8] / (680 - T)
    
def parameter_fun_moller(a, T):
    # eq_num = 3
    return a[0] + a[1] * T + a[2] / T + a[3] * np.log(T) + a[4] / (T - 263) + a[5]* T ** 2 + a[6] / (680 - T) + a[7] / (T - 227)

def parameter_fun_holmes(a, T):
    # eq_num = 4
    Tr = 298.15
    return a[1] + a[2]*(1/T - 1/Tr) + a[3]*np.log(T/Tr) + a[4]*(T - Tr) + a[5]*(T**2 - Tr**2) + a[6]*np.log(T - 260)
     
def parameter_fun_eq5(a, T):
    # eq_num = 5
    log_k = a[0] + a[1] * T + a[2] / T + a[3] * np.log10(T) + a[4] / T**2
    return log_k * 2.303


def get_parameter_pypitzer(a, T, eq_num):
    if eq_num == 0:
        return parameter_fun_spencer(a, T)
    elif eq_num == 1:
        return parameter_fun_marion(a, T)
    elif eq_num == 2:
        return parameter_fun_eq2(a, T)
    elif eq_num == 3:
        return parameter_fun_moller(a, T)
    elif eq_num == 4:
        return parameter_fun_holmes(a, T)
    elif eq_num == 5:
        return parameter_fun_eq5(a, T)



def a_phi_moller(T):
    a = [
         3.36901532e-01,
        -6.32100430e-04,
         9.142523590e00,
        -1.35143986e-02,
         2.26089488e-03,
         1.92118597e-06,
         4.525864640e01,
                      0,
    ]
    return parameter_fun_moller(a, T)


def a_phi_spencer(T):
    """
    Calculate the A(φ) (Debye-Hukel constant) with temperature-dependent equation.
    :param T: temperature of solution, [K]
    :return: A(φ)
    """
    a = [
            8.66836498e1,
            8.48795942e-2,
           -8.88785150e-5,
            4.88096393e-8,
           -1.32731477e3,
           -1.76460172e1
        ]
    return parameter_fun_spencer(a, T)


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


import itertools

def order_tuple(t):
    """
    Order elements of a tuple alphabetically.
    Example: ("Na+", "Cl-") -> ("Cl-", "Na+")
    """
    return tuple(sorted([str(x).strip() for x in t]))

def group_components(components):
    """
    Find groups from components of ions and neutral species.
    :param components: list of cations, anions, and neutral species.
    :return: dictionary of grouped components with tuples
    """
    cations = [c for c in components if '+' in c]
    anions = [a for a in components if '-' in a]
    neutrals = [n for n in components if '+' not in n and '-' not in n]

    # Basic pairs
    cation_anion_pairs = [order_tuple(p) for p in itertools.product(cations, anions)]
    cation_pairs = [order_tuple(p) for p in itertools.combinations(cations, 2)] if len(cations) >= 2 else []
    anion_pairs = [order_tuple(p) for p in itertools.combinations(anions, 2)] if len(anions) >= 2 else []
    neutral_pairs = [order_tuple(p) for p in itertools.combinations(neutrals, 2)]
    neutral_ion_pairs = [order_tuple(p) for p in itertools.product(neutrals, cations + anions)]
    neutral_cation_anion_pairs = [order_tuple((n, c, a)) for n in neutrals for c, a in itertools.product(cations, anions)]

    # Additional triplets
    cca = [order_tuple((c1, c2, a)) for c1, c2 in itertools.combinations(cations, 2) for a in anions]
    aac = [order_tuple((a1, a2, c)) for a1, a2 in itertools.combinations(anions, 2) for c in cations]

    # Neutral-cation / neutral-anion pairs
    nc = [order_tuple((n, c)) for n in neutrals for c in cations]
    na = [order_tuple((n, a)) for n in neutrals for a in anions]
    nca = [order_tuple((n, c, a)) for n in neutrals for c, a in itertools.product(cations, anions)]

    return {
        'cations': tuple(cations),
        'anions': tuple(anions),
        'neutrals': tuple(neutrals),
        'cation_anion_pairs': tuple(cation_anion_pairs),
        'cation_pairs': tuple(cation_pairs),
        'anion_pairs': tuple(anion_pairs),
        'nn': tuple(neutral_pairs),
        'neutral_ion_pairs': tuple(neutral_ion_pairs),
        'neutral_cation_anion_pairs': tuple(neutral_cation_anion_pairs),
        'cca': tuple(cca),
        'aac': tuple(aac),
        'nc': tuple(nc),
        'na': tuple(na),
        'nca': tuple(nca)
    }




def ternary_parameter_cal(T, parameters):
    para_psi  = parameters['psi']
    para_zeta = parameters['zeta'] if 'zeta' in parameters else None
    eq_num = parameters['eq_num']

    psi = get_parameter_pypitzer(para_psi, T, eq_num)
    zeta = get_parameter_pypitzer(para_zeta, T, eq_num) if para_zeta else 0
    return (psi, zeta)


def get_beta_012(ion_pair, T, parameters):
    para_b0 = parameters['b0']
    para_b1 = parameters['b1']
    para_b2 = parameters['b2'] if "b2" in parameters else None
    eq_num = parameters['eq_num']

    b0 = get_parameter_pypitzer(para_b0, T, eq_num) 
    b1 = get_parameter_pypitzer(para_b1, T, eq_num)
    b2 = get_parameter_pypitzer(para_b2, T, eq_num) if para_b2 else 0
    return (b0, b1, b2)


def beta_calculate(ion_pair, ionic_strength, T, parameters):
    beta_012 = get_beta_012(ion_pair, T, parameters)
    b_phi = get_beta(beta_012, ion_pair, ionic_strength)
    return b_phi


def beta_phi_calculate(ion_pair, ionic_strength, T, parameters):
    beta_012 = get_beta_012(ion_pair, T, parameters)
    b_phi = get_beta_phi(beta_012, ion_pair, ionic_strength)
    return b_phi

def beta_prime_calculate(ion_pair, ionic_strength, T, parameters):
    beta_012 = get_beta_012(ion_pair, T, parameters)
    b_prime = get_beta_prime(beta_012, ion_pair, ionic_strength)
    return b_prime

def c_calculate(ion_pair, T, parameters):
    para_cphi = parameters['c_phi']

    eq_num = parameters['eq_num']
    c0 = get_parameter_pypitzer(para_cphi, T, eq_num)

    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])

    c = get_c(c0, z1, z2)
    return c

def lambda_cal(pair,T, parameters):
    # parameters = binary_query(pair)

    para_lambda = parameters['lambda']
    eq_num = parameters['eq_num']
    return get_parameter_pypitzer(para_lambda, T, eq_num)


def zeta_cal(pair,T, parameters):
    # parameters = binary_query(pair)
    para_zeta = parameters['zeta']
    eq_num = parameters['eq_num']
    return get_parameter_pypitzer(para_zeta, T, eq_num)

def cc_phi_calculate(ion_pair, a_phi, ionic_strength, T, parameters):
    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])

    para_theta = parameters['theta']
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


def aa_phi_calculate(ion_pair, a_phi, ionic_strength, T, parameters):
    z1 = get_charge_number(ion_pair[0])
    z2 = get_charge_number(ion_pair[1])

    # parameters = binary_query(ion_pair)

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


# class Reaction:
#     def __init__(self, reaction_key):
#         self.reaction_key = reaction_key

#     def reaction_data(self):
#         if self.reaction_key not in reaction_database:
#             raise ValueError(f"{self.reaction_key} not found in reactions dictionary")
#         return reaction_database[self.reaction_key]

#     @property
#     def stoich(self):
#         return get_solid_stoichiometry(self.reaction_data())
    
#     @property
#     def analytic(self):
#         return self.reaction_data()['analytic']
    
#     @property
#     def eq_num(self):
#         return self.reaction_data()['eq_num']



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
