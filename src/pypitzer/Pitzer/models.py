# -*- coding: utf-8 -*-
# Author: Yiping Liu
# Description:
# Version: 1.0
# Last Modified: Nov 03, 2025

import time

import numpy as np
from scipy.optimize import minimize

from pypitzer.database.db_query import PypitzerDB
from pypitzer.public.low_level import type_of_species, get_charge_number
from pypitzer.public.icemelting import (clegg_and_brimblecombe, spencer, monnin)

import pypitzer.Pitzer.methods as pm

db = PypitzerDB()


class FluidPitzer:
    def __init__(self, x0, species, equilibrium, t=25.0, neutral=False, ):
        """
        Initiate the solution.
        :param x0: initiate values of x(x1, x2), namely value of mNa and mCl, [mol/kg].
        :param species: a dict containing aqueous spcies and their relative concentrations to Na, [-].
        :param t: melting temperature, [°C].
        :param equilibrium: str, a chemical expression indicating the solubility equilibrium, e.g., "NaCl(s) = Na+(aq) + Cl-(aq)".
        """
        self.x0 = x0
        self.t = t
        self.T = t + 273.16
        self.input_species = species
        self.equilibrium = equilibrium
        self.neutral = neutral  # to determine if a neutral species is considered in this fluid.
        self.solid_disso = db.get_reaction(equilibrium, self.t)

        self._species = None
        self._component_groups = None
        self._param_db = None

    @property
    def species(self):
        if self._species is None :
            species = self.input_species
            species['Cl-'] = 0 # add Cl into the dict and set a default value of 0.
            self._species = species
        return  self._species
    
    @property
    def component_groups(self):
        if self._component_groups is None:
            components = tuple(self.species.keys())
            self._component_groups = pm.group_components(components)
        return self._component_groups
    
    @staticmethod
    def order_str(tuple_pair):
        """
        Order elements of a tuple alphabetically and join them as a string.
        Example: ("Na+", "Cl-") -> "Cl-,Na+"
        """
        sorted_list = sorted([str(x).strip() for x in tuple_pair])
        return ",".join(sorted_list)

    @property
    def param_db(self):
        result = {}
        if self._param_db is None:
            if self.component_groups['cation_anion_pairs']:
                for ca in self.component_groups['cation_anion_pairs']:
                    ca_str = self.order_str(ca)
                    result[ca_str] = db.get_binary(ca_str, self.t)
            if self.component_groups['cation_pairs']:
                for cc in self.component_groups['cation_pairs']:
                    cc_str = self.order_str(cc)
                    result[cc_str] = db.get_binary(cc_str, self.t)
            if self.component_groups['anion_pairs']:
                for aa in self.component_groups['anion_pairs']:
                    aa_str = self.order_str(aa)
                    result[aa_str] = db.get_binary(aa_str, self.t)
            if self.component_groups['nc']:
                for pair in self.component_groups['nc']:
                    pair_str = self.order_str(pair)
                    result[pair_str] = db.get_binary(pair_str, self.t)
            if self.component_groups['na']:
                for pair in self.component_groups['na']:
                    pair_str = self.order_str(pair)
                    result[pair_str] = db.get_binary(pair_str, self.t)
            if self.component_groups['nn']:
                for pair in self.component_groups['nn']:
                    pair_str = self.order_str(pair)
                    result[pair_str] = db.get_binary(pair_str, self.t)
            if self.component_groups['nca']:
                for pair in self.component_groups['nca']:
                    pair_str = self.order_str(pair)
                    result[pair_str ] = db.get_ternary(pair_str, self.t)
            if self.component_groups['cca']:
                for pair in self.component_groups['cca']:
                    pair_str = self.order_str(pair)
                    result[pair_str] = db.get_ternary(pair_str, self.t)
            if self.component_groups['aac']:
                for pair in self.component_groups['aac']:
                    pair_str = self.order_str(pair)
                    result[pair_str] = db.get_ternary(pair_str, self.t)
            self._param_db = result 
        return self._param_db 


    def get_molalities(self, x):
        """
        x = (x1,x2)
        x1: molality of Na
        x2: molality of Cl
        """
        molalities = pm.calculate_molality(tuple(x), self.species)
        return molalities

    def charge_balance(self, x):
        # x = tuple(x)
        # balance = pm.calculate_charge_balance(x, self.get_molalities(x))
        # return balance
        balance = 0
        molalities = self.get_molalities(x)
        for species in molalities.keys():
            balance += get_charge_number(species) * molalities[species]
        return balance

    def get_ionic_strength(self, x):
        """
        For calculating ionic strength. I = (1/2) ∑ mᵢ·zᵢ²
        :param x: a tuple (x1,x2).
        :return: ionic strength.
        """
        molalities = self.get_molalities(x)
        sum_value = 0

        for ion in molalities.keys():
            charge_number = get_charge_number(ion)
            sum_value += molalities[ion] * charge_number ** 2
        return sum_value / 2
    

    def get_z(self, x):
        molalities = self.get_molalities(x)
        ions = molalities.keys()
        z_value = 0
        for ion in ions:
            charge_number = get_charge_number(ion)
            z_value += molalities[ion] * abs(charge_number)
        return z_value


    def get_a_phi(self):
        return pm.a_phi_spencer(self.T)

    def get_b(self, x):
        pair_parameters = {}
        cation_anion_pairs = self.component_groups['cation_anion_pairs']
        for cap in cation_anion_pairs:
            parameters = self.param_db[self.order_str(cap)]
            pair_parameters[cap] = pm.beta_calculate(
                ion_pair=cap,
                ionic_strength=self.get_ionic_strength(x),
                T=self.T,
                parameters = parameters
            )
        return pair_parameters

    def get_b_phi(self, x):
        pair_parameters = {}
        cation_anion_pairs = self.component_groups['cation_anion_pairs']
        for cap in cation_anion_pairs:
            parameters = self.param_db[self.order_str(cap)]
            pair_parameters[cap] = pm.beta_phi_calculate(
                ion_pair=cap,
                ionic_strength=self.get_ionic_strength(x),
                T=self.T,
                parameters = parameters
            )

        return pair_parameters

    def get_b_prime(self, x):

        """
        Get Betas ready for calculating the "F" function.
        :return: betas for each type of salt
        """

        pair_parameters = {}
        cation_anion_pairs = self.component_groups['cation_anion_pairs']

        for cap in cation_anion_pairs:
            parameters = self.param_db[self.order_str(cap)]
            pair_parameters[cap] = pm.beta_prime_calculate(
                ion_pair=cap,
                ionic_strength=self.get_ionic_strength(x),
                T=self.T,
                parameters = parameters
            )
        return pair_parameters

    def get_c(self):
        pair_parameters = {}
        cation_anion_pairs = self.component_groups['cation_anion_pairs']
        for cap in cation_anion_pairs:
            parameters = self.param_db[self.order_str(cap)]
            pair_parameters[cap] = pm.c_calculate(
                ion_pair=cap,
                T=self.T,
                parameters = parameters
            )
        return pair_parameters

    def get_cc_phi(self, x):
        cation_pairs = self.component_groups['cation_pairs']
        phi_dict = {}
        for cap in cation_pairs:
            parameters = self.param_db[self.order_str(cap)]
            phi_dict[cap] = pm.cc_phi_calculate(
                ion_pair=cap,
                ionic_strength=self.get_ionic_strength(x),
                a_phi=self.get_a_phi(),
                T=self.T,
                parameters = parameters
            )
        return phi_dict

    def get_cc_phi_prime_phi(self, x):
        i = self.get_ionic_strength(x)
        dict = {}
        phis = self.get_cc_phi(x)
        for pair in phis.keys():
            dict[pair] = phis[pair]["phi"] + i * phis[pair]["phi_prime"]
        return dict

    def get_aa_phi(self, x):
        if self.component_groups['anion_pairs']:
            anion_pairs = self.component_groups['anion_pairs']
            phi_dict = {}
            for aap in anion_pairs:
                parameters = self.param_db[self.order_str(aap)]
                phi_dict[aap] = pm.aa_phi_calculate(
                    ion_pair=aap,
                    ionic_strength=self.get_ionic_strength(x),
                    a_phi=self.get_a_phi(),
                    T=self.T,
                    parameters=parameters
                )
            return phi_dict
        return {}

    def get_aa_phi_prime_phi(self, x):
        if self.component_groups['anion_pairs']:
            i = self.get_ionic_strength(x)
            dict = {}
            phis = self.get_aa_phi(x)
            for pair in phis.keys():
                dict[pair] = phis[pair]["phi"] + i * phis[pair]["phi_prime"]
            return dict
        else: 
            return {}

    def get_cca_psi(self):
        
        cation_pairs = self.component_groups['cation_pairs']
        anions = self.component_groups['anions']
        cca_pairs = {}
        for cation_pair in cation_pairs:
            cation1 = cation_pair[0]
            cation2 = cation_pair[1]
            for anion in anions:
                parameters = self.param_db[self.order_str((cation1, cation2, anion))]
                ternary_parameters = pm.ternary_parameter_cal(self.T, parameters)
                cca_pairs[cation1, cation2, anion] = ternary_parameters[0]
        return cca_pairs

    def get_aac_psi(self):
        anion_pairs = self.component_groups['anion_pairs']
        if anion_pairs:
            cations = self.component_groups['cations']
            aac_pairs = {}
            for anion_pair in anion_pairs:
                anion_list = list(anion_pair)
                a1 = anion_list[0]
                a2 = anion_list[1]
                for cation in cations:
                    parameters = self.param_db[self.order_str((a1, a2, cation))]
                    ternary_parameters = pm.ternary_parameter_cal(self.T, parameters)
                    aac_pairs[a1, a2, cation] =ternary_parameters[0]
            return aac_pairs
        else:
            return {}

    def get_lambdas(self):
        neutral_ion_pairs = self.component_groups['neutral_ion_pairs']
        lambdas = {}
        for pair in neutral_ion_pairs:
            parameters = self.param_db[self.order_str(pair)]
            lambdas[pair] = pm.lambda_cal(pair, self.T, parameters)
        return lambdas

    def get_zetas(self):
        nca_pairs = self.component_groups['neutral_cation_anion_pairs']
        zetas = {}
        for pair in nca_pairs:
            parameters = self.param_db[self.order_str(pair)]
            zetas[pair] = pm.zeta_cal(pair, self.T, parameters)
        return zetas

    def get_osmotic_coefficient(self, x):
        molalities = self.get_molalities(x)
        m_sum = sum(molalities.values())
        i = self.get_ionic_strength(x)
        a_phi = self.get_a_phi()
        z = self.get_z(x)
        b = self.get_b_phi(x)
        c = self.get_c()
        cations = self.component_groups['cations']
        cation_pairs = self.component_groups['cation_pairs']
        anions = self.component_groups['anions']
        anion_pairs = self.component_groups['anion_pairs']
        cc_phis = self.get_cc_phi_prime_phi(x)
        cca_psis = self.get_cca_psi()

        neutrals = self.component_groups['neutrals']
        lambdas = self.get_lambdas()
        zetas = self.get_zetas()
        nca_pairs = self.component_groups['neutral_cation_anion_pairs']

        item0 = (2 / m_sum)

        item1 = -(a_phi * i ** 1.5) / (1 + 1.2 * i ** 0.5)

        item2 = 0
        for cap in b.keys():
            m_1 = molalities[cap[0]]
            m_2 = molalities[cap[1]]

            b_phi = 0
            for key in b.keys():
                if set(cap) == set(key):
                    b_phi = b[key]

            c_value = 0
            for key in c.keys():
                if set(cap) == set(key):
                    c_value = c[key]

            item2 += m_1 * m_2 * (b_phi + z * c_value)

        item3 = 0
        for cp in cation_pairs:
            c1 = cp[0]
            c2 = cp[1]

            m_c1 = molalities[c1]
            m_c2 = molalities[c2]

            # find value of cc_phi
            cc_phi = 0
            for key in cc_phis.keys():
                if {cp[0], cp[1]} == set(key):
                    cc_phi = cc_phis[key]

            item3_subitem1 = 0
            for anion in anions:
                m_a = molalities[anion]

                # find value of cca_psi
                cca_psi = 0
                for key in cca_psis.keys():
                    if set(key) == {c1, c2, anion}:
                        cca_psi = cca_psis[key]

                item3_subitem1 += m_a * cca_psi

            item3 += m_c1 * m_c2 * (cc_phi + item3_subitem1)

        item4 = 0

        if len(anions) > 1:
            aa_phis = self.get_aa_phi_prime_phi(x)
            aac_psis = self.get_aac_psi()

            for ap in anion_pairs:
                m_a1 = molalities[ap[0]]
                m_a2 = molalities[ap[1]]

                # find value of aa_phi
                aa_phi = 0
                for key in aa_phis.keys():
                    if {ap[0], ap[1]} == set(key):
                        aa_phi = aa_phis[key]

                item4_subitem1 = 0
                for cation in cations:
                    m_c = molalities[cation]

                    aac_psi = 0

                    # find value of aac_psi
                    for key in cca_psis.keys():
                        if {ap[0], ap[1], cation} == set(key):
                            aac_psi = aac_psis[key]

                    item4_subitem1 += m_c * aac_psi
                item4 += m_a1 * m_a2 * (aa_phi + item4_subitem1)

        item5 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]
            for cation in cations:
                m_c = molalities[cation]

                lambda_nc = 0
                for key in lambdas.keys():
                    if {neutral, cation} == set(key):
                        lambda_nc = lambdas[key]
                item5 += m_n * m_c * lambda_nc

        item6 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]
            for anion in anions:
                m_a = molalities[anion]

                lambda_na = 0
                for key in lambdas.keys():
                    if {neutral, anion} == set(key):
                        lambda_na = lambdas[key]
                item6 += m_n * m_a * lambda_na

        item7 = 0
        for nca_pair in nca_pairs:
            m1 = molalities[nca_pair[0]]
            m2 = molalities[nca_pair[1]]
            m3 = molalities[nca_pair[2]]

            zeta_nca = 0
            for key in zetas.keys():
                if set(nca_pair) == set(key):
                    zeta_nca = zetas[key]

            item7 += m1 * m2 * m3 * zeta_nca

        osmotic_coefficient = item0 * (item1 + item2 + item3 + item4 + item5 + item6 + item7) + 1
        return osmotic_coefficient

    def get_water_activity(self, x):
        m_sum = sum(self.get_molalities(x).values())
        phi = self.get_osmotic_coefficient(x)
        ln_activity = -phi * m_sum / 55.50844
        return ln_activity

    def get_f_uppercase(self, x):
        a_phi = self.get_a_phi()
        i = self.get_ionic_strength(x)
        component_groups = self.component_groups
        cation_pairs = component_groups['cation_pairs']
        anion_pairs = component_groups['anion_pairs']
        molalities = self.get_molalities(x)
        ca_beta_primes = self.get_b_prime(x)
        cation_anion_pairs = component_groups['cation_anion_pairs']
        cc_phis = self.get_cc_phi(x)
        aa_phis = self.get_aa_phi(x)
        f_gamma = pm.get_f_gamma(a_phi, i)

        item0 = f_gamma

        item1 = 0
        for cap in cation_anion_pairs:
            m_p1 = molalities[cap[0]]
            m_p2 = molalities[cap[1]]
            ca_beta_prime = 0

            # find the value of beta_prime of cap
            for key in ca_beta_primes.keys():
                if {cap[0], cap[1]} == set(key):
                    ca_beta_prime = ca_beta_primes[key]

            item1 += m_p1 * m_p2 * ca_beta_prime

        item2 = 0
        for cp in cation_pairs:
            m_p1 = molalities[cp[0]]
            m_p2 = molalities[cp[1]]

            cc_phi_prime = 0

            # find the value of beta_prime of cap
            for key in cc_phis.keys():
                if {cp[0], cp[1]} == set(key):
                    cc_phi_prime = cc_phis[key]['phi_prime']

            item2 += m_p1 * m_p2 * cc_phi_prime

        item3 = 0
        if anion_pairs:
            for ap in anion_pairs:
                m_a1 = molalities[ap[0]]
                m_a2 = molalities[ap[1]]

                aa_phi_prime = 0
                for key in aa_phis.keys():
                    if {ap[0], ap[1]} == set(key):
                        aa_phi_prime = aa_phis[key]['phi_prime']
                item3 += m_a1 * m_a2 * aa_phi_prime

        f_uppercase = item0 + item1 + item2 + item3
        return f_uppercase

    def get_cation_activity_coefficients(self, target_cation, x):

        cations = self.component_groups['cations']
        anions = self.component_groups['anions']
        anion_pairs =self.component_groups['anion_pairs']
        cation_anion_pairs = self.component_groups['cation_anion_pairs']
        neutrals = self.component_groups['neutrals']
        molalities = self.get_molalities(x)
        target_cation = target_cation
        betas = self.get_b(x)
        cs = self.get_c()
        z = self.get_z(x)
        lambdas = self.get_lambdas()
        zetas = self.get_zetas()
        cc_phis = self.get_cc_phi(x)
        cation_psis = self.get_cca_psi()
        anion_psis = self.get_aac_psi()

        f_uppercase = self.get_f_uppercase(x)
        charge_number = get_charge_number(target_cation)

        item0 = charge_number ** 2 * f_uppercase

        item1 = 0
        for anion in anions:
            m_a = molalities[anion]
            beta_ca = 0
            for key in betas.keys():
                if {target_cation, anion} == set(key):
                    beta_ca = betas[key]
            c_ca = 0
            for key in cs.keys():
                if {target_cation, anion} == set(key):
                    c_ca = cs[key]
            item1 += m_a * (2 * beta_ca + z * c_ca)

        item2 = 0
        for cation in cations:
            if cation != target_cation:
                m_c = molalities[cation]
                phi_mc = 0
                for key in cc_phis.keys():
                    if set(key) == {cation, target_cation}:
                        phi_mc = cc_phis[key]['phi']

                item2_subitem1 = 0
                for anion in anions:
                    m_a = molalities[anion]
                    for tkey in cation_psis.keys():
                        if {target_cation, cation, anion} == set(tkey):
                            item2_subitem1 += m_a * cation_psis[tkey]
                item2 += m_c * (2 * phi_mc + item2_subitem1)

        item3 = 0
        if anion_pairs:
            for ap in anion_pairs:
                m_a1 = molalities[ap[0]]
                m_a2 = molalities[ap[1]]
                psi_maa = 0

                # find match for psi
                for key in anion_psis.keys():
                    if {ap[0], ap[1], target_cation} == set(key):
                        psi_maa = anion_psis[key]

                item3 += m_a1 * m_a2 * psi_maa

        item4 = 0
        for cap in cation_anion_pairs:
            ion1 = cap[0]
            ion2 = cap[1]
            c_ca = 0

            # find match for c_ca
            for key in cs.keys():
                if set(key) == {ion1, ion2}:
                    c_ca = cs[key]

            item4 += molalities[ion1] * molalities[ion2] * c_ca

        item5 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]

            lambda_nm = 0
            for key in lambdas.keys():
                if {neutral, target_cation} == set(key):
                    lambda_nm = lambdas[key]
            item5 += m_n * lambda_nm

        item6 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]
            for anion in anions:
                m_a = molalities[anion]

                zeta_nam = 0
                for key in zetas.keys():
                    if {target_cation, anion, neutral} == set(key):
                        zeta_nam = zetas[key]
                item6 += m_n * m_a * zeta_nam

        ln_coefficient = item0 + item1 + item2 + item3 + abs(charge_number) * item4 + 2 * item5 + item6
        return ln_coefficient

    def get_anion_activity_coefficients(self, target_anion, x):
        groups = self.component_groups
        cations = groups['cations']
        anions = groups['anions']
        cation_pairs = groups['cation_pairs']
        cation_anion_pairs = groups['cation_anion_pairs']
        neutrals = groups['neutrals']
        lambdas = self.get_lambdas()
        zetas = self.get_zetas()
        molalities = self.get_molalities(x)
        target_anion = target_anion
        betas = self.get_b(x)
        cs = self.get_c()
        z = self.get_z(x)

        aa_phis = self.get_aa_phi(x)

        aac_psis = self.get_aac_psi()
        cca_psis = self.get_cca_psi()

        f_uppercase = self.get_f_uppercase(x)
        charge_number = get_charge_number(target_anion)

        item0 = charge_number ** 2 * f_uppercase

        item1 = 0
        for cation in cations:
            m_a = molalities[cation]
            beta_ca = 0
            for key in betas.keys():
                if {target_anion, cation} == set(key):
                    beta_ca = betas[key]
            c_ca = 0
            for key in cs.keys():
                if {target_anion, cation} == set(key):
                    c_ca = cs[key]
            item1 += m_a * (2 * beta_ca + z * c_ca)

        item2 = 0
        for anion in anions:
            if anion != target_anion:
                m_a = molalities[anion]
                phi_xa = 0
                for key in aa_phis.keys():
                    if set(key) == {anion, target_anion}:
                        phi_xa = aa_phis[key]['phi']

                item2_subitem1 = 0
                for cation in cations:
                    m_c = molalities[cation]
                    for tkey in aac_psis.keys():
                        if {target_anion, cation, anion} == set(tkey):
                            item2_subitem1 += m_c * aac_psis[tkey]
                item2 += m_a * (2 * phi_xa + item2_subitem1)

        item3 = 0
        for cp in cation_pairs:
            m_c1 = molalities[cp[0]]
            m_c2 = molalities[cp[1]]
            psi_ccx = 0

            # find match for psi
            for key in cca_psis.keys():
                if {cp[0], cp[1], target_anion} == set(key):
                    psi_ccx = cca_psis[key]
            item3 += m_c1 * m_c2 * psi_ccx

        item4 = 0
        for cap in cation_anion_pairs:
            ion1 = cap[0]
            ion2 = cap[1]
            m1 = molalities[ion1]
            m2 = molalities[ion2]

            c_ca = 0
            for key in cs.keys():
                if set(key) == {ion1, ion2}:
                    c_ca = cs[key]

            item4 += m1 * m2 * c_ca

        item5 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]
            lambda_nx = 0

            for key in lambdas.keys():
                if {neutral, target_anion} == set(key):
                    lambda_nx = lambdas[key]

            item5 += m_n * lambda_nx

        item6 = 0
        for neutral in neutrals:
            m_n = molalities[neutral]
            for cation in cations:
                m_c = molalities[cation]

                zeta_nam = 0

                for key in zetas.keys():
                    if {target_anion, cation, neutral} == set(key):
                        zeta_nam = zetas[key]
                item6 += m_n * m_c * zeta_nam

        ln_coefficient = item0 + item1 + item2 + item3 + abs(charge_number) * item4 + 2 * item5 + item6

        return ln_coefficient

    def get_neutral_activity_coefficients(self, target_neutral, x):
        groups = self.component_groups
        cations = groups['cations']
        anions = groups['anions']
        cation_anion_pairs = groups['cation_anion_pairs']

        molalities = self.get_molalities(x)

        lambdas = self.get_lambdas()
        zetas = self.get_zetas()

        item0 = 0
        for cation in cations:
            m_c = molalities[cation]
            lambda_nc = 0
            for key in lambdas.keys():
                if {cation, target_neutral} == set(key):
                    lambda_nc = lambdas[key]
            item0 += m_c * lambda_nc

        item1 = 0
        for anion in anions:
            m_a = molalities[anion]
            lambda_na = 0
            for key in lambdas:
                if {anion, target_neutral} == set(key):
                    lambda_na = lambdas[key]
            item1 += m_a * lambda_na

        item2 = 0
        for pair in cation_anion_pairs:
            m1 = molalities[pair[0]]
            m2 = molalities[pair[1]]

            zeta_nca = 0
            for key in zetas.keys():
                if {target_neutral, pair[0], pair[1]} == set(key):
                    zeta_nca = zetas[key]

            item2 += m1 * m2 * zeta_nca

        ln_coefficient = 2 * item0 + 2 * item1 + item2

        return ln_coefficient

    def total_gibbs_energy(self, x):
        r = 8.314
        molalities = self.get_molalities(x)
        total_g = 0
        for species in molalities.keys():
            species_type = type_of_species(species)
            m_i = molalities[species]
            if m_i != 0:
                ln_gamma = 0
                if species_type == 0:
                    ln_gamma = self.get_neutral_activity_coefficients(species, x)
                elif species_type == 1:
                    ln_gamma = self.get_cation_activity_coefficients(species, x)
                elif species_type == -1:
                    ln_gamma = self.get_anion_activity_coefficients(species, x)
                total_g += m_i * (0 + r * self.T * (np.log(m_i) + ln_gamma))

        # gibbs energy of water
        ln_a_w = self.get_water_activity(x)

        water_data = db.get_reaction("H2(g) + 0.5O2(g) = H2O(l)", self.t)
        water_parameters = water_data['analytic']

        std_cp_water = pm.get_chemical_potential(water_parameters, self.T, water_data['eq_num'])

        cp_water = (1000 / 18.015) * (std_cp_water + r * self.T * ln_a_w)

        total_g = total_g + cp_water
        return total_g

    def solubility_equilibrium(self, x):
        # print(x)
        lna_pitzer = self.get_water_activity(x)
        if self.equilibrium == 'H2O(s) = H2O(l)':
            lnk_ice = clegg_and_brimblecombe(self.T)
            # lnk_ice = monnin(self.t)
            # lnk_ice = spencer(self.t)
            f = lnk_ice - lna_pitzer
        else:
            solid_data = self.solid_disso['stoich']
            parameters = self.solid_disso["analytic"]

            molalities = self.get_molalities(x)
            lnk_potential = pm.get_chemical_potential(parameters, self.T, self.solid_disso["eq_num"])

            lnk_activity = 0
            for species in solid_data.keys():
                # get the stochiometric number of this species first
                sto = solid_data[species]['value']

                if solid_data[species]['type'] == 'cation':
                    m_c = molalities[species]
                    ln_gamma_c = self.get_cation_activity_coefficients(species, x)
                    lnk_activity += sto * (np.log(m_c) + ln_gamma_c)
                elif solid_data[species]['type'] == 'anion':
                    m_a = molalities[species]
                    ln_gamma_a = self.get_anion_activity_coefficients(species, x)
                    lnk_activity += sto * (np.log(m_a) + ln_gamma_a)
                else:
                    lnk_activity += sto * lna_pitzer

            # if a neutral species should be considered.
            if self.neutral:
                ln_gamma_n = self.get_neutral_activity_coefficients(target_species, x)

                lnk_activity = (np.log(molalities[target_species]) + ln_gamma_n) - lnk_activity
            f = lnk_potential - lnk_activity
        return f

    def solve(self):
        # Start the timer
        start_time = time.time()

        bounds = ((0, None), (0, None))
        # bounds = ((0, None))
        cons1 = {'type': 'eq', 'fun': self.charge_balance}
        cons2 = {'type': 'eq', 'fun': self.solubility_equilibrium}
        res = minimize(
            fun=self.total_gibbs_energy,
            x0=self.x0,
            method='SLSQP',
            bounds=bounds,
            constraints=(cons1, cons2)
            # constraints=(cons1)
        )

        # end of timer
        elapsed_time = time.time() - start_time
        # print('times used:', elapsed_time)

        if res.success:
            return res.x
        else:
            return [None, None]
