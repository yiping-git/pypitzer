# -*- coding: utf-8 -*-
# Author: Yiping Liu
# Description:
# Version: 1.0
# Last Modified: May 7, 2023


def find_pair(x, pairs):
    """
    :param x: a tuple
    :param pairs: a list of tuples
    :return: a dict contains a boolean value and the match (if there is).
    """
    pair_length = len(pairs)
    i = 0
    has = False
    target = ''
    while i < pair_length:
        pair = pairs[i]
        if set(x) == set(pair):
            has = True
            target = pair
        i += 1
    return {
        "has": has,
        "target": target
    }


def get_charge_number(ion):
    """
    get the charge number of an ion (str)
    :param ion: ion name, a string with "+" or "-" sign followed by number of charge
    :return: charge number
    """
    if "+" in ion:
        lis = ion.split("+")
        if lis[1]:
            result = lis[1]
        else:
            result = 1
    elif "-" in ion:
        lis = ion.split("-")
        if lis[1]:
            lis[1] = '-' + lis[1]
            result = lis[1]
        else:
            result = -1
    else:
        result = 0
    return int(result)


def type_of_species(ion):
    """
    determine the type of a ion (string).
    :param ion:
    :return: 0-neutral, 1-cation, -1-anion
    """
    type_value = 0

    charge = get_charge_number(ion)
    if charge == 0:
        type_value = 0
    if charge > 0:
        type_value = 1
    if charge < 0:
        type_value = -1
    return type_value
