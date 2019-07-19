#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare to 1702.00408.

Author: Jack Runburg
Date: 10-07-2019 12:58


"""
import numpy as np


def convert_to_gev_cm5(num):
    """Perform the unit change."""
    c = 299792.458  # km/s
    kpc_to_cm = float(3.1E+19) * 100
    msun_to_gev = float(2E+30) * (c * 1000)**2 / float(1.602E-19) / 10**9
    # changeunits_to_gevcm = 37.96**2 * 10**(-18) * 1000 * 3.086 * 10**18
    return num * msun_to_gev**2 / kpc_to_cm**5


def kpc_to_cm(num):
    """Convert from kpc to cm."""
    return num * float(3.1E+19) * 100


def rmax_to_r_s(num):
    """Convert from rmax to rs."""
    r_s = num/2.16 * float(3.1E+19) * 100  # cm
    return r_s


def vmax_to_rho_s(num, rs, c=299792.458):
    """Convert from vmax to rhos."""
    g_n = float(6.7E-39)  # GeV^-2
    gev_to_cm = (1/0.197)*10000000000000  # cm^-1
    rho_s = (num/.465/c)**2/(4*np.pi*rs**2*g_n)*gev_to_cm
    return rho_s


def calculate_j_factor(list, wave):
    """Try to calculate j-factors with good units."""
    if wave == 's':
        n = 0
        cj = 1/3
    if wave == 'som':
        n = -1
        conversion_factor = 2 * np.pi * 0.01
        cj = 0.99 * conversion_factor
    c = 299792.458  # km/s
    # g_n = float(4.303E-6)  # kpc * (km/s)^2 / M_solar
    g_n = float(6.7E-39)
    gev_to_cm = (1/0.197)*10000000000000
    d = kpc_to_cm(list[1].astype(np.float))  # cm
    r_s = rmax_to_r_s(list[2].astype(np.float))  # cm
    rho_s = vmax_to_rho_s(list[3].astype(np.float), r_s)  # GeV/cm^3

    j = (4*np.pi*rho_s**2*r_s**3/d**2) * (4*np.pi*g_n*rho_s*r_s**2/gev_to_cm)**(n/2.) * cj  # (M_solar)^2/kpc^5
    return j


if __name__ == '__main__':
    dwarfs = np.array([['coma', 44, .38, 9.8],
                       ['ursa minor', 76, 1.32, 24.1],
                       ['draco', 76, 0.86, 17.7],
                       ['segue 1', 23, 0.76, 16.2],
                       ['reticulum 2', 32, 0.28, 7.6]
                       ])

    for dwarf in dwarfs:
        print('{}: {}'.format(dwarf[0], calculate_j_factor(dwarf, 's')))
