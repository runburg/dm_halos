#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 19:44:24 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp

# set working precision of decimal points
# mp.dps = 25

'''s-wave annihilation'''


def gs_wave(file):
    """Perform integration of velocity integrals."""
    with np.load(file+"_nounits.txt", 'rb', allow_pickle=True) as npzfile:
        r = npzfile['r']
        v = npzfile['v'].astype(np.float)
        fe = npzfile['fe'].astype(np.float)

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r, return_index=True)
    # initial arrays for grabbing parts of the data
    v_temp = []
    integrand = []
    g_swave = []
    # loop through all of the unique values of r
    for i, (rad, j) in enumerate(np.array(r_unique).T, start=1):
        # for each set of (v,fe) that correspond to the given r, create [x]
        # and [y]
        # for num. int.
        if i == len(r_unique[0]):
            i=None
        else:
            i=r_unique[1][i]

        v_temp = v[j:i]
        integrand = v_temp**2 * fe[j:i]
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        try:
            func = interpolate.interp1d(v_temp, integrand, kind='cubic', fill_value='extrapolate')
            g_swave.append(16*mp.pi**2*integrate.quad(func, 0, v_temp[-1], points=v_temp, limit=10000)[0]**2)
            # g_swave.append(16*mp.pi**2*(integrate.simps(integrand, v_temp))**2)
        except ValueError:
            print(i, j, integrand, v_temp)
            g_swave.append(16*mp.pi**2*(integrate.trapz(integrand, v_temp))**2)


    with open("./dimensionless_j_factors/df_nfw/df_nfw_g_s.txt", 'wb') as outfile:
        np.savez(outfile, g_s=np.array(g_swave), r=np.array(r_unique[0]))


if __name__ == '__main__':
    gs_wave('./dimensionless_j_factors/df_nfw/df_nfw')
