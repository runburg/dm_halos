#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Jan 10 20:56:58 2019.

@author: runburg
"""

from scipy import integrate
import numpy as np
import mpmath as mp

# set working precision of decimal points
# mp.dps = 25

'''p-wave annihilation'''


def gp_wave(file):
    with np.load(file+"_nounits.txt", 'rb', allow_pickle=True) as npzfile:
        r = npzfile['r']
        v = npzfile['v']
        fe = npzfile['fe']

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r, return_index=True)

    def rho(x):
        return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)

    # initial arrays for grabbing parts of the data
    v_temp = []
    integrand = []
    g_pwave = []
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
        integrand = v_temp**4 * fe[j:i]
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        # extra 4 pi from redfinition in paper
        g_pwave.append(1/(4*np.pi)*8*mp.pi*rho(rad)*(integrate.simps(integrand, v_temp)))

    with open("./dimensionless_j_factors/df_nfw/df_nfw_g_p.txt", 'wb') as outfile:
        np.savez(outfile, g_p=np.array(g_pwave), r=np.array(r_unique[0]))


if __name__ == "__main__":
    gp_wave('./dimensionless_j_factors/df_nfw/df_nfw')
