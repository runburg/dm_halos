#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 13:01:55 2019

@author: runburg
"""
from scipy import integrate
import numpy as np
import mpmath as mp

# set working precision of decimal points
mp.dps = 25

'''d-wave annihilation'''


def gd_wave(file):
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
    integrand1 = []
    integrand2 = []
    g_dwave = []
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
        integrand1 = 2*integrate.simps(4*np.pi*v_temp**6 * fe[j:i]/rho(rad), v_temp)
        integrand2 = 5/3*integrate.simps(4*np.pi*v_temp**4 * fe[j:i]/rho(rad), v_temp)**2
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        g_dwave.append(1/(4*np.pi)**2*rho(rad)**2*(integrand1+integrand2))

    with open("./dimensionless_j_factors/df_nfw/df_nfw_g_d.txt", 'wb') as outfile:
        np.savez(outfile, g_d=np.array(g_dwave), r=np.array(r_unique[0]))


if __name__ == '__main__':
    gd_wave('./dimensionless_j_factors/df_nfw/df_nfw')
