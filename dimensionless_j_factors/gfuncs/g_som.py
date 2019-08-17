#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Som. enh. velocity integration for dm project.

The annihilation scheme of sommerfeld enhancement.
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp

# set working precision of decimal points
mp.dps = 50

'''sommerfeld-enhanced annihilation'''


def gsom_wave(file):
    """Integrate over velocity distribution for som. enh. dm."""
    with np.load(file+"_nounits.txt", 'rb', allow_pickle=True) as npzfile:
        r = npzfile['r']
        v = npzfile['v']
        fe = npzfile['fe']

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r, return_index=True)

    print(len(r_unique[0]))
    # initial arrays for grabbing parts of the data
    v_temp = []
    fe_temp = []
    integrand = []
    rf = []
    g_som = []
    # loop through all of the unique values of r
    for i, (rad, j) in enumerate(np.array(r_unique).T, start=1):
        print(i)
        # for each set of (v,fe) that correspond to the given r, create [x]
        # and [y]
        # for num. int.
        if i == len(r_unique[0]):
            i=None
        else:
            i=r_unique[1][i]

        v_temp = v[j:i].astype(np.float)
        fe_temp = fe[j:i].astype(np.float)
        # integrand = [vel**2]

        if len(v_temp)>2:
            func1 = interpolate.interp1d(v_temp, fe_temp*v_temp**2, kind='quadratic')
            func2 = interpolate.interp1d(v_temp, fe_temp*v_temp, kind='quadratic')
            g_som.append(np.sqrt(4*np.pi)*32*np.pi**2*integrate.dblquad(lambda v2, v1: func1(v1)*func2(v2), v_temp[0], v_temp[-1], lambda v1: v1, lambda v1: v_temp[-1], epsabs=1e-04, epsrel=1e-04)[0])
        else:
            g_som.append(0)

    # save in g_som.txt
    with open("./dimensionless_j_factors/df_nfw/df_nfw_g_som.txt", 'wb') as outfile:
        np.savez(outfile, g_som=np.array(g_som), r=np.array(r_unique[0]))


if __name__ == '__main__':
    gsom_wave('./dimensionless_j_factors/df_nfw/df_nfw')
