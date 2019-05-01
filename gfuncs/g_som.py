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
    with np.load(file+"_nounits.txt", 'rb') as npzfile:
        r = npzfile['r']
        v = npzfile['v']
        fe = npzfile['fe']

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r)

    i = 0
    # initial arrays for grabbing parts of the data
    v_temp = []
    fe_temp = []
    func = []
    rf = []
    g_som = []
    # loop through all of the unique values of r
    for rad in r_unique:
        # for each set of (v,fe) that correspond to the given r, create [x]
        # and [y]
        # for num. int.
        while rad == r[i]:
            # [x] for integration
            v_temp.append(v[i])
            # [y] for integration
            fe_temp.append(fe[i])
            i += 1
            # abort final loop to avoid out of bounds error
            if i >= len(r):
                break
        for vv, ffee in zip(v_temp, fe_temp):
            x = np.array(v_temp).astype(float)
            y = np.array([vvv**2*fffeee for (vvv, fffeee)
                          in zip(v_temp, fe_temp)])
            if len(x) > 1 and len(y) > 1:
                v2int = interpolate.interp1d(
                    x, y, kind='linear', fill_value='extrapolate')
            func.append(ffee*vv*integrate.quad(v2int, vv, v_temp[-1])[0])
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        g_som.append(32*mp.pi**2*integrate.simps(func, v_temp))
        fe_temp.clear()
        v_temp.clear()
        func.clear()
        rf.append(rad)

    # save in g_som.txt
    with open("g_som.txt", 'wb') as outfile:
        np.savez(outfile, g_som=np.array(g_som), r=np.array(rf))
