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
    """Perform transformation to function of velocity and radius."""
    with np.load(file+"_nounits.txt", 'rb') as npzfile:
        r = npzfile['r']
        v = npzfile['v']
        fe = npzfile['fe']

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r)

    i = 0
    # initial arrays for grabbing parts of the data
    v_temp = []
    func = []
    rf = []
    g_swave = []
    # loop through all of the unique values of r
    for rad in r_unique:
        # for each set of (v,fe) that correspond to the given r, create [x]
        # and [y]
        # for num. int.
        while rad == r[i]:
            # [x] for integration
            v_temp.append(v[i])
            # [y] for integration
            func.append(fe[i]*v[i]*v[i])
            i += 1
            # abort final loop to avoid out of bounds error
            if i >= len(r):
                break
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        vels = np.array(v_temp).astype(np.float)
        if len(v_temp) > 3:
            g_sfunc = interpolate.interp1d(vels, func, kind='cubic', fill_value='extrapolate')
            g_swave.append(16*np.pi**2*integrate.quad(g_sfunc, 0, vels[-1], points=vels, limit=100*len(v_temp))[0]**2)
        # print('{}\t{}\t{}'.format(len(v_temp), v_temp[0], v_temp[-1]))
        else:
            g_swave.append(16*mp.pi**2*(integrate.simps(func, v_temp))**2)
        v_temp.clear()
        func.clear()
        rf.append(rad)

    with open("/Users/runburg/github/dm_halos/df_nfw/df_nfw_g_s.txt", 'wb') as outfile:
        np.savez(outfile, g_s=np.array(g_swave), r=np.array(rf))


if __name__ == '__main__':
    gs_wave('/Users/runburg/github/dm_halos/df_nfw/df_nfw')
