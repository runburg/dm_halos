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
    with np.load(file+"_nounits.txt", 'rb') as npzfile:
        r = npzfile['r']
        v = npzfile['v']
        fe = npzfile['fe']

    # create a list of unique r values and how often they occur
    r_unique = np.unique(r)

    def rho(x):
        return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)

    i = 0
    # initial arrays for grabbing parts of the data
    v_temp = []
    func1 = []
    func2 = []
    rf = []
    g_dwave = []
    # loop through all of the unique values of r
    for rad in r_unique:
        # for each set of (v,fe) that correspond to the given r, create [x]
        # and [y]
        # for num. int.
        while rad == r[i]:
            # [x] for integration
            v_temp.append(v[i])
            # [y] for integration
            func1.append(fe[i]*v[i]**6)
            func2.append(fe[i]*v[i]**4)
            i += 1
            # abort final loop to avoid out of bounds error
            if i >= len(r):
                break
        # stores the value of the velocity integration
        # this will change when not doing s-wave
        g_dwave.append(8*mp.pi*rho(rad)*integrate.simps(func1, v_temp)
                       + 160*mp.pi**2 / 3*integrate.simps(func2, v_temp)**2)
        v_temp.clear()
        func1.clear()
        func2.clear()
        rf.append(rad)

    with open("g_d.txt", 'wb') as outfile:
        np.savez(outfile, gd=np.array(g_dwave), r=np.array(rf))
