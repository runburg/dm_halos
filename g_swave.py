#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 19:44:24 2019

@author: runburg
"""

from scipy import integrate
import numpy as np
import mpmath as mp

# set working precision of decimal points
# mp.dps = 25

'''s-wave annihilation'''


def gs_wave(file):
    infile = open(file+"_nounits.txt", 'rb')
    npzfile = np.load(infile)
    r = npzfile['r']
    v = npzfile['v']
    fe = npzfile['fe']
    infile.close()

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
        g_swave.append(16*mp.pi**2*(integrate.simps(func, v_temp))**2)
        v_temp.clear()
        func.clear()
        rf.append(rad)

    with open("g_s.txt", 'wb') as outfile:
        np.savez(outfile, gs=np.array(g_swave), r=np.array(rf))
