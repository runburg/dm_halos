#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 20:56:58 2019

@author: runburg
"""

from scipy import integrate
import numpy as np
import mpmath as mp

# set working precision of decimal points
# mp.dps = 25

'''p-wave annihilation'''
infile = open("fe_GC_NFW_nounits.txt", 'rb')
npzfile = np.load(infile)
r = npzfile['r']
v = npzfile['v']
fe = npzfile['fe']
infile.close()

# create a list of unique r values and how often they occur
r_unique = np.unique(r)


def rho(x):
    return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)


i = 0
# initial arrays for grabbing parts of the data
v_temp = []
func = []
rf = []
g_pwave = []
# loop through all of the unique values of r
for rad in r_unique:
    # for each set of (v,fe) that correspond to the given r, create [x] and [y]
    # for num. int.
    while rad == r[i]:
        # [x] for integration
        v_temp.append(v[i])
        # [y] for integration
        func.append(fe[i]*v[i]**4)
        i += 1
        # abort final loop to avoid out of bounds error
        if i >= len(r):
            break
    # stores the value of the velocity integration
    # this will change when not doing s-wave
    g_pwave.append(8*np.pi*rho(rad)*integrate.simps(func, v_temp))
    v_temp.clear()
    func.clear()
    rf.append(rad)

outfile = open("g_p.txt", 'wb')
np.savez(outfile, gp=np.array(g_pwave), r=np.array(rf))
outfile.close()
