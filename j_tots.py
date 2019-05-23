#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Check J_tot values.

Author: Jack Runburg
Date: 22-05-2019 16:14


"""
from scipy import interpolate, integrate
import numpy as np


def weighted_average(listofvalues):
    """Compute weighted average of list (a,b) with 'a' as the weights and 'b' as the values."""
    sum = 0
    weights = 0
    for [w, v] in listofvalues:
        sum += w*v
        weights += w
    return sum/weights


# loading j function values, function of tilde theta!
with np.load('./j_factors/j_thetas.txt', 'rb') as infile:
    jst = infile['jst']
    js = infile['js']
    jp = infile['jp']
    jd = infile['jd']
    jsom = infile['jsom']
    tildetheta = infile['theta']

jsfunc = interpolate.interp1d(tildetheta, js, kind='cubic', fill_value='extrapolate')
jpfunc = interpolate.interp1d(tildetheta, jp, kind='cubic', fill_value='extrapolate')
jdfunc = interpolate.interp1d(tildetheta, jd, kind='cubic', fill_value='extrapolate')
jsomfunc = interpolate.interp1d(tildetheta, jsom, kind='cubic', fill_value='extrapolate')

# dwarf name

dwarflist = ['draco', 'segue']
for dwarf in dwarflist:
    data = np.load('./j_factors/dwarfs/'+dwarf+'1_chain.npy')
    for wave in ['s', 'p', 'd', 'som']:
        ave_jvals = []
        for x in data:
            n = 0
            if wave == 's':
                n = 0
                jtot = integrate.quad(lambda y: y * jsfunc(0.5/180*np.pi*x[1]/10**x[2]), 0, 50)[0]
            if wave == 'p':
                n = 2
                jtot = integrate.quad(lambda y: y * jpfunc(0.5/180*np.pi*x[1]/10**x[2]), 0, 50)[0]
            if wave == 'd':
                n = 4
                jtot = integrate.quad(lambda y: y * jdfunc(0.5/180*np.pi*x[1]/10**x[2]), 0, 50)[0]
            if wave == 'som':
                n = -1
                jtot = integrate.quad(lambda y: y * jsomfunc(0.5/180*np.pi*x[1]/10**x[2]), 0, 50)[0]

            ave_jvals.append([x[0], 4*np.pi*(10**x[2])**3*(10**x[3])**2/x[1]**2*(4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792458)**2)**(n/2.)*jtot])
        print(dwarf+'\t'+wave+'\t'+str(np.log10(weighted_average(ave_jvals))))
