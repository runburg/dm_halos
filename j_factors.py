#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate dimensionful J-factors.

Author: Jack Runburg
Date: 20-05-2019 10:11


"""

from scipy import interpolate
import numpy as np


def weighted_average(listofvalues):
    """Compute weighted average of list (a,b) with 'a' as the weights and 'b' as the values."""
    sum = 0
    weights = 0
    for [w, v] in listofvalues:
        sum += w*v
        weights += w
    return sum/weights


def weighted_sd(listofvalues):
    """Compute weighted standard deviation of list (a,b) with 'a' as the weights and 'b' as the values."""
    ave = weighted_average(listofvalues)
    sum = 0
    weights = 0
    for [w, v] in listofvalues:
        sum += w*(v-ave)**2
        weights += w
    return np.sqrt(sum*len(listofvalues)/((len(listofvalues)-1)*weights))


def sigma_upper(listofvalues, ave):
    """Compute the standard deviation of values above the mean."""
    above = []
    for [w, v] in listofvalues:
        if v > ave:
            above.append(v)
    above.sort()
    i = int(np.floor(len(listofvalues)*.34134))
    if i >= len(above):
        # return weighted_sd(listofvalues)
        return above[-1]
    else:
        return above[i]


def sigma_lower(listofvalues, ave):
    """Compute the standard deviation of values below the mean."""
    lower = []
    for [w, v] in listofvalues:
        if v < ave:
            lower.append(v)
    lower.sort()
    lower.reverse()
    i = int(np.floor(len(listofvalues)*.34134))
    if i >= len(lower):
        return weighted_sd(listofvalues)
    else:
        return lower[i]


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
dwarflist = ['draco1', 'segue1']
for dwarf in dwarflist:
    data = np.load('./j_factors/dwarfs/'+dwarf+'_chain.npy')

    ave_jvals, up_jvals, low_jvals = [], [], []
    thetas = np.logspace(-5, np.log10(np.pi/180*3), num=50)
    for wave in ['s', 'p', 'd', 'som']:
        ave_jvals, low_jvals, up_jvals = [], [], []
        for t in thetas:
            temp_j = []
            weightsum = 0
            for x in data:
                if wave == 's':
                    temp_j.append([x[0], jsfunc(t*x[1]/(10**x[2]))*2*10**x[2]*(10**x[3])**2])
                if wave == 'p':
                    temp_j.append([x[0], jpfunc(t*x[1]/(10**x[2]))*2*10**x[2]*(10**x[3])**2*4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792.458)**2])
                if wave == 'd':
                    temp_j.append([x[0], jdfunc(t*x[1]/(10**x[2]))*2*10**x[2]*(10**x[3])**2*(4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792.458)**2)**2])
                if wave == 'som':
                    temp_j.append([x[0], jsomfunc(t*x[1]/(10**x[2]))*2*10**x[2]*(10**x[3])**2/np.sqrt(4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792.458)**2)])
                weightsum += x[0]
            temp_j = np.array(temp_j)
            temp_j = temp_j[temp_j[:, 1].argsort()]
            temp_j[:, 0] = np.cumsum(temp_j, axis=0)[:, 0]
            temp_j = np.divide(temp_j, np.array([weightsum, 1]))
            ave_jvals.append(temp_j[np.searchsorted(temp_j[:, 0], 0.5, side='left')][1])
            low_jvals.append(temp_j[np.searchsorted(temp_j[:, 0], .5-.34, side='left')][1])
            up_jvals.append(temp_j[np.searchsorted(temp_j[:, 0], .5+.34, side='right')][1])
        #     temp_ave = weighted_average(temp_j)
        #     ave_jvals.append(temp_ave)
        #     low_jvals.append(sigma_lower(temp_j, temp_ave))
        #     up_jvals.append(sigma_upper(temp_j, temp_ave))
        #     # sigma_upper(temp_j, temp_ave)
        #
        with open("./j_factors/j_ave/j_ave_thetas_"+dwarf+"_"+wave+".txt", 'wb') as outfile:
            np.savez(outfile, ave_jvals=np.array(ave_jvals), up_jvals=np.array(up_jvals), low_jvals=np.array(
                low_jvals), thetas=np.array(thetas))
