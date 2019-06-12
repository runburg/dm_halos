#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check J_tot values.

Author: Jack Runburg
Date: 22-05-2019 16:14

Compute dimensionful total j_factors inside solid angle given by ub.
Returns a weighted median and the asymmetric sigma values over given halo parameters
for each dwarf in /j_factors/dwarfs.
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


def OLD_integrated_j_factor(list, upper_bound, wave):
    """Return the value of the integrated j-factor for wave annihilation to angle upper_bound in degrees."""
    temp_j = []
    changeunits_to_gevcm = 37.96**2 * 10**(-18) * 1000 * 3.086 * 10**18
    for x in list:
        n = 0
        ub = upper_bound/180*np.pi*x[1]/10**x[2]
        if wave == 's':
            n = 0
            jtot = integrate.quad(lambda y: y * jsfunc(y), 0, ub, full_output=1)[0]
        if wave == 'p':
            n = 2
            jtot = integrate.quad(lambda y: y * jpfunc(y), 0, ub, full_output=1)[0]
        if wave == 'd':
            n = 4
            jtot = integrate.quad(lambda y: y * jdfunc(y), 0, ub, full_output=1)[0]
        if wave == 'som':
            n = -1
            jtot = integrate.quad(lambda y: y * jsomfunc(y), 0, ub, full_output=1)[0]

        temp_j.append([x[0], 4*np.pi*(10**x[2])**3*(10**x[3])**2/(x[1]**2)*(4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792.458)**2)**(n/2.)*jtot * changeunits_to_gevcm])

    return weighted_average(temp_j)


def weighted_median(values):
    """Compute the weighted median and return the +/- sigma values."""
    low_jvals, ave_jvals, up_jvals = [], [], []
    temp_j = np.asarray(values)
    weightsum = temp_j[:, 0].sum()
    temp_j = temp_j[temp_j[:, 1].argsort()]
    temp_j[:, 0] = np.cumsum(temp_j, axis=0)[:, 0]
    temp_j = np.divide(temp_j, np.array([weightsum, 1]))
    low_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5-.34, side='left')][1]
    ave_jvals = temp_j[np.searchsorted(temp_j[:, 0], 0.5, side='left')][1]
    up_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5+.34, side='right')][1]

    return np.array([low_jvals, ave_jvals, up_jvals])


def integrated_j_factor(list, upper_bound, wave):
    """Return the value of the integrated j-factor for wave annihilation to angle upper_bound in degrees."""
    temp_j = []
    changeunits_to_gevcm = 37.96**2 * 10**(-18) * 1000 * 3.086 * 10**18
    for x in list:
        n = 0
        ub = upper_bound/180*np.pi*x[1]/10**x[2]
        if wave == 's':
            n = 0
            jtot = jsfunc(ub)
        if wave == 'p':
            n = 2
            jtot = jpfunc(ub)
        if wave == 'd':
            n = 4
            jtot = jdfunc(ub)
        if wave == 'som':
            n = -1
            jtot = jsomfunc(ub)

        temp_j.append((x[0], 4*np.pi*(10**x[2])**3*(10**x[3])**2/(x[1]**2)*(4 * np.pi * float(4.325E-6) * 10**x[3] * (10**x[2])**2 / (299792.458)**2)**(n/2.)*jtot * changeunits_to_gevcm))

    return weighted_median(temp_j)


# loading j function values, function of tilde theta!
with np.load('./j_factors/int_dimless_j_fac.npy', 'rb') as infile:
    tildetheta = infile['angle']
    js = infile['js']
    jp = infile['jp']
    jd = infile['jd']
    jsom = infile['jsom']

jsfunc = interpolate.interp1d(tildetheta, js, kind='cubic', fill_value='extrapolate')
jpfunc = interpolate.interp1d(tildetheta, jp, kind='cubic', fill_value='extrapolate')
jdfunc = interpolate.interp1d(tildetheta, jd, kind='cubic', fill_value='extrapolate')
jsomfunc = interpolate.interp1d(tildetheta, jsom, kind='cubic', fill_value='extrapolate')


# dwarf names to compute j_factors for
dwarflist_exclude = ['cetus', 'eridanus2', 'leot', 'and1', 'and3', 'and5', 'and7', 'and14', 'and18', 'tucana3', 'triangulum2', 'segue2', 'hydra2', 'leo4', 'leo5', 'pegasus3', 'pisces2', 'draco2', 'grus1']

dwarflist = np.load('./j_factors/data2_jfac_extra_v1.npy')['name_short']
dwarflist = np.setdiff1d(dwarflist, dwarflist_exclude)

bounds = [0.5, 10]
waves = ['s', 'p', 'd', 'som']

with open('tot_j_factors.txt', 'w') as jvals:
    jvals.write('# median j factors for given wave given as [j-sigma, j, j+sigma]\n')
    jvals.write('# the bounds used are ' + str(bounds) + '\n')
    jvals.write('# dwarf\t\twave\tJ(0.5)\tJ_tot\n')
    for dwarf in dwarflist:
        data = np.load('./j_factors/dwarfs/'+dwarf+'_chain.npy')
        for wave in waves:
            jvals.write(dwarf+'\t'+wave)
            jfac = [np.log10(integrated_j_factor(data, ub, wave)) for ub in bounds]
            jvals.write('\t'+str(jfac)+'\n')
