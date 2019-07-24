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
import matplotlib.pyplot as plt


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


def plot_parameter_values(values, dwarf=None):
    """Plot the parameters as a function of the weights."""
    p = plt.figure()
    plt.plot(values[:, 0], np.divide(values[:, 3], values[:, 2]), lw=.1)
    plt.yscale('log')
    plt.ylabel(r"$r_s/D$")
    plt.xlabel(r"weight")
    p.savefig('param_rs_by_d.pdf', bbox_inches='tight')
    p = plt.figure()
    plt.plot(values[:, 0], values[:, 3], lw=.1)
    plt.yscale('log')
    plt.ylabel(r"$r_s$")
    plt.xlabel(r"weight")
    p.savefig('param_rs.pdf', bbox_inches='tight')
    p = plt.figure()
    plt.plot(values[:, 0], values[:, 4], lw=.1)
    plt.yscale('log')
    plt.ylabel(r"$\rho_s$")
    plt.xlabel(r"weight")
    p.savefig('param_rhos.pdf', bbox_inches='tight')
    p = plt.figure()
    plt.plot(values[:, 0], values[:, 1], lw=.1)
    plt.yscale('log')
    plt.ylabel(r"$J_s$")
    plt.xlabel(r"weight")
    p.savefig('param_js.pdf', bbox_inches='tight')
    p = plt.figure()
    ax1 = plt.subplot(311)
    plt.plot(values[:, 0], values[:, 3], lw=.1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.yscale('log')
    plt.ylabel(r"$r_s$")
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(values[:, 0], values[:, 4], lw=.1)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.yscale('log')
    plt.ylabel(r"$\rho_s$")
    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(values[:, 0], values[:, 1], lw=.1)
    plt.yscale('log')
    plt.ylabel(r"$J_s$")
    plt.xlabel(r"weight")
    p.savefig('param_stacked.pdf', bbox_inches='tight')


def weighted_median(values, plot=False, dwarf=None):
    """Compute the weighted median and return the +/- sigma values."""
    low_jvals, ave_jvals, up_jvals = [], [], []
    temp_j = np.asarray(values)
    weightsum = temp_j[:, 0].sum()
    temp_j = temp_j[temp_j[:, 1].argsort()]
    temp_j[:, 0] = np.cumsum(temp_j, axis=0)[:, 0]
    temp_j = np.divide(temp_j, np.array([weightsum, 1, 1, 1, 1]))
    low_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5-.34, side='left')][1]
    ave_jvals = temp_j[np.searchsorted(temp_j[:, 0], 0.5, side='left')]
    up_jvals = temp_j[np.searchsorted(temp_j[:, 0], .5+.34, side='right')][1]

    if plot:
        # index = np.searchsorted(temp_j[:, 0], 0.5, side='left')
        # print(temp_j[index-10:index+10][:, [0, 1, 3, 4]])
        plot_parameter_values(temp_j, dwarf)

    # if dwarf == 'sagittarius2' or dwarf == 'horologium1':
    #     index = np.searchsorted(temp_j[:, 0], 0.5, side='left')
    #     print(f'{dwarf}: low: {low_jvals} ave: {ave_jvals[1]} high: {up_jvals}')

    return np.array([low_jvals, ave_jvals[1], up_jvals, ave_jvals[2], ave_jvals[3], ave_jvals[4]])


def kpc_to_cm(num):
    """Convert from kpc to cm."""
    return num * float(3.1E+19) * 100


def integrated_j_factor(list, upper_bound, wave, plot=False, dwarf=None):
    """Return the value of the integrated j-factor for wave annihilation to angle upper_bound in degrees."""
    temp_j = []
    # changeunits_to_gevcm = 37.96**2 * 10**(-18) * 1000 * 3.086 * 10**18
    changeunits_to_gevcm = float(3.1E+21)**(-5) * (float(3E+8)**2 * float(1.98E+30) * float(6.26E+9))**2
    for x in list:
        n = 0
        ub = upper_bound/180*np.pi*x[1]/10**x[2]
        if ub > 30:
            ub = 30
        if wave == 's':
            n = 0
            jtot = jsfunc(ub)
            j_tot = 1/3
        if wave == 'p':
            n = 2
            jtot = jpfunc(ub)
            j_tot = 0.15
        if wave == 'd':
            n = 4
            jtot = jdfunc(ub)
            j_tot = 0.12
        if wave == 'som':
            n = -1
            jtot = jsomfunc(ub)
            # don't have a working angular distribution of jsom due to instabilities in the numerical integrsation
            jtot = 1.01

        c = 299792.458  # km/s
        g_n = float(4.303E-6)  # kpc * (km/s)^2 / M_solar
        # g_n = float(6.7E-39)  # GeV^-2
        d = x[1]  # kpc
        r_s = 10**x[2]  # kpc
        rho_s = 10**x[3]  # M_solar / kpc^3
        # if (dwarf == 'sagittarius2' or dwarf == 'horologium1') and wave == 'd':
        #     print('{}:: angle: {} jtot: {}'.format(dwarf, ub, jtot))
        temp_j.append((x[0], 4*np.pi*r_s**3*rho_s**2/(d**2)*(4*np.pi * g_n * rho_s * r_s**2 / c**2)**(n/2.) * jtot * changeunits_to_gevcm, d, r_s, rho_s))
    if len(list) == 1:
        return(temp_j)
    else:
        results = weighted_median(temp_j, plot=plot, dwarf=dwarf)
        # if (dwarf == 'tucana2_des' or dwarf == 'tucana2_k15'):
        #     print('{}::low: {}, ave: {}, up: {}, d: {}, r_s: {}, rho_s: {}'.format(dwarf, *results))
        return results


def check_integrated_js():
    """Return plot of j fucntions."""
    load_funcs()
    x = np.logspace(-3, 1, num=100)
    plt.clf()
    p = plt.figure()
    plt.plot(x, [jsfunc(xx) for xx in x], label='s-wave')
    plt.plot(x, [jpfunc(xx) for xx in x], label='p-wave')
    plt.plot(x, [jdfunc(xx) for xx in x], label='d-wave')
    plt.plot(x, [jsomfunc(xx) for xx in x], label='som enh')
    plt.xscale('log')
    plt.legend()

    p.savefig('integrated_j_factors_by_theta.pdf', bbox_inches='tight')


# loading j function values, function of tilde theta!
# from integrated_j_factors.py
def load_funcs():
    """Create functions for computing j-factors."""
    global jsfunc, jpfunc, jdfunc, jsomfunc
    with np.load('/Users/runburg/github/dm_halos/j_factors/int_dimless_j_fac.npy', 'rb') as infile:
        tildetheta = infile['angle']
        js = infile['js']
        jp = infile['jp']
        jd = infile['jd']
        jsom = infile['jsom']

    jsfunc = interpolate.interp1d(tildetheta, js, kind='cubic', fill_value='extrapolate')
    jpfunc = interpolate.interp1d(tildetheta, jp, kind='cubic', fill_value='extrapolate')
    jdfunc = interpolate.interp1d(tildetheta, jd, kind='cubic', fill_value='extrapolate')
    jsomfunc = interpolate.interp1d(tildetheta, jsom, kind='quadratic', fill_value='extrapolate')


def main():
    """Calculate the median J-factors for each dwarf."""
    load_funcs()
    # dwarf names to compute j_factors for
    dwarflist_exclude = ['cetus', 'eridanus2', 'leot', 'and1', 'and3', 'and5', 'and7', 'and14', 'and18', 'tucana3', 'triangulum2', 'segue2', 'hydra2', 'leo4', 'leo5', 'pegasus3', 'pisces2', 'draco2', 'grus1', 'horologium1_des', 'reticulum2_des', 'tucana2_k15']

    dwarflist = np.load('./j_factors/data2_jfac_extra_v1.npy')['name_short']
    dwarflist = np.concatenate((np.setdiff1d(dwarflist, dwarflist_exclude), np.array(['crater2', 'hydrus1', 'sagittarius2'])))

    bounds = [0.5, 10]
    waves = ['s', 'p', 'd', 'som']

    angledict = {
        0.5: 'half',
        10: '10',
        25: '25',
        50: '50'
    }
    ub = 0.5
    # dwarflist = ['comaberenices']
    for ub in bounds:
        jfac = []
        # dwarflist = dwarflist_exclude
        for dwarf in dwarflist:
            data = np.load('./j_factors/dwarfs/'+dwarf+'_chain.npy')
            for wave in waves:
                plot = False
                # if dwarf == 'comaberenices' and wave == 'som' and ub == 10:
                #     plot = True
                l, a, u, d, r_s, rho_s = integrated_j_factor(data, ub, wave, plot=plot, dwarf=dwarf)
                # if dwarf == 'reticulum2_k15' or dwarf == 'reticulum2_des':
                #     print('{}, {}: {}'.format(dwarf, wave, r_s/d))
                jfac.append([dwarf, wave, np.log10(l), np.log10(a), np.log10(u)])
                # if (dwarf == 'comaberenices' or dwarf == 'draco1' or dwarf == 'reticulum2_k15' or dwarf == 'ursaminor' or dwarf == 'segue1') and wave == 'som' and ub == 10:
                #     print('{}, {}, {}, {}, {}'.format(dwarf, wave, d, r_s, rho_s))

        outfile = "/Users/runburg/github/dm_halos/j_factors/tot_j_factors/tot_j_fac_inclusive_" + angledict[ub]

        np.save(outfile + ".npy", np.array(jfac))
        np.savetxt(outfile + ".txt", np.array(jfac), fmt='%s')


if __name__ == '__main__':
    main()
    # check_integrated_js()

    # load_funcs()
    # crosscheck = np.array([
    #     [[1, 44, np.log10(0.176), np.log10(264052679.5)]],
    #     [[2, 76, np.log10(0.611), np.log10(132340217.3)]],
    #     [[3, 76, np.log10(0.398), np.log10(168172604.8)]],
    #     [[4, 23, np.log10(0.352), np.log10(180388341.3)]],
    #     [[5, 32, np.log10(0.13), np.log10(292493852.4)]]
    # ])
    # for x in crosscheck:
    #     print(integrated_j_factor(x, 10, 'som')[0][1]*2*np.pi*.01)
