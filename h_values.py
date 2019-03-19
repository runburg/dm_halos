#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:21:01 2018

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import matplotlib as mpl
import h_funcs as h
mpl.use('agg')


def hval():
    g = {}
    gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]
    for key in ['g_s', 'g_p', 'g_d', 'g_som']:
        with np.load(key+'.txt', 'rb') as infile:
            g[key] = infile[key]
            g['r'] = infile['r']

    gsfunc = interpolate.interp1d(
        g['r'], g['g_s'], kind='cubic', fill_value='extrapolate')

    gpfunc = interpolate.interp1d(
        g['r'], g['g_p'], kind='cubic', fill_value='extrapolate')

    gdfunc = interpolate.interp1d(
        g['r'], g['g_d'], kind='cubic', fill_value='extrapolate')

    gsomfunc = interpolate.interp1d(
        g['r'], g['g_som'], kind='cubic', fill_value='extrapolate')
    '''h func computation'''
    # normalizations and f_values
    dst = integrate.quad(h.hts, 0, 50)[0]
    nst = integrate.quad(lambda x: x**2*h.hts(x), 0, 50)[0]
    ds = integrate.quad(h.hs, 0, 50)[0]
    ns = integrate.quad(lambda x: x**2*h.hs(x), 0, 50)[0]
    dp = integrate.quad(h.hp, 0, 50)[0]
    nP = integrate.quad(lambda x: x**2*h.hp(x), 0, 50)[0]
    dd = integrate.quad(h.hd, 0, 50)[0]
    nd = integrate.quad(lambda x: x**2*h.hd(x), 0, 50)[0]
    dsom = integrate.quad(h.hsom, 0, 50)[0]
    nsom = integrate.quad(lambda x: x**2*h.hsom(x), 0, 50)[0]

    with open('f_values.txt', 'a') as fvals:
        fvals.write("F-values\n-----------------------------\n")
        fvals.write("Theor. s-wave\t" + str(np.sqrt(nst/dst)))
        fvals.write("\ns-wave\t" + str(np.sqrt(ns/ds)))
        fvals.write("\np_wave\t" + str(np.sqrt(nP/dp)))
        fvals.write("\nd-wave\t" + str(np.sqrt(nd/dd)))
        fvals.write("\nsom.enh.\t" + str(np.sqrt(nsom/dsom)))

    yval = np.logspace(-2.85, 2, num=500)

    # normalized h values
    hst = [h.hts(yy)/dst for yy in yval]
    hs = [h.hs(yy)/ds for yy in yval]
    hp = [h.hp(yy)/dp for yy in yval]
    hd = [h.hd(yy)/dd for yy in yval]
    hsom = [h.hsom(yy)/dsom for yy in yval]

    with open("h_values.txt", 'wb') as outfile:
        np.savez(outfile, hst=np.array(hst), hs=np.array(hs), hp=np.array(
            hp), hd=np.array(hd), hsom=np.array(hsom), radius=np.array(yval))


if __name__ == '__main__':
    hval()
