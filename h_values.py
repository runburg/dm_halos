#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:21:01 2018

@author: runburg
"""

from scipy import integrate
import numpy as np
import matplotlib as mpl
import h_funcs as h
mpl.use('agg')

g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def hval():
    h.g_import()
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
        fvals.write("F-values\t|\th_n\n-----------------------------\n")
        fvals.write("Theor. s-wave\t" + str(np.sqrt(nst/dst)) + "\t" + str(dst))
        fvals.write("\ns-wave\t\t" + str(np.sqrt(ns/ds)) + "\t" + str(ds))
        fvals.write("\np-wave\t\t" + str(np.sqrt(nP/dp)) + "\t" + str(dp))
        fvals.write("\nd-wave\t\t" + str(np.sqrt(nd/dd)) + "\t" + str(dd))
        fvals.write("\nsom.enh.\t" + str(np.sqrt(nsom/dsom))
                    + "\t" + str(dsom) + "\n")

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
