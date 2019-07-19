#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute h and f values for dm_halos project."""

from scipy import integrate
import numpy as np
import matplotlib as mpl
import h_funcs as h
mpl.use('agg')

g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def hval():
    """H func computation."""
    h.g_import()
    # normalizations and f_values
    limit = 1000
    dst = integrate.quad(lambda x: x * h.hts(x), 0, 50, limit=limit)[0]
    nst = integrate.quad(lambda x: x**2*h.hts(x), 0, 50, limit=limit)[0]
    ds = integrate.quad(lambda x: x * h.hs(x), 0, 50, limit=limit)[0]
    ns = integrate.quad(lambda x: x**2*h.hs(x), 0, 50, limit=limit)[0]
    dp = integrate.quad(lambda x: x * h.hp(x), 0, 50, limit=limit)[0]
    nP = integrate.quad(lambda x: x**2*h.hp(x), 0, 50, limit=limit)[0]
    dd = integrate.quad(lambda x: x * h.hd(x), 0, 50, limit=limit)[0]
    nd = integrate.quad(lambda x: x**2*h.hd(x), 0, 50, limit=limit)[0]
    dsom = integrate.quad(lambda x: x * h.hsom(x), 0, 50, limit=limit)[0]
    nsom = integrate.quad(lambda x: x**2*h.hsom(x), 0, 50, limit=limit)[0]

    with open('f_values.txt', 'a') as fvals:
        fvals.write("ang. size\t|\th_n\n-----------------------------\n")
        fvals.write("Theor. s-wave\t" + str(nst/dst) + "\t" + str(dst))
        fvals.write("\ns-wave\t\t" + str(ns/ds) + "\t" + str(ds))
        fvals.write("\np-wave\t\t" + str(nP/dp) + "\t" + str(dp))
        fvals.write("\nd-wave\t\t" + str(nd/dd) + "\t" + str(dd))
        fvals.write("\nsom.enh.\t" + str(nsom/dsom) + "\t" + str(dsom) + "\n")

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

    # hst = [yy * h.hts(yy) for yy in yval]
    # hs = [yy * h.hs(yy) for yy in yval]
    # hp = [yy * h.hp(yy) for yy in yval]
    # hd = [yy * h.hd(yy) for yy in yval]
    # hsom = [yy * h.hsom(yy) for yy in yval]

    with open("j_thetas.txt", 'wb') as outfile:
        np.savez(outfile, jst=np.array(hst)*dst, js=np.array(hs)*hs, jp=np.array(
            hp)*dp, jd=np.array(hd)*dd, jsom=np.array(hsom)*dsom, theta=np.array(yval))


if __name__ == '__main__':
    hval()
