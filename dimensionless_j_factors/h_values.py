#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute h and f values for dm_halos project."""

from scipy import integrate
import numpy as np
import matplotlib as mpl
import h_funcs as h

g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def hval():
    """H func computation."""
    h.g_import()
    limit = 10000
    # normalizations and f_values
    ds = integrate.quad(lambda x: x * h.hs(x), 0, 25, limit=limit, epsabs=1e-05, epsrel=1e-05)[0]
    print('ds, ', ds)

    yval = np.logspace(-3, 1.5, num=500)

    hst = [h.hts(yy) for yy in yval]
    hs = [h.hs(yy)for yy in yval]
    hp = [h.hp(yy) for yy in yval]
    hd = [h.hd(yy) for yy in yval]
    hsom = [h.hsom(yy) for yy in yval]
    print('done')

    with open("./total_j_factors/j_factors/j_thetas.txt", 'wb') as outfile:
        np.savez(outfile, jst=np.array(hst), js=np.array(hs), jp=np.array(
            hp), jd=np.array(hd), jsom=np.array(hsom), theta=np.array(yval))

    dst = integrate.quad(lambda x: x * h.hts(x), 0.001, 25, limit=limit)[0]
    dp = integrate.quad(lambda x: x * h.hp(x), 0.001, 25, limit=limit)[0]
    print('dp, ', dp)
    dd = integrate.quad(lambda x: x * h.hd(x), 0.001, 25, limit=limit)[0]
    print('dd, ', dd)
    dsom = integrate.quad(lambda x: x * h.hsom(x), 0.001, 25, limit=limit)[0]
    print('dsom, ', dsom)
    nst = integrate.quad(lambda x: x**2*h.hts(x), 0, 25, limit=limit)[0]
    ns = integrate.quad(lambda x: x**2*h.hs(x), 0, 25, limit=limit)[0]
    nP = integrate.quad(lambda x: x**2*h.hp(x), 0, 25, limit=limit)[0]
    nd = integrate.quad(lambda x: x**2*h.hd(x), 0, 25, limit=limit)[0]
    nsom = integrate.quad(lambda x: x**2*h.hsom(x), 0, 25, limit=limit)[0]

    with open('./dimensionless_j_factors/f_values.txt', 'a') as fvals:
        fvals.write("ang. size\t|\th_n\n-----------------------------\n")
        fvals.write("Theor. s-wave\t" + str(nst/dst) + "\t" + str(dst))
        fvals.write("\ns-wave\t\t" + str(ns/ds) + "\t" + str(ds))
        fvals.write("\np-wave\t\t" + str(nP/dp) + "\t" + str(dp))
        fvals.write("\nd-wave\t\t" + str(nd/dd) + "\t" + str(dd))
        fvals.write("\nsom.enh.\t" + str(nsom/dsom) + "\t" + str(dsom) + "\n")

    # normalized h values


    with open("./dimensionless_j_factors/h_values.txt", 'wb') as outfile:
        np.savez(outfile, hst=np.array(hst)/dst, hs=np.array(hs)/ds, hp=np.array(
            hp), hd=np.array(hd)/dd, hsom=np.array(hsom)/dsom, radius=np.array(yval))

    # hst = [yy * h.hts(yy) for yy in yval]
    # hs = [yy * h.hs(yy) for yy in yval]
    # hp = [yy * h.hp(yy) for yy in yval]
    # hd = [yy * h.hd(yy) for yy in yval]
    # hsom = [yy * h.hsom(yy) for yy in yval]


if __name__ == '__main__':
    hval()
