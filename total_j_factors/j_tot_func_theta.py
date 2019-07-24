#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute j tot as a function of theta for dm_halos project."""

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import h_funcs as h

g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def jval():
    """H func computation."""
    h.g_import()
    # normalizations and f_values
    limit = 1000
    yval = np.logspace(-1, np.log10(50), num=20)
    jst = [integrate.quad(lambda x: x * h.hts(x), 0, y, limit=limit)[0] for y in yval]
    js = [integrate.quad(lambda x: x * h.hs(x), 0, y, limit=limit)[0] for y in yval]
    jp = [integrate.quad(lambda x: x * h.hp(x), 0, y, limit=limit)[0] for y in yval]
    jd = [integrate.quad(lambda x: x * h.hd(x), 0, y, limit=limit)[0] for y in yval]
    jsom = [integrate.quad(lambda x: x * h.hsom(x), 0, y, limit=limit)[0] for y in yval]

    with open("j_func_thetas.txt", 'wb') as outfile:
        np.savez(outfile, jst=np.array(jst), js=np.array(js), jp=np.array(
            jp), jd=np.array(jd), jsom=np.array(jsom), theta=np.array(yval))


if __name__ == '__main__':
    # jval()
    with np.load('j_func_thetas.txt', 'rb') as infile:
        theta = infile['theta']
        js = infile['js']
        jp = infile['jp']
        jd = infile['jd']
        jsom = infile['jsom']
        # plt.plot(theta, js)
        # plt.plot(theta, jp)
        # plt.plot(theta, jd)
    h.g_import()
    # yval = np.logspace(-1, 1, num=20)
    # plt.plot(theta, [h.hsom(y) for y in yval])
    # plt.show()
    print([(y, h.hsom(y)) for y in np.logspace(-4, 1, num=50)])
    print(integrate.quad(lambda x: x * h.hsom(x), 0, 0.1, limit=50))
    print(integrate.quad(lambda x: x * h.hsom(x), 0, 1, limit=50))
