#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute dimensionless j values as a function of tilde theta for dm_halos project."""

import numpy as np
from scipy import integrate
import matplotlib as mpl
import h_funcs as h
import matplotlib.pyplot as plt
mpl.use('agg')

g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def jval():
    """H func computation."""
    h.g_import()
    # yval = np.logspace(-4, 2, num=1000)

    # unnormalized h values
    # hst = [h.hts(yy) for yy in yval]
    # hs = [h.hs(yy) for yy in yval]
    # hp = [h.hp(yy) for yy in yval]
    # hd = [h.hd(yy) for yy in yval]
    # hsom = [h.hsom(yy) for yy in yval]
    #
    # with open("./j_factors/j_thetas.txt", 'wb') as outfile:
    #     np.savez(outfile, jst=np.array(hst), js=np.array(hs), jp=np.array(
    #         hp), jd=np.array(hd), jsom=np.array(hsom), theta=np.array(yval))

    # limitvals = np.logspace(0, 4, num=8)
    # for limit in limitvals:
    #     print('limit={}, result={}'.format(limit, integrate.quad(lambda x: x * h.hsom(x), 0, 50, limit=int(limit))[0]))
    theta = np.logspace(-1, 1.5, num=10)
    plt.plot(theta, [integrate.quad(lambda x: x * h.hsom(x), 0, angle, limit=50)[0] for angle in theta])
    plt.show()

if __name__ == '__main__':
    jval()
