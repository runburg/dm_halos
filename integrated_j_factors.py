#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Return integrated J-factors as a function of angle.

Author: Jack Runburg
Date: 07-06-2019 11:47


"""

from scipy import interpolate, integrate
import numpy as np


def dimensionless_integrated_j_factor():
    """Integrates dimensionless angular distribution and returns it as a function of angle."""
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

    jtot_s, jtot_p, jtot_d, jtot_som = [], [], [], []
    ub = np.logspace(-3, 1, num=100)
    for angle in ub:
        jtot_s.append(integrate.quad(lambda y: y * jsfunc(y), 0, angle, full_output=1)[0])
        jtot_p.append(integrate.quad(lambda y: y * jpfunc(y), 0, angle, full_output=1)[0])
        jtot_d.append(integrate.quad(lambda y: y * jdfunc(y), 0, angle, full_output=1)[0])
        jtot_som.append(integrate.quad(lambda y: y * jsomfunc(y), 0, angle, full_output=1)[0])

    return np.array([ub, jtot_s, jtot_p, jtot_d, jtot_som])


if __name__ == '__main__':
    funcs = dimensionless_integrated_j_factor()
    with open('int_dimless_j_fac.npy', 'wb') as outfile:
        np.savez(outfile, angle=funcs[0], js=funcs[1], jp=funcs[2], jd=funcs[3], jsom=funcs[4])