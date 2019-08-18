#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Define h functs for dm_halos project."""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

# requires g_*.txt
# run g_*.py before hand
#
# defines the h functions given the computed g-gvalues
# from g_*.py in g_*.txt
# functional forms and motivations given in SubstructureNotes
# in ~/Desktop/dm_halos/substructure notes

# exact precision of mp functions
mp.dps = 50

# import relevant files
g = {}
gsfunc, gpfunc, gdfunc, gsomfunc = [0, 0, 0, 0]


def g_import():
    """Globally import all g functions."""
    global gsfunc, gpfunc, gdfunc, gsomfunc
    for key in ['g_s', 'g_p', 'g_d', 'g_som']:
        with np.load('./dimensionless_j_factors/df_nfw/df_nfw_'+key+'.txt', 'rb', allow_pickle=True) as infile:
            g[key] = infile[key].astype(np.float)
            g['r'] = infile['r'].astype(np.float)
        # g[key] = np.insert(g[key], 0, g[key][0]+g[key][1])

    # g['r'] = np.insert(g['r'], 0, 0)

    # g['g_s'][0] = g['g_s'][1]+g['g_s'][2]
    gsfunc = interpolate.interp1d(
        g['r'], g['g_s'], kind='cubic', fill_value='extrapolate')
    print(gsfunc(100))

    gpfunc = interpolate.interp1d(
        g['r'], g['g_p'], kind='cubic', fill_value='extrapolate')

    gdfunc = interpolate.interp1d(
        g['r'], g['g_d'], kind='cubic', fill_value='extrapolate')

    # g['g_som'][0] = g['g_som'][1]+g['g_som'][2]
    gsomfunc = interpolate.interp1d(
        g['r'], g['g_som'], kind='cubic', fill_value='extrapolate')


if __name__ == '__main__':
    g_import()
    plt.plot(g['r'], g['g_s'])
    # plt.plot(g['r'], g['g_p'])
    # plt.plot(g['r'], g['g_d'])
    # plt.plot(g['r'], g['g_som'])
    plt.xscale('log')
    plt.show()

# -----------------------------------------------------------------------------
# THEORETICAL S-WAVE
# -----------------------------------------------------------------------------


def rho(x):
    """Return the density function value."""
    return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)


def gtfunc(x):
    """Return the theoretical g_s."""
    return rho(x) * rho(x)


def hts(y):
    """Return the theoretical h_s."""
    return mp.quad(lambda x: float(gtfunc(x))/(mp.sqrt(1-(y/x)**2)), [y,  g['r'][-1]])


def h_integral(g_function, theta, divs=None):
    return integrate.quad(lambda x: g_function(x)/(np.sqrt(1-(theta/x)**2)), theta*1.000001, 35, limit=10000, points=divs, epsabs=1e-05, epsrel=1e-05)[0]


# -----------------------------------------------------------------------------
# S-WAVE
# -----------------------------------------------------------------------------


def hs(y):
    """Return the computed h_s."""
    # divs gives the divergences (at interpolation points)
    # divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    # return integrate.quad(lambda x: gsfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']), epsabs=1e-03, epsrel=1e-03)[0]
    return h_integral(gsfunc, y)


# -----------------------------------------------------------------------------
# P-WAVE
# -----------------------------------------------------------------------------


def hp(y):
    """Return h_p."""
    # divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    # return integrate.quad(lambda x: gpfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']), epsabs=1e-03, epsrel=1e-03)[0]
    return h_integral(gpfunc, y)


# -----------------------------------------------------------------------------
# D-WAVE
# -----------------------------------------------------------------------------


def hd(y):
    """Return h_d."""
    # divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    # return integrate.quad(lambda x: gdfunc(x)/(np.sqrt(1-(y/(x))**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']), epsabs=1e-03, epsrel=1e-03)[0]
    return h_integral(gdfunc, y)


# -----------------------------------------------------------------------------
# Sommerfeld-enhancement
# -----------------------------------------------------------------------------


def hsom(y):
    """Return h_som."""
    # divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    # return integrate.quad(lambda x: gsomfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000000001, g['r'][-1], points=divs, limit=100*len(g['r']), epsabs=1e-03, epsrel=1e-03)[0]
    return h_integral(gsomfunc, y)
