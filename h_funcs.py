#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Define h functs for dm_halos project."""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp

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
        with np.load('./df_nfw/df_nfw_'+key+'.txt', 'rb') as infile:
            g[key] = infile[key]
            g['r'] = infile['r']

    gsfunc = interpolate.interp1d(
        g['r'], g['g_s'], kind='cubic', fill_value='extrapolate')

    gpfunc = interpolate.interp1d(
        g['r'], g['g_p'], kind='cubic', fill_value='extrapolate')

    gdfunc = interpolate.interp1d(
        g['r'], g['g_d'], kind='cubic', fill_value='extrapolate')

    g['g_som'][0] = g['g_som'][1]+g['g_som'][2]
    gsomfunc = interpolate.interp1d(
        g['r'], g['g_som'], kind='cubic', fill_value='extrapolate')


if __name__ == '__main__':
    g_import()


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


# -----------------------------------------------------------------------------
# S-WAVE
# -----------------------------------------------------------------------------


def hs(y):
    """Return the computed h_s."""
    # divs gives the divergences (at interpolation points)
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return integrate.quad(lambda x: gsfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']))[0]


# -----------------------------------------------------------------------------
# P-WAVE
# -----------------------------------------------------------------------------


def hp(y):
    """Return h_p."""
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return integrate.quad(lambda x: gpfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']))[0]


# -----------------------------------------------------------------------------
# D-WAVE
# -----------------------------------------------------------------------------


def hd(y):
    """Return h_d."""
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return integrate.quad(lambda x: gdfunc(x)/(np.sqrt(1-(y/(x))**2)), y*1.00000001, g['r'][-1], points=divs, limit=100*len(g['r']))[0]


# -----------------------------------------------------------------------------
# Sommerfeld-enhancement
# -----------------------------------------------------------------------------


def hsom(y):
    """Return h_som."""
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return integrate.quad(lambda x: gsomfunc(x)/(np.sqrt(1-(y/x)**2)), y*1.00000000001, g['r'][-1], points=divs, limit=10*len(g['r']))[0]
