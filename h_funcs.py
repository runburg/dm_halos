#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:29:06 2019

@author: runburg
"""
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

if __name__ == '__main__':
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

# -----------------------------------------------------------------------------
# THEORETICAL S-WAVE
# -----------------------------------------------------------------------------


def rho(x):
    return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)


def gtfunc(x):
    return rho(x) * rho(x)


def hts(y):
    return y * mp.quad(lambda x: float(gtfunc(x))/(mp.sqrt(1-(y/x)**2)),
                       [y, g['r'][-1]])


# -----------------------------------------------------------------------------
# S-WAVE
# -----------------------------------------------------------------------------


def hs(y):
    # divs gives the divergences (at interpolation points)
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return y * integrate.quad(lambda x: gsfunc(x)/(np.sqrt(1-(y/x)**2)), y,
                              g['r'][-1], points=divs, limit=len(g['r']))[0]


# -----------------------------------------------------------------------------
# P-WAVE
# -----------------------------------------------------------------------------


def hp(y):
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return y * integrate.quad(lambda x: gpfunc(x)/(np.sqrt(1-(y/(x))**2)),
                              y, g['r'][-1], points=divs, limit=len(g['r']))[0]


# -----------------------------------------------------------------------------
# D-WAVE
# -----------------------------------------------------------------------------


def hd(y):
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return y * integrate.quad(lambda x: gdfunc(x)/(np.sqrt(1-(y/(x))**2)),
                              y, g['r'][-1], points=divs, limit=len(g['r']))[0]


# -----------------------------------------------------------------------------
# Sommerfeld-enhancement
# -----------------------------------------------------------------------------


def hsom(y):
    divs = g['r'][np.searchsorted(g['r'], y, side="right"):-2]
    return y * integrate.quad(lambda x: gsomfunc(x)/(np.sqrt(1-(y/x)**2)),
                              y, g['r'][-1], points=divs, limit=len(g['r']))[0]
