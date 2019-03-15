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

# -----------------------------------------------------------------------------
# THEORETICAL S-WAVE
# -----------------------------------------------------------------------------


def rho(x):
    return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)


def gtfunc(x):
    return rho(x) * rho(x)


def hts(y):
    return y * mp.quad(lambda x: float(gtfunc(x))/(mp.sqrt(1-(y/x)**2)),
                       [y, rs[-1]])


# -----------------------------------------------------------------------------
# S-WAVE
# -----------------------------------------------------------------------------

infile = open("g_s.txt", 'rb')
npzfile = np.load(infile)
g_s = npzfile['gs']
rs = npzfile['r']

gsfunc = interpolate.interp1d(rs, g_s, kind='cubic', fill_value='extrapolate')


def hs(y):
    # divs gives the divergences (at interpolation points)
    divs = rs[np.searchsorted(rs, y, side="right"):-2]
    return y * integrate.quad(lambda x: gsfunc(x)/(np.sqrt(1-(y/x)**2)), y,
                              rs[-1], points=divs, limit=len(rs))[0]


# -----------------------------------------------------------------------------
# P-WAVE
# -----------------------------------------------------------------------------

infile = open("g_p.txt", 'rb')
npzfile = np.load(infile)
g_p = npzfile['gp']
rp = npzfile['r']

gpfunc = interpolate.interp1d(rp, g_p, kind='cubic', fill_value='extrapolate')


def hp(y):
    divs = rp[np.searchsorted(rp, y, side="right"):-2]
    return y * integrate.quad(lambda x: gpfunc(x)/(np.sqrt(1-(y/(x))**2)),
                              y, rp[-1], points=divs, limit=len(rp))[0]


# -----------------------------------------------------------------------------
# D-WAVE
# -----------------------------------------------------------------------------

infile = open("g_d.txt", 'rb')
npzfile = np.load(infile)
g_d = npzfile['gd']
rd = npzfile['r']

gdfunc = interpolate.interp1d(rd, g_d, kind='cubic', fill_value='extrapolate')


def hd(y):
    divs = rd[np.searchsorted(rd, y, side="right"):-2]
    return y * integrate.quad(lambda x: gdfunc(x)/(np.sqrt(1-(y/(x))**2)),
                              y, rd[-1], points=divs, limit=len(rp))[0]


# -----------------------------------------------------------------------------
# Sommerfeld-enhancement
# -----------------------------------------------------------------------------

infile = open("g_som.txt", 'rb')
npzfile = np.load(infile)
g_som = npzfile['gsom']
rsom = npzfile['r']

g_som = g_som * mp.pi/2
gsomfunc = interpolate.interp1d(
    rsom, g_som, kind='cubic', fill_value='extrapolate')


def hsom(y):
    divs = rd[np.searchsorted(rsom, y, side="right"):-2]
    return y * integrate.quad(lambda x: gsomfunc(x)/(np.sqrt(1-(y/x)**2)),
                              y, rsom[-1], points=divs, limit=len(rp))[0]
