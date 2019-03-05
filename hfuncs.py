#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:29:06 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp

mp.dps = 100

def rho(x):
    return 1/(mp.mpf(x)*(1+mp.mpf(x))**2)

gtfunc = lambda x: rho(x) * rho(x)

def hts(y):
    return y * integrate.quad(lambda x: gtfunc(x)/(mp.sqrt(1-(y/x)**2)), y, rs[-1])[0]

infile = open("g_s.txt", 'rb')
npzfile = np.load(infile)
g_s = npzfile['gs']
rs = npzfile['r']

gfunc = interpolate.interp1d(rs, g_s, kind = 'cubic', fill_value='extrapolate')

def hs(y):
#    hsint = interpolate.interp1d(rs, np.vectorize(lambda x: gfunc(x)/np.sqrt(1.0000000001-(y/x)**2))(rs), kind = 'cubic')
#    return y * integrate.quad(interpolate.interp1d(rs, np.vectorize(lambda x: gfunc(x)/(mp.re(mp.sqrt(1.00000000001-(y/(x))**2)+0.00000000000000001)))(rs), kind = 'cubic'), y, rs[-1])[0]
    return y * integrate.quad(lambda x: gfunc(x)/(mp.sqrt(1-(y/(x))**2)), y, rs[-1])[0]

infile = open("g_p.txt", 'rb')
npzfile = np.load(infile)
g_p = npzfile['gp']
rp = npzfile['r']

gpfunc = interpolate.interp1d(rp, g_p, kind = 'cubic', fill_value='extrapolate')

def hp(y):
#    hpints = [gpfunc(x)/np.sqrt(1-(y/x)**2) for x in rp]
#    hpint = interpolate.interp1d(rp, hpints, kind = 'cubic')
#    hpint = interpolate.interp1d(rp[:-1], np.vectorize(lambda x: gpfunc(x)/np.sqrt(1-((y-.00000001)/x)**2))(rp[:-1]), kind = 'cubic', fill_value = 'extrapolate')
#    return y * integrate.quad(hpint, y, rp[-1])[0]
    return y * integrate.quad(lambda x: gpfunc(x)/(mp.sqrt(1-(y/(x))**2)), y, rp[-1])[0]

infile = open("g_d.txt", 'rb')
npzfile = np.load(infile)
g_d = npzfile['gd']
rd = npzfile['r']

gdfunc = interpolate.interp1d(rd, g_d, kind = 'cubic', fill_value='extrapolate')

def hd(y):
#    hint = interpolate.interp1d(rd, np.vectorize(lambda x: gdfunc(x)/np.sqrt(1-((y-0.000000001)/x)**2))(rd), kind = 'cubic', fill_value='extrapolate')
#    return y * integrate.quad(hint, y, rd[-1])[0]
    return y * integrate.quad(lambda x: gdfunc(x)/(mp.sqrt(1-(y/(x))**2)), y, rd[-1])[0]

infile = open("g_som.txt", 'rb')
npzfile = np.load(infile)
g_som = npzfile['gsom']
rsom = npzfile['r']

g_som = g_som * mp.pi/2
gsomfunc = interpolate.interp1d(rsom, g_som, kind = 'cubic', fill_value='extrapolate')

def hsom(y):
#    hint = interpolate.interp1d(rd, np.vectorize(lambda x: gdfunc(x)/np.sqrt(1-((y-0.000000001)/x)**2))(rd), kind = 'cubic', fill_value='extrapolate')
#    return y * integrate.quad(hint, y, rd[-1])[0]
    return y * integrate.quad(lambda x: gsomfunc(x)/(mp.sqrt(1-(mp.fdiv(y,x))**2)), y, rsom[-1])[0]
