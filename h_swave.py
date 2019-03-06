#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 11:59:45 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp
from text import sendtext

# set working precision of decimal points
mp.dps = 30

infile = open("g_s.txt", 'rb')
npzfile = np.load(infile)
g_s = npzfile['gs']
rs = npzfile['r']

gfunc = interpolate.interp1d(rs, g_s, kind = 'cubic', fill_value='extrapolate')

def hs(y):
#    hsint = interpolate.interp1d(rs, np.vectorize(lambda x: gfunc(x)/np.sqrt(1.0000000001-(y/x)**2))(rs), kind = 'cubic')
#    return y * integrate.quad(interpolate.interp1d(rs, np.vectorize(lambda x: gfunc(x)/(mp.re(mp.sqrt(1.00000000001-(y/(x))**2)+0.00000000000000001)))(rs), kind = 'cubic'), y, rs[-1])[0]
    return y * integrate.quad(lambda x: mp.fdiv(gfunc(x),mp.sqrt(1-(y/x)**2)), y, rs[-1])[0]

#l = np.linspace(0, 15, num=500)
#h_s = np.array([h(yy) for yy in l])
#hstot = integrate.quad(h, 0, 5) + integrate.quad(h, 5, 10) + integrate.quad(h, 10, 30) + integrate.quad(h, 30, 60)

#ns=integrate.quad(lambda x: x**2 * hs(x), 0, 50)[0]
#ds=integrate.quad(hs,0,50)[0]
#
#F_s = mp.sqrt(ns/ds)
sendtext(str(__file__), "nan")
