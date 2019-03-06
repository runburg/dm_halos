#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 21:01:21 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp
from text import sendtext

# set working precision of decimal points
mp.dps = 50

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

#l = np.linspace(0, 15, num=500)
#h_s = np.array([h(yy) for yy in l])
#hstot = integrate.quad(h, 0, 5) + integrate.quad(h, 5, 10) + integrate.quad(h, 10, 30) + integrate.quad(h, 30, 60)

nsom=integrate.quad(lambda x: x**2 * hsom(x), 0, 50)[0]
dsom=integrate.quad(hsom,0,50)[0]

F_som = mp.sqrt(nsom/dsom)
sendtext(str(__file__), 'nan')
