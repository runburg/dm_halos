#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:51:05 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp
from matplotlib import pyplot as plot
from text import sendtext

# set working precision of decimal points
mp.dps = 30

infile = open("g_d.txt", 'rb')
npzfile = np.load(infile)
g_d = npzfile['gd']
rd = npzfile['r']

gdfunc = interpolate.interp1d(rd, g_d, kind = 'cubic', fill_value='extrapolate')

def hd(y):
#    hint = interpolate.interp1d(rd, np.vectorize(lambda x: gdfunc(x)/np.sqrt(1-((y-0.000000001)/x)**2))(rd), kind = 'cubic', fill_value='extrapolate')
#    return y * integrate.quad(hint, y, rd[-1])[0]
    return y * integrate.quad(lambda x: gdfunc(x)/(mp.sqrt(1-(y/(x))**2)), y, rd[-1])[0]

#l = np.linspace(0, 15, num=500)
#h_s = np.array([h(yy) for yy in l])
#hstot = integrate.quad(h, 0, 5) + integrate.quad(h, 5, 10) + integrate.quad(h, 10, 30) + integrate.quad(h, 30, 60)

nd=integrate.quad(lambda x: x**2 * hd(x), 0, 50)[0]
dd=integrate.quad(hd,0,50)[0]

F_d = mp.sqrt(nd/dd)
sendtext(str(__file__), 'nan')