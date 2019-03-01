#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:50:12 2019

@author: runburg
"""

from scipy import integrate, interpolate
import numpy as np
import mpmath as mp
from matplotlib import pyplot as plot
from text import sendtext

# set working precision of decimal points
mp.dps = 30

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

#l = np.linspace(0, 15, num=500)
#h_s = np.array([h(yy) for yy in l])
#hstot = integrate.quad(h, 0, 5) + integrate.quad(h, 5, 10) + integrate.quad(h, 10, 30) + integrate.quad(h, 30, 60)

pn=integrate.quad(lambda x: x**2 * hp(x), rp[0], 50)[0]
dp=integrate.quad(hp,0,50)[0]

F_p = np.sqrt(pn/dp)
sendtext(str(__file__), "nan")