#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 12 2019

@author: runburg
"""

from scipy import integrate
import numpy as np
import matplotlib as mpl
import hfuncs as h
mpl.use('agg')

yval = np.logspace(-1, 2, num=200)

dst = integrate.quad(h.hts, 0, 50)[0]
# nst = integrate.quad(lambda x: x**2*h.hts(x), 0, 50)[0]
ds = integrate.quad(h.hs, 0, 50)[0]
# ns = integrate.quad(lambda x: x**2*h.hs(x), 0, 50)[0]
dp = integrate.quad(h.hp, 0, 50)[0]
# nP = integrate.quad(lambda x: x**2*h.hp(x), 0, 50)[0]
dd = integrate.quad(h.hd, 0, 50)[0]
# nd = integrate.quad(lambda x: x**2*h.hd(x), 0, 50)[0]
dsom = integrate.quad(h.hsom, 0, 50)[0]

isst = [h.hts(yy)/dst for yy in yval]
# iss = [h.hs(yy)/ds for yy in yval]
ipp = [h.hp(yy)/dp for yy in yval]
idd = [h.hd(yy)/dd for yy in yval]
isom = [h.hsom(yy)/dsom for yy in yval]

with open("h_funcs_zoom.txt", 'wb') as outfile:
    np.savez(outfile, radius=yval, hsn=np.array(isst), hpn=np.array(
        ipp), hdn=np.array(idd), hdsomn=np.array(isom))
