#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:21:01 2018

@author: runburg
"""

from scipy import integrate
from scipy import interpolate
import numpy as np
import mpmath as mp
from matplotlib import pyplot as plot
import hfuncs as h
import send_textmessage


#infile = open("fine_detail_jS_arrays.txt", 'rb')
#npzfile = np.load(infile)
#j_tot_s = npzfile['j_tot_s']
#var = npzfile['j_s_var']
#D = npzfile['D']
#infile.close()
#sd = np.sqrt(var)
#b = 1/D
#F = sd/b
#
#k=1000
#l=50

infile = open("g_p.txt", 'rb')
npzfile = np.load(infile)
g_p = npzfile['gp']
r = npzfile['r']
infile.close()

infile = open("g_s.txt", 'rb')
npzfile = np.load(infile)
g_s = npzfile['gs']
infile.close()

infile = open("g_d.txt", 'rb')
npzfile = np.load(infile)
g_d = npzfile['gd']
infile.close()

infile = open("g_som.txt", 'rb')
npzfile = np.load(infile)
g_som = npzfile['gsom']
infile.close()


def rho2(x):
    return (1/(float(x)*(1+float(x))**2))**2

vrho2 = np.vectorize(rho2)

o=250
p = plot.figure()
theor = plot.plot(r[:o], vrho2(r[:o]), 'm-', markersize = .2, label = r"$s$-wave")
#swave = plot.plot(r[:o], g_s[:o], 'r--', markersize = .2, label = r"$s$-wave")
pwave = plot.plot(r[:o], g_p[:o], 'g-', markersize = .1, label = r"$p$-wave")
dwave = plot.plot(r[:o], g_d[:o], 'b-', markersize = .2, label = r"$d$-wave")
som = plot.plot(r[:o], g_som[:o], 'y-', markersize = .2, label = r"Sommerfeld")
plot.legend(markerscale=25)
plot.xlabel(r"$\tilde{r}$")
plot.ylabel(r"$G(\tilde{r})$")
plot.yscale('log')

p.savefig("gvalues.pdf", bbox_inches="tight")


'''h func plotting'''
dst = integrate.quad(h.hts, 0, 50)[0]
ds=integrate.quad(h.hs,0,50)[0]
dp=integrate.quad(h.hp,0,50)[0]
dd=integrate.quad(h.hd,0,50)[0]
dsom = integrate.quad(h.hsom,0,50)[0]

yval = np.logspace(-2.85, 0.4, num = 250)

isst = [h.hts(yy)/dst for yy in yval]
#iss = [h.hs(yy)/ds for yy in yval]
ipp = [h.hp(yy)/dp for yy in yval]
idd = [h.hd(yy)/dd for yy in yval]
isom = [h.hsom(yy)/dsom for yy in yval]

p = plot.figure()

swavet = plot.plot(yval, isst, 'm-', markersize = .2, label = r"$s$-wave")
#swave = plot.plot(yval, iss, 'r--', markersize = .2, label = r"$s$-wave")
pwave = plot.plot(yval, ipp, 'g-', markersize = .2, label = r"$p$-wave")
dwave = plot.plot(yval, idd, 'b-', markersize = .2, label = r"$d$-wave")
som = plot.plot(yval, isom, 'y-', markersize = .2, label = r"Sommerfeld")
plot.xlabel(r"$y$")
plot.ylabel(r"$H(y)$")
plot.xscale('log')
plot.legend(markerscale=25)

p.savefig("hfuncs.pdf", bbox_inches='tight')
text.sendtext(str(__file__), "hfuncs.pdf")



''''''



#o=250
#p = plot.figure()
#pwave = plot.plot(r[:o], g_p[:o], 'r+', markersize = .2, label = r"$G_p$")
#
#plot.legend(markerscale=25)
#plot.xlabel(r"$\tilde{r}$")
#plot.ylabel(r"$G_p(\tilde{r})$")
#
#plot.show()
#
#p.savefig("gpvalues.pdf", bbox_inches="tight")

#p=plot.figure()
#plot.plot(b[:l], sd[:l], marker="+")
#plot.xlabel("β")
#plot.ylabel(r"$\sqrt{\langle\theta^2\rangle_\beta}$")
#plot.show()
##
#p.savefig("swavevariance.pdf", bbox_inches="tight")
#
#p=plot.figure()
#plot.plot(b[:k], F[:k], marker = "+")
#plot.xlabel("β")
#plot.ylabel(r"$\sqrt{\langle\theta^2\rangle_\beta}\,/ \beta$")
#plot.show()
##
#p.savefig("swavevariancesqrt.pdf", bbox_inches="tight")
#
#theor = []
#comp = []
#with open('fine_detail_jS_results.txt', 'r') as f:
#    content = f.readlines()
#    for x in content:
#        row = x.split()
#        theor.append(row[1])
#        comp.append(row[2])
#
#theor.pop(0)
#theor = np.array(theor)
#theor = theor.astype(np.float)
#p=plot.figure()
#plot.plot(1/D, j_tot_s, linewidth=0.5)
#plot.plot(1/D, theor, linewidth=0.5)
#plot.xlabel("β")
#plot.ylabel(r"Total effective J-factor, $\overline{J}_S$")
#plot.show()
#
#p.savefig("swavejtotscomp.pdf", bbox_inches="tight")

#infile = open("dm_halos_jS_nparrays.txt", 'rb')
#npzfile = np.load(infile)
#j_s = npzfile['j_s']
#theta = npzfile['theta']
#
#o=1000
#p=plot.figure()
#plot.plot(theta[100:o], j_s[100:o], 'ob')
#plot.xlabel(r"$\theta$")
#plot.ylabel(r"$J_S$")
#plot.show()
#
#p.savefig("jfactortheta.pdf", bbox_inches="tight")
