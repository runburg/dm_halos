#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11, 2019 15:18 HST

@author: runburg
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot

mpl.rcParams['lines.linewidth'] = 1
plot.rc('text', usetex=True)
plot.rc('font', family='serif')

with np.load("g_s.txt") as infile:
    g_s = infile['g_s']
    r = infile['r']

with np.load("g_p.txt") as infile:
    g_p = infile['g_p']
    r = infile['r']

with np.load("g_d.txt") as infile:
    g_d = infile['g_d']
    r = infile['r']

with np.load("g_som.txt") as infile:
    g_som = infile['g_som']
    r = infile['r']

with np.load("h_values.txt") as infile:
    hsn = infile['hs']
    hpn = infile['hp']
    hdn = infile['hd']
    hsomn = infile['hsom']
    radius = infile['radius']


# g plot
scolor = 'xkcd:azure'
pcolor = 'xkcd:coral'
dcolor = 'xkcd:peach'
somcolor = 'xkcd:light turquoise'

o = 250
p = plot.figure()
swave = plot.plot(r[:o], g_s[:o], '-',
                  color=scolor, label=r"$s$-wave")
pwave = plot.plot(r[:o], g_p[:o], '-',
                  color=pcolor, label=r"$p$-wave")
dwave = plot.plot(r[:o], g_d[:o], '-', color=dcolor, label=r"$d$-wave")
som = plot.plot(r[:o], g_som[:o], '-', color=somcolor, label=r"Sommerfeld")
plot.xlabel(r"$\tilde{r}$")
plot.ylabel(r"$G(\tilde{r})$")
plot.yscale('log')
plot.ylim(bottom=float(g_som[o]), top=float(g_d[0]))
plot.xlim(left=r[0], right=r[o-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig("gvalues.pdf", bbox_inches="tight")

# h plot
p = plot.figure()
swave = plot.plot(radius, hsn, '-', color=scolor, label=r"$s$-wave")
pwave = plot.plot(radius, hpn, '-', color=pcolor, label=r"$p$-wave")
dwave = plot.plot(radius, hdn, '-', color=dcolor, label=r"$d$-wave")
som = plot.plot(radius, hsomn, '-', color=somcolor, label=r"Sommerfeld")

plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')

plot.ylim(bottom=float(hdn[-1])-0.04, top=float(hsomn[0]))
plot.xlim(left=radius[0], right=radius[-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig("hfuncs.pdf", bbox_inches='tight')

# h plot end behavior

swave = plot.plot(radius, hsn, '-', color=scolor, label=r"$s$-wave")
pwave = plot.plot(radius, hpn, '-', color=pcolor, label=r"$p$-wave")
dwave = plot.plot(radius, hdn, '-', color=dcolor, label=r"$d$-wave")
som = plot.plot(radius, hsomn, '-', color=somcolor, label=r"Sommerfeld")
plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')
plot.yscale('log')

plot.ylim(bottom=float(hdn[-1]), top=float(hdn[0]))
plot.xlim(left=radius[0], right=radius[-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig("hfuncs_zoom.pdf", bbox_inches='tight')
