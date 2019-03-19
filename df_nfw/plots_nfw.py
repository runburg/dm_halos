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

file = 'df_nfw'
g = {}
for key in ['g_s', 'g_p', 'g_d', 'g_som']:
    with np.load(file+'/'+file+'_'+key+'.txt', 'rb') as infile:
        g[key] = infile[key]
        g['r'] = infile['r']

with np.load(file+'/'+file+"_h_values.txt") as infile:
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

o = 280
p = plot.figure()
swave = plot.plot(g['r'][:o], g['g_s'][:o], '-',
                  color=scolor, label=r"$s$-wave")
pwave = plot.plot(g['r'][:o], g['g_p'][:o], '-',
                  color=pcolor, label=r"$p$-wave")
dwave = plot.plot(g['r'][:o], g['g_d'][:o], '-',
                  color=dcolor, label=r"$d$-wave")
som = plot.plot(g['r'][:o], g['g_som'][:o], '-',
                color=somcolor, label=r"Sommerfeld")
plot.xlabel(r"$\tilde{r}$")
plot.ylabel(r"$G(\tilde{r})$")
plot.yscale('log')
plot.ylim(bottom=float(g['g_s'][o])-100, top=float(g['g_som'][0]))
plot.xlim(left=g['r'][0], right=g['r'][o-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_gvalues.pdf", bbox_inches="tight")

# h plot
p = plot.figure()
swave = plot.plot(radius, hsn, '-', color=scolor, label=r"$s$-wave")
pwave = plot.plot(radius, hpn, '-', color=pcolor, label=r"$p$-wave")
dwave = plot.plot(radius, hdn, '-', color=dcolor, label=r"$d$-wave")
som = plot.plot(radius, hsomn, '-', color=somcolor, label=r"Sommerfeld")

plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')
plot.yscale('log')

o = -100
plot.ylim(bottom=float(hdn[o]), top=float(hsomn[0]))
plot.xlim(left=radius[0], right=radius[o])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_hfuncs_log.pdf", bbox_inches='tight')

# h plot end behavior
p = plot.figure()
swave = plot.plot(radius, hsn, '-', color=scolor, label=r"$s$-wave")
pwave = plot.plot(radius, hpn, '-', color=pcolor, label=r"$p$-wave")
dwave = plot.plot(radius, hdn, '-', color=dcolor, label=r"$d$-wave")
som = plot.plot(radius, hsomn, '-', color=somcolor, label=r"Sommerfeld")

plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')
plot.yscale('log')

o = 200
plot.ylim(bottom=float(hdn[-1]), top=float(hdn[o])+1)
plot.xlim(left=radius[o], right=radius[-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_hfuncs_zoom.pdf", bbox_inches='tight')
