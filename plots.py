#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plotting for dm_halos."""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'xx-large'

mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 2.5
mpl.rcParams['xtick.minor.width'] = 0.5

mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 2.5
mpl.rcParams['ytick.minor.width'] = 0.5

mpl.rcParams['legend.fontsize'] = 'x-large'

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

slabel = r"$s$-wave, $n=1$"
plabel = r"$p$-wave, $n=2$"
dlabel = r"$d$-wave, $n=4$"
somlabel = r"Sommerfeld, $n=-1$"
o = 280
p = plot.figure()
swave = plot.plot(g['r'][:o], g['g_s'][:o]/np.sqrt(4*np.pi), '-',
                  color=scolor, label=slabel)
pwave = plot.plot(g['r'][:o], g['g_p'][:o]/np.sqrt((4*np.pi)**2), '-',
                  color=pcolor, label=dlabel)
dwave = plot.plot(g['r'][:o], g['g_d'][:o]/np.sqrt((4*np.pi)**4), '-',
                  color=dcolor, label=plabel)
som = plot.plot(g['r'][:o], g['g_som'][:o]*np.sqrt(4*np.pi), '-',
                color=somcolor, label=somlabel)
plot.xlabel(r"$\tilde{r}$")
plot.ylabel(r"$P^2_n(\tilde{r})$")
plot.xscale('log')
plot.yscale('log')
plot.ylim(bottom=float(g['g_s'][o])-100, top=float(g['g_som'][0]))
plot.xlim(left=g['r'][0], right=g['r'][o-1])
leg = plot.legend(frameon=False, markerscale=50)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_p2values.pdf", bbox_inches="tight")

# h plot
p = plot.figure()
swave = plot.plot(radius, hsn/radius, '-', color=scolor, label=slabel)
pwave = plot.plot(radius, hpn/radius, '-', color=pcolor, label=plabel)
dwave = plot.plot(radius, hdn/radius, '-', color=dcolor, label=dlabel)
som = plot.plot(radius, hsomn/radius, '-', color=somcolor, label=somlabel)

plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')
plot.yscale('log')

o = -100
plot.ylim(bottom=float(hdn[o]), top=float(hsomn[0]))
plot.xlim(left=radius[0], right=radius[o])
leg = plot.legend(frameon=False, markerscale=50, loc=3)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_jfuncs_log.pdf", bbox_inches='tight')

# h plot end behavior
p = plot.figure()
swave = plot.plot(radius, hsn/radius, '-', color=scolor, label=slabel)
pwave = plot.plot(radius, hpn/radius, '-', color=pcolor, label=plabel)
dwave = plot.plot(radius, hdn/radius, '-', color=dcolor, label=dlabel)
som = plot.plot(radius, hsomn/radius, '-', color=somcolor, label=somlabel)

plot.xlabel(r"$\tilde{\theta}$")
plot.ylabel(r"$H_n(\tilde{\theta})/h_n$")
plot.xscale('log')
plot.yscale('log')

o = 205
plot.ylim(bottom=float(hdn[-1]), top=float(hdn[o])+1)
plot.xlim(left=radius[o], right=radius[-1])
leg = plot.legend(frameon=False, markerscale=50, loc=3)
for line in leg.get_lines():
    line.set_linewidth(6)
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)

p.savefig(file+'/'+file+"_jfuncs_zoom.pdf", bbox_inches='tight')
