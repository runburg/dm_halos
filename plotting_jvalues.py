#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plotting j-value averages.

Author: Jack Runburg
Date: 21-05-2019 12:10


"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot

# Plotting parameters
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'xx-large'

mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 7.5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3.75
mpl.rcParams['xtick.minor.width'] = 0.5

mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.size'] = 7.5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3.75
mpl.rcParams['ytick.minor.width'] = 0.5

mpl.rcParams['legend.fontsize'] = 'x-large'

plot.rc('text', usetex=True)
plot.rc('font', family='serif')

dwarflist = ["draco", "segue"]

colordict = {
    'draco': 'xkcd:azure',
    'segue': 'xkcd:coral'
    }

labeldict = {
    'draco': r"Draco",
    'segue': r"Segue 1"
}

wavedict = {
    's': 'xkcd:azure',
    'p': 'xkcd:coral',
    'd': 'xkcd:peach',
    'som': 'xkcd:light turquoise'
}

wavelabeldict = {
    's': r"$s$-wave, $n=0$",
    'p': r"$s$-wave, $n=2$",
    'd': r"$d$-wave, $n=4$",
    'som': r"Sommerfeld, $n=-1$"
}
# dcolor = 'xkcd:peach'
# somcolor = 'xkcd:light turquoise'

# dlabel = r"$d$-wave, $n=4$"
# somlabel = r"Sommerfeld, $n=-1$"

for wave in ['s', 'p', 'd', 'som']:
    ymax = 0
    ymin = 1000000000000

    ave_jvals, up_jvals, low_jvals, thetas = [], [], [], []
    p = plot.figure()
    for dwarf in dwarflist:
        with np.load("./j_factors/j_ave/j_ave_thetas_"+dwarf+"_"+wave+".txt", 'rb') as infile:
            ave_jvals = infile['ave_jvals']
            up_jvals = infile['up_jvals']
            low_jvals = infile['low_jvals']
            thetas = infile['thetas']

        thetas = np.array(thetas*180/np.pi, dtype=float)
        ave_jvals = np.array(ave_jvals, dtype=float)
        up_jvals = np.array(up_jvals, dtype=float)
        low_jvals = np.array(low_jvals, dtype=float)

        if up_jvals[0] > ymax:
            ymax = up_jvals[0]
        if low_jvals[-1] < ymin:
            ymin = low_jvals[-1]

        plot.plot(thetas, ave_jvals, '-', color=colordict[dwarf], label=labeldict[dwarf])
        plot.plot(thetas, up_jvals, '-', linewidth=.01, color=colordict[dwarf])
        plot.plot(thetas, low_jvals, '-', linewidth=.01, color=colordict[dwarf])
        plot.fill_between(thetas, up_jvals, ave_jvals, facecolor=colordict[dwarf], alpha=.4)
        plot.fill_between(thetas, ave_jvals, low_jvals, facecolor=colordict[dwarf], alpha=.4)

    plot.xlabel(r"$\theta[^{\circ}]$")
    plot.ylabel(r"$J_{"+wave+r"}(\theta)$")
    plot.xscale('log')
    plot.yscale('log')
    plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
    plot.tick_params('y', which='both', direction='in', left=True, right=True)
    leg = plot.legend(frameon=False, markerscale=50, loc=3)
    for line in leg.get_lines():
        line.set_linewidth(6)

    plot.ylim(bottom=.95*low_jvals[-1], top=1.05*up_jvals[0])
    plot.xlim(left=thetas[1], right=thetas[-1])

    p.savefig('./j_factors/plots/dwarf_jfactors_'+wave+'.pdf', bbox_inches='tight')
    p.clf()

for dwarf in dwarflist:
    ymax = 0
    ymin = 1000000000000

    ave_jvals, up_jvals, low_jvals, thetas = [], [], [], []
    p = plot.figure()
    for wave in ['s', 'p', 'd', 'som']:
        with np.load("./j_factors/j_ave/j_ave_thetas_"+dwarf+"_"+wave+".txt", 'rb') as infile:
            ave_jvals = infile['ave_jvals']
            up_jvals = infile['up_jvals']
            low_jvals = infile['low_jvals']
            thetas = infile['thetas']

        thetas = np.array(thetas*180/np.pi, dtype=float)
        ave_jvals = np.array(ave_jvals, dtype=float)
        up_jvals = np.array(up_jvals, dtype=float)
        low_jvals = np.array(low_jvals, dtype=float)

        if up_jvals[0] > ymax:
            ymax = up_jvals[0]
        if low_jvals[-1] < ymin:
            ymin = low_jvals[-1]

        plot.plot(thetas, ave_jvals, '-', color=wavedict[wave], label=wavelabeldict[wave])
        plot.plot(thetas, up_jvals, '-', linewidth=.01, color=wavedict[wave])
        plot.plot(thetas, low_jvals, '-', linewidth=.01, color=wavedict[wave])
        plot.fill_between(thetas, up_jvals, ave_jvals, facecolor=wavedict[wave], alpha=.4)
        plot.fill_between(thetas, ave_jvals, low_jvals, facecolor=wavedict[wave], alpha=.4)

    plot.xlabel(r"$\theta[^{\circ}]$")
    plot.ylabel(r"$J(\theta)$ for " + dwarf)
    plot.xscale('log')
    plot.yscale('log')
    plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
    plot.tick_params('y', which='both', direction='in', left=True, right=True)
    leg = plot.legend(frameon=False, markerscale=50, loc=3)
    for line in leg.get_lines():
        line.set_linewidth(6)

    plot.ylim(bottom=.95*ymin, top=1.05*ymax)
    plot.xlim(left=thetas[1], right=thetas[-1])

    p.savefig('./j_factors/plots/dwarf_jfactors_'+dwarf+'.pdf', bbox_inches='tight')
    p.clf()
