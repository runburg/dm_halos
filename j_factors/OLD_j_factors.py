#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate dimensionful J-factors.

Author: Jack Runburg
Date: 20-05-2019 10:11


"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot

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

scolor = 'xkcd:azure'
pcolor = 'xkcd:coral'
dcolor = 'xkcd:peach'
somcolor = 'xkcd:light turquoise'

slabel = r"$s$-wave, $n=0$"
plabel = r"$p$-wave, $n=2$"
dlabel = r"$d$-wave, $n=4$"
somlabel = r"Sommerfeld, $n=-1$"

with np.load('./j_factors/j_thetas.txt', 'rb') as infile:
    jst = infile['jst']
    js = infile['js']
    jp = infile['jp']
    jd = infile['jd']
    jsom = infile['jsom']
    theta = infile['theta']

data = np.load('./j_factors/draco1_chain.npy')

jvals = []
rs_ave, rhos_ave, d_ave, weightsum = 0, 0, 0, 0
for x in data:
    rs_ave += x[1]*x[0]
    rhos_ave += x[2]*x[0]
    d_ave += x[3]*x[0]
    weightsum += x[0]

rs_ave = rs_ave/weightsum
rhos_ave = rhos_ave/weightsum
d_ave = d_ave/weightsum

rsigma, rhosigma, dsigma = 0, 0, 0
for x in data:
    rsigma = x[0]*(x[1]-rs_ave)**2
    rhosigma = x[0]*(x[2]-rhos_ave)**2
    dsigma = x[0]*(x[3]-d_ave)**2

rsigma = np.sqrt(rsigma/((len(data)-1)*weightsum/len(data)))
rhosigma = np.sqrt(rhosigma/((len(data)-1)*weightsum/len(data)))
dsigma = np.sqrt(dsigma/((len(data)-1)*weightsum/len(data)))

factor = 2 * rs_ave * rhos_ave**2
# np.sqrt(4 * np.pi * float(4.325E-6) * x[2] * x[1]**2 / (299792458)**2)

z = 1.96/np.sqrt(len(data))
print('ave: ' + str(factor))
print('ur: ' + str(2 * (rs_ave+z*rsigma) * (rhos_ave+z*rhosigma)**2))
print('br: ' + str(2 * (rs_ave+z*rsigma) * (rhos_ave-z*rhosigma)**2))
print('ul: ' + str(2 * (rs_ave-z*rsigma) * (rhos_ave+z*rhosigma)**2))
print('bl: ' + str(2 * (rs_ave-z*rsigma) * (rhos_ave-z*rhosigma)**2))

theta = rs_ave/d_ave*theta

p = plot.figure()

swave = plot.plot(theta, js, '-', color=scolor, label=slabel)
som = plot.plot(theta, jst, '-', color=somcolor, label=somlabel)

plot.xlabel(r"$\theta$")
plot.ylabel(r"$J_s(\theta)$")
plot.xscale('log')
plot.yscale('log')
plot.tick_params('x', which='both', direction='in', bottom=True, top=True)
plot.tick_params('y', which='both', direction='in', left=True, right=True)
leg = plot.legend(frameon=False, markerscale=50, loc=3)
for line in leg.get_lines():
    line.set_linewidth(6)


p.savefig('./j_factors/dwarf_jfactors.pdf', bbox_inches='tight')
