#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make plots to compare to other computed total J-factors.

Author: Jack Runburg
Date: 05-07-2019 11:49


"""

import numpy as np

# ## 1408.0002
# s wave
# 0.5 deg

geringer2015 = np.array([
    ['CarI', 's', 17.87+0.10, 17.87, 17.87-0.09],
    ['DraI', 's', 18.84+0.12, 18.84, 18.84-0.13],
    ['FnxI', 's', 17.83+0.12, 17.83, 17.83-0.06],
    ['LeoI', 's', 17.84+0.20, 17.84, 17.84-0.16],
    ['LeoII', 's', 17.97+0.20, 17.97, 17.97-0.18],
    ['SclI', 's', 18.54+0.06, 18.54, 18.54-0.05],
    ['SxtI', 's', 17.52+0.28, 17.52, 17.52-0.18],
    ['UMiI', 's', 18.93+0.27, 18.93, 18.93-0.19],
    ['BooI', 's', 18.24+0.40, 18.24, 18.24-0.37],
    ['CBerI', 's', 19.02+0.37, 19.02, 19.02-0.41],
    ['CVenI', 's', 17.43+0.37, 17.43, 17.43-0.28],
    ['CVenII', 's', 17.65+0.45, 17.65, 17.65-0.43],
    ['HerI', 's', 16.86+0.74, 16.86, 16.8-0.68],
    ['Seg1', 's', 19.36+0.32, 19.36, 19.36-0.35],
    ['UMaI', 's', 17.87+0.56, 17.87, 17.87-0.33],
    ['UMaII', 's', 19.42+0.44, 19.42, 19.42-0.42]
])
geringer2015[:, [2, 3, 4]] = geringer2015[:, [4, 3, 2]]


# ## 1712.03188
# s and som
# total
bergstrom2017 = np.array([
    ['BooI', 's', 'som', 17.95+0.54, 17.95, 17.95-0.74, 21.13+0.40, 21.13, 21.13-0.48],
    ['UMaII', 's', 'som', 19.87+0.27, 19.87, 19.87-0.18, 22.76+0.25, 22.76, 22.76-0.14],
    ['CVenII', 's', 'som', 18.49+0.31, 18.49, 18.49-0.70, 21.27+0.23, 21.27, 21.27-0.46],
    ['HerI', 's', 'som', 18.12+0.27, 18.12, 18.12-0.35, 21.35+0.22, 21.35, 21.35-0.31],
    ['UMaI', 's', 'som', 18.22+0.95, 18.22, 18.22-0.58, 21.52+0.66, 21.52, 21.52-0.70],
    ['Wil1', 's', 'som', 19.69+0.31, 19.69, 19.69-0.52, 22.54+0.29, 22.54, 22.54-0.23],
    ['CBerI', 's', 'som', 19.42+0.28, 19.42, 19.42-0.45, 22.35+0.21, 22.35, 22.35-0.31],
    ['Seg1', 's', 'som', 19.26+0.48, 19.26, 19.26-0.46, 22.72+0.42, 22.72, 22.72-0.44],
    ['UMiI', 's', 'som', 19.57+0.08, 19.57, 19.57-0.25, 22.62+0.06, 22.62, 22.62-0.27],
    ['CVenI', 's', 'som', 18.01+0.28, 18.01, 18.01-0.29, 21.11+0.29, 21.11, 21.11-0.25],
    ['LeoI', 's', 'som', 17.68+0.23, 17.68, 17.68-0.17, 20.56+0.29, 20.56, 20.56-0.13],
    ['DraI', 's', 'som', 18.78+0.21, 18.78, 18.78-0.26, 21.65+0.23, 21.65, 21.65-0.16],
    ['SxtI', 's', 'som', 18.73+0.22, 18.73, 18.73-0.19, 21.86+0.16, 21.86, 21.86-0.18],
    ['CarI', 's', 'som', 17.71+0.79, 17.71, 17.71-0.02, 20.84+0.86, 20.84, 20.84-0.02],
    ['SclI', 's', 'som', 18.92+0.10, 18.92, 18.92-0.14, 21.94+0.12, 21.94, 21.94-0.15],
    #['SgrII', 's', 'som', 20.25+0.09, 20.25, 20.25-0.12, 23.16+0.09, 23.16, 23.16-0.11],
    ['FnxI', 's', 'som', 18.94+0.08, 18.94, 18.94-0.07, 21.88+0.12, 21.88, 21.88-0.11]
])
bergstrom2017_s = bergstrom2017[:, [0, 1, 5, 4, 3]]
bergstrom2017_som = bergstrom2017[:, [0, 2, 8, 7, 6]]

bergstrom2017_som[:, [2, 3, 4]] = bergstrom2017_som[:, [2, 3, 4]].astype(np.float) + np.log10(100/(2*np.pi))
# bergstrom2017_som = bergstrom2017_som[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16]]

# ## 1802.06811
# s wave
# 0.5 degree
pace2018 = np.array([
    ['CVenI', 's', 17.42+0.17, 17.42, 17.42-0.15],
    ['CarI', 's', 17.83+0.10, 17.83, 17.83-0.09],
    ['DraI', 's', 18.83+0.12, 18.83, 18.83-0.12],
    ['FnxI', 's', 18.09+0.10, 18.09, 18.09-0.10],
    ['LeoI', 's', 17.64+0.14, 17.64, 17.64-0.12],
    ['LeoII', 's', 17.76+0.22, 17.76, 17.76-0.18],
    ['SclI', 's', 18.58+0.05, 18.58, 18.58-0.05],
    ['SxtI', 's', 17.73+0.13, 17.73, 17.73-0.12],
    ['UMiI', 's', 18.75+0.12, 18.75, 18.75-0.12],
    ['AquII', 's', 18.27+0.66, 18.27, 18.27-0.58],
    ['BooI', 's', 18.17+0.31, 18.17, 18.17-0.29],
    ['CVenII', 's', 17.82+0.47, 17.82, 17.82-0.47],
    ['CarII', 's', 18.25+0.55, 18.25, 18.25-0.54],
    ['CBerI', 's', 19.00+0.36, 19.00, 19.00-0.35],
    ['HerI', 's', 17.37+0.53, 17.37, 17.37-0.53],
    # ['Horologium I, des', 's', 18.79+0.90, 18.79, 18.79-0.86],
    ['HorI', 's', 19.27+0.77, 19.27, 19.27-0.71],
    # ['Reticulum II, des', 's', 18.88+0.39, 18.88, 18.88-0.37],
    ['RetII', 's', 18.96+0.38, 18.96, 18.96-0.37],
    ['Seg1', 's', 19.12+0.49, 19.12, 19.12-0.58],
    ['TucII', 's', 18.84+0.55, 18.84, 18.84-0.5],
    ['UMaI', 's', 18.26+0.29, 18.26, 18.26-0.27],
    ['UMaII', 's', 19.44+0.41, 19.44, 19.44-0.39],
    ['Wil1', 's', 19.53+0.50, 19.53, 19.53-0.50]
])
pace2018[:, [2, 4]] = pace2018[:, [4, 2]]

# ## 1804.05052
# s and Sommerfeld
# total
petac2018 = np.array([
    ['CarI', 's', 'som', 17.63+0.16, 17.63, 17.63-0.11, 20.76+0.09, 20.76, 20.76-0.18],
    ['DraI', 's', 'som', 18.67+0.12, 18.67, 18.67-0.15, 21.51+0.15, 21.51, 21.51-0.12],
    ['FnxI', 's', 'som', 17.95+0.37, 17.95, 17.95-0.29, 20.99+0.39, 20.99, 20.99-0.24],
    ['LeoI', 's', 'som', 17.63+0.20, 17.63, 17.63-0.13, 20.57+0.16, 20.57, 20.57-0.10],
    ['LeoII', 's', 'som', 17.65+0.21, 17.65, 17.65-0.31, 20.58+0.30, 20.58, 20.58-0.17],
    ['SclI', 's', 'som', 18.35+0.13, 18.35, 18.35-0.13, 21.37+0.16, 21.37, 21.37-0.13],
    ['SxtI', 's', 'som', 17.30+1.45, 17.30, 17.30-0.36, 20.54+1.37, 20.54, 20.54-0.34],
    ['UMiI', 's', 'som', 18.65+0.14, 18.65, 18.65-0.13, 21.63+0.21, 21.63, 21.63-0.13]
])
petac2018_s = petac2018[:, [0, 1, 5, 4, 3]]
petac2018_som = petac2018[:, [0, 2, 8, 7, 6]]
petac2018_som[:, [2, 3, 4]] = petac2018_som[:, [2, 3, 4]].astype(np.float) + np.log10(100/np.pi)

# ## 1702.00408
# Sommerfeld
# 0.5
kumar2017_som = np.array([
    ['RetII', 'som', np.log10(2.1*10**21*100/2/np.pi), np.log10(4.8*10**21*100/2/np.pi), np.log10(9.5*10**21*100/2/np.pi)],
    ['CBerI', 'som', np.log10(2*10**21*100/2/np.pi), np.log10(4*10**21*100/2/np.pi), np.log10(7.2*10**21*100/2/np.pi)],
    ['Seg1', 'som', np.log10(4.1*10**21*100/2/np.pi), np.log10(1.8*10**22*100/2/np.pi), np.log10(4*10**22*100/2/np.pi)],
    ['DraI', 'som', np.log10(1.8*10**21*100/2/np.pi), np.log10(3.3*10**21*100/2/np.pi), np.log10(6*10**21*100/2/np.pi)],
    ['UMiI', 'som', np.log10(2.2*10**21*100/2/np.pi), np.log10(5*10**21*100/2/np.pi), np.log10(9*10**21*100/2/np.pi)]
])

# ## 1511.06296
# s
# 0.5
walker2015 = np.array(['TucII', 's', 18.7-0.7, 18.7, 18.7+0.9])

# ## 1612.06398
# s
# 0.5
caldwell2016 = np.array(['CraII', 's', 15.7-0.25, 15.7, 15.7+0.25])

# ## 1804.06430
# s
# 0.5
koposov2018 = np.array(['HyiI', 's', 18.33-0.34, 18.33, 18.33+0.38])

# ## 1601.02181
# s, p, som
vels = np.array([9.1, 4.3, 4.6, 4.3, 6.7, 9.5])
c = 2.99e+05
zhao2016_s = np.array([
    ['DraI', 's', 18.8 - 0.16, 18.8, 18.8 + 0.16],
    ['Seg1', 's', 19.5 - 0.29, 19.5 , 19.5 + 0.29],
    ['CBerI', 's', 19.0 - 0.25, 19.0, 19.0 + 0.25],
    ['Wil1', 's', 19.1 - 0.31, 19.1, 19.1 + 0.31],
    ['UMaII', 's', 19.3 - 0.28, 19.3, 19.3 + 0.28],
    ['UMiI', 's', 18.8 - .19, 18.8 , 18.8 + .19]
])

zhao2016_p = np.copy(zhao2016_s)
zhao2016_p[:, 1] = ['p'] * len(zhao2016_p[:, 1])
for i in range(2,5):
    zhao2016_p[:, i] = zhao2016_p[:, i].astype(np.float) + np.log10(6 / c**2 * vels**2)

zhao2016_som = np.copy(zhao2016_s)
zhao2016_som[:, 1] = ['som'] * len(zhao2016_som[:, 1])
for i in range(2,5):
    zhao2016_som[:, i] = zhao2016_som[:, i].astype(np.float) + np.log10(np.sqrt(np.pi) * c / vels)

# ## 1711.04696
# p
zhao2017 = np.array([
    ['Wil1', 'p', (19.5-0.6), 19.5, (19.5+1.3)],
    ['RetII', 'p', (19.5-0.8), 19.5, (19.5+1.1)]
])

zhao2017[:, 2:] = zhao2017[:, 2:].astype(np.float)-np.array(np.log10(c**2))


np.savez('./total_j_factors/j_factors/compare_values_j_factors.npz', geringer2015=geringer2015, bergstrom2017_s=bergstrom2017_s, bergstrom2017_som=bergstrom2017_som, pace2018=pace2018, petac2018_s=petac2018_s, petac2018_som=petac2018_som, kumar2017_som=kumar2017_som, walker2015=walker2015, caldwell2016=caldwell2016, koposov2018=koposov2018, zhao2016_s=zhao2016_s, zhao2016_p=zhao2016_p, zhao2016_som=zhao2016_som, zhao2017=zhao2017)
