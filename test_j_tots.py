#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test to see if the j_tots are smooth.

Author: Jack Runburg
Date: 28-06-2019 19:35


"""

import numpy as np
import matplotlib.pyplot as plt

with np.load('./j_factors/j_thetas.txt', 'rb') as infile:
    js = infile['js']
    jp = infile['jp']
    jd = infile['jd']
    jsom = infile['jsom']
    tildetheta = infile['theta']

with np.load('./j_factors/int_dimless_j_fac.npy', 'rb') as infile:
    tildetheta = infile['angle']
    js = infile['js']
    jp = infile['jp']
    jd = infile['jd']
    jsom = infile['jsom']

fig, axs = plt.subplots()

for d, c, l in zip([js, jp, jd, jsom], ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', 'xkcd:light turquoise'], [r"$s$-wave, ", r"$p$-wave, ", r"$d$-wave, ", r"Sommerfeld, "]):
    axs.plot(tildetheta, d, color=c, label=l)
print(js)
plt.legend()
plt.show()
