#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Jack Runburg
# Date: 15-03-2019 11:17
# main.py for running dark matter halo distribution calculation
#
import undimensionalize_changeunits as undim
import gfuncs.g_swave as gs
import gfuncs.g_pwave as gp
import gfuncs.g_dwave as gd
import gfuncs.g_som as gsom
import h_values as hval

# file without extension
file = "fe_GC_NFW"
# units of kpc/solar mass * (km/s)^2
g_n = float(4.325E-6)
# units of kpc
r_s = float(20.0)
# units of solar mass/kpc^3
rho_s = float(8E6)

undim.undim(file, g_n, r_s, rho_s)
gs.gs_wave(file)
gp.gp_wave(file)
gd.gd_wave(file)
gsom.gsom_wave(file)
hval.hval()

print(file)
