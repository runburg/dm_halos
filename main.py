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
file = "df_nfw"
# units of kpc/solar mass * (km/s)^2
g_n = float(4.325E-6)
# units of kpc
r_s = float(20.0)
# units of solar mass/kpc^3
rho_s = float(8E6)

undim.undim(file, g_n, r_s, rho_s)
print("firnished undim")
gs.gs_wave(file)
print("finished s wave")
gp.gp_wave(file)
print("finsished p wave")
gd.gd_wave(file)
print("finsihed d wave")
gsom.gsom_wave(file)
print("finished som enh")
hval.hval()

with open("_parameters.txt") as outfile:
    outfile.write("r_s: " + r_s)
    outfile.write("rho_s: " + rho_s)
    outfile.write("input file is: " + file)
