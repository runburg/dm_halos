# dm_halos
For given dimensionful stellar dispersion data, the work flow goes as:
  - Units are removed in *undimensionalize_changeunits.py*
  - The _gfuncs/*.py_ change the dispersion to be a function of r and v
  - Then *h_funcs.py* and *h_values.py* compute the dimensionless J-factors as functions of dimensionless radius/angle
  - *plots.py* is written to make plots with the dimensionless data
  - *integrated_j_factors.py* produces a file with the integrated J-factor as a function of (dimensionless) angle
  - *j_tots.py* uses this file to calculate J-factors of real dwarfs using the mcmc chain data in */j_factors/dwarfs/*
  - Then *tot_j_factor_plotting.py* makes plots of all dwarfs for all annihilation modes it can also make plots to compare this analysis to known results from the data in *compare_to_other_values.py*
