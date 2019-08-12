# dm_halos
For given dimensionful stellar dispersion data, the work flow goes as:
  - Units are removed in *undimensionalize_changeunits.py*
  - The _gfuncs/*.py_ change the dispersion to be a function of r and v
  - Then *h_funcs.py* and *h_values.py* compute the dimensionless J-factors as functions of dimensionless radius/angle
    - *main.py* can complete the above steps
  - *plots.py* is written to make plots with the dimensionless data
To compute dimensionful J-factors,
  - *integrated_j_factors.py* produces a file with the integrated (dimensionless) J-factor as a function of (dimensionless) angle
  - *j_tots.py* uses this file to calculate J-factors of real dwarfs using the mcmc chain data in */j_factors/dwarfs/*
  - Then *tot_j_factor_plotting.py* makes plots of all dwarfs for all annihilation modes it can also make plots to compare this analysis to known results from the data in *compare_to_other_values.py*
To get DM annihilation bounds,
  - Using data from PPPC4DMID, *integrate_gamma_spectra.py* provides the integrated spectra for various _channels/_
  - *madhat_formatting.py* makes the files with J-factors for MADHAT in _models/_
  - *batch_madhat.sh* computes bounds for all _models/*.txt_ and _channels/*.txt_
  - *madhat_plotting.py* produces the cross-section bound plots
