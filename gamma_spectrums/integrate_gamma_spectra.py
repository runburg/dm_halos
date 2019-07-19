#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Integrate the gamma spectra for various channels.

Author: Jack Runburg
Date: 15-07-2019 17:40


"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate


def integrate_spectra_table(data, names):
    """Sort through table and integrate all mass values."""
    result_labels = names.tolist()
    result_labels.pop(1)
    result_labels = np.char.strip(result_labels, '\\[]')
    m_unique, indices = np.unique(data[:, 0], return_index=True)
    # color = iter(plt.cm.rainbow(np.linspace(0, 1, 28)))
    results_dtype = np.dtype({'names': result_labels, 'formats': [np.float]*len(result_labels)})
    results_index = 0
    results = np.zeros(len(m_unique), dtype=results_dtype)
    for i, j, m in zip(indices, np.concatenate((indices[1:], [-1])), m_unique):
        x = np.log10(data[i:j, 1])
        temp_row = [m]
        for column, label in zip(data[i:j, 2:].T, names[2:]):
            integrand = interpolate.interp1d(x, column, kind='cubic', fill_value='extrapolate')
            temp_row.append(integrate.quad(integrand, x[0], x[-1])[0])
        results[results_index] = tuple(temp_row)
        results_index += 1

    return results


def main(file="/Users/runburg/github/dm_halos/gamma_spectrums/gammma_spectra.dat", outfile='/Users/runburg/github/dm_halos/gamma_spectrums/integrated_spectra.npy'):
    """Load and create data table."""
    data = np.loadtxt(file, dtype='str')
    columns = data[0]
    spectra = data[1:].astype(np.float)

    # only keep energies between 1-100 GeV
    spectra[:, 1] = np.power(10, spectra[:, 1])*spectra[:, 0]
    spectra = spectra[[energy >= 1 and energy <= 100 for energy in spectra[:, 1]]]

    # integrate the spectra for each channels
    results = integrate_spectra_table(spectra, columns)
    np.save(outfile, results)
    np.savetxt(outfile[:-4]+'.txt', results, fmt='%.5g')
    return results


if __name__ == '__main__':
    results = main()
