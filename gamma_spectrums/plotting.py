#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the spectra of the integrated gamma spectra.

Author: Jack Runburg
Date: 18-07-2019 09:55


"""
import numpy as np
import matplotlib.pyplot as plt


def main():
    """Plot the spectra for each decay."""
    data = np.load('/Users/runburg/github/dm_halos/gamma_spectrums/integrated_spectra.npy')

    # \phi_pp values from 1802.03826
    phi_pp = [3.5E-30, 3.86E-30, 1.62E-30, 2.76E-30, 7.41E-33]
    masses = data['mDM']
    particle_physics_scenario_1_factor = np.array(8 * np.pi * phi_pp[0] * masses**2)
    print(particle_physics_scenario_1_factor)
    channels = ['Tau', 'b', 'Mu', 'W']
    y_vals = data[channels]
    for channel, color in zip(channels, ['b', 'r', 'g', 'k']):
        y_vals[channel] = particle_physics_scenario_1_factor / y_vals[channel]
        plt.plot(masses, y_vals[channel], color=color, label=channel)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig('./gamma_spectrums/integrated_spectra_plot.pdf')
    plt.show()

if __name__ == '__main__':
    main()
