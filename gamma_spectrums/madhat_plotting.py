#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot spectra from MADHAT runs.

Author: Jack Runburg
Date: 22-07-2019 11:37


"""
import numpy as np
import matplotlib.pyplot as plt
import glob


def import_madhat_output(infile):
    """Created structured arrays from madhat data files."""
    with open(infile, 'r') as input:
        # read in and format columns for array
        first_line = input.readline().split('  ')
        first_line = list(map(str.strip, list(filter(lambda x: x != '', first_line))))
        first_line = [line.partition('(')[0] for line in first_line]

        #create structured array of madhat data
        array_from_file = []

        # read in lines and create list of floats
        for line in input.readlines():
            array_from_file.append(tuple(np.array(line.split()).astype(np.float)))

        # define the names and formats of the strcutred array and initialize array
        dt = np.dtype({'names': first_line, 'formats': ['f8']*len(first_line)})
        array_from_file = np.array(array_from_file, dtype=dt)

    return array_from_file


def style_ax(ax):
    """Style subplot axes."""

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(left=7.5, right=5500)
    ax.set_xlabel(r'$m_\chi$ [GeV]')
    ax.set_ylabel(r'$(\sigma v)_0$ [$\mathrm{cm}^3\mathrm{s}^{-1}$]')


def cross_section_plots(files):
    """Plot the s and som cutoffs of cross sections."""
    # get masses for plotting
    masses = np.load('./gamma_spectrums/integrated_spectra.npy')['mDM']
    # set up figure
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    figsize = [12, 5]
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=figsize)

    # define colors for the different channels
    color_dict = dict(zip(['Tau', 'b', 'Mu', 'W'],['b', 'r', 'g', 'k']))
    for filename in files:
        # get channel and wave names
        channel = filename.partition('_')[2].partition('_')
        wave = channel[2].strip('.out')
        channel = channel[0]

        # plot s and som on different subplots
        if wave == 's':
            ax = axs[0]
        if wave == 'som':
            ax = axs[1]

        # get cross section values
        cs = import_madhat_output(filename)['cs']
        ax.plot(masses[-len(cs):], cs, label=channel, color = color_dict[channel], lw=5)
        ax.legend(loc='lower right')

    for ax in axs:
        style_ax(ax)
    fig.savefig('./gamma_spectrums/cross_section_plots.pdf', bbox_inches='tight')

def main():
    """Create plots using madhat data."""
    path_to_files = '/Users/runburg/Desktop/MADHAT/Output/channel_*.out'
    cross_section_plots(glob.glob(path_to_files))


if __name__ == '__main__':
    main()
