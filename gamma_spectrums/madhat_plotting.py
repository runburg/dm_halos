#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot spectra from MADHAT runs.

Author: Jack Runburg
Date: 22-07-2019 11:37


"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
sys.path.insert(0, './total_j_factors/')
from tot_j_factor_plotting import lighten_color


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


def style_ax(ax, axs):
    """Style subplot axes."""

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(left=7.5, right=5500)
    ax.set_xlabel(r'$m_\chi$ [GeV]')
    ax.set_ylabel(r'$(\sigma v)_0$ [$\mathrm{cm}^3\mathrm{s}^{-1}$]')
    if ax is not axs:
        if axs.shape == (2,):
            if ax == axs[0]:
                ax.set_ylim(top=10**(-21))
            if ax == axs[1]:
                ax.set_ylim(top=10**(-23))
        else:
            axs[0,0].set_ylim(bottom=10**(-28), top=10**(-21))
            axs[0,1].set_ylim(bottom=10**(-21), top=10**(-14))
            axs[1,0].set_ylim(bottom=10**(-14), top=10**(-7))
            axs[1,1].set_ylim(bottom=10**(-32), top=10**(-25))


def smooth_plot(x, y):
    from scipy.interpolate import splrep, splev

    x_new = np.logspace(np.log10(x.min()), np.log10(x.max()), num=200)
    spline = splrep(np.log10(x),np.log10(y),s=5)
    spliney = splev(np.log10(x), spline)

    return x, 10**spliney


def cross_section_channel_plots(files, outfile):
    """Plot the s and som cutoffs of cross sections."""
    # get masses for plotting
    masses = np.load('./gamma_spectrums/integrated_spectra.npy')['mDM']
    # set up figure
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    lw=2

    figsize = [12, 5]
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=figsize)

    # define colors for the different channels
    color_dict = dict(zip(['Tau', 'b', 'Mu', 'W'],['b', 'r', 'g', 'k']))
    for filename in files:
        # get channel and wave names
        channel = filename.partition('_')[2].partition('_')
        wave = channel[2].partition('_')[2].strip('.txt')
        channel = channel[0]
        linestyle='-'

        # plot s and som on different subplots
        if wave == 's':
            ax = axs[0]
        if wave == 'som':
            ax = axs[1]

        # get cross section values and check for finiteness
        data =  import_madhat_output(filename)
        indices = np.isfinite(data['cs'])
        xmasses = masses[indices]
        data = data[indices]
        cs = data['cs']
        dcs_up = data['+dcs']
        dcs_down = data['-dcs']

        # pick color of plots
        color = color_dict[channel]

        # plot s and som on different subplots
        if wave == 's':
            ax = axs[0]
        if wave == 'som':
            ax = axs[1]
        label = channel

        central = smooth_plot(xmasses, cs)
        up = smooth_plot(xmasses, dcs_up)
        down = smooth_plot(xmasses, dcs_down)
        ax.plot(*central, label=label, color=color, lw=lw, linestyle=linestyle)
        ax.plot(*up, color=lighten_color(color, amount=0.7), lw=lw*.2, linestyle=linestyle)
        ax.plot(*down, color=lighten_color(color, amount=0.7), lw=lw*.2, linestyle=linestyle)
        ax.fill_between(*central, up[-1], color=lighten_color(color, amount=0.3))
        ax.fill_between(*central, down[-1], color=lighten_color(color, amount=0.3))

    for ax in axs:
        style_ax(ax, axs)
        x = np.logspace(np.log10(10), np.log10(5000), num=30)
        ax.plot(x, [3*10**(-26)]*len(x), lw=lw, linestyle='--', color='gray')

    fig.savefig(f'./gamma_spectrums/madhat_plots/{outfile}', bbox_inches='tight')


def cross_section_wave_plots(files, outfile, nrows=1, ncols=1):
    """Plot the s and som cutoffs of cross sections."""
    # get masses for plotting
    masses = np.load('./gamma_spectrums/integrated_spectra.npy')['mDM']

    # set up figure
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    lw=2

    figsize = [2+ncols*5, nrows*5]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize=figsize)

    ax = axs
    colors = iter(['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), lighten_color('xkcd:purple', amount=0.8), 'xkcd:hot pink'])
    for filename in files:
        # get channel and wave name
        wave = filename.partition('_')[2].partition('_')
        channel = wave[0]
        wave = wave[2].partition('_')[2].strip('.txt')

        # get cross section values and check for finiteness
        data =  import_madhat_output(filename)
        indices = np.isfinite(data['cs'])
        xmasses = masses[indices]
        data = data[indices]
        cs = data['cs']
        dcs_up = data['+dcs']
        dcs_down = data['-dcs']

        # pick color of plots
        color = next(colors)

        # labels
        labels = {'s': r"$s$-wave", 'p': r"$p$-wave", 'd': r"$d$-wave", 'som': r"Som. enh."}
        label = labels[wave]

        # ax.errorbar(*smooth_plot(xmasses, cs), yerr=[dcs_down, dcs_up], label=label, color=color, lw=lw)
        central = smooth_plot(xmasses, cs)
        up = smooth_plot(xmasses, dcs_up)
        down = smooth_plot(xmasses, dcs_down)
        ax.plot(*central, label=label, color=color, lw=lw)
        ax.plot(*up, color=lighten_color(color, amount=0.7), lw=lw*.2)
        ax.plot(*down, color=lighten_color(color, amount=0.7), lw=lw*.2)
        ax.fill_between(*central, up[-1], color=lighten_color(color, amount=0.3))
        ax.fill_between(*central, down[-1], color=lighten_color(color, amount=0.3))
        ax.legend(loc='lower right')


    style_ax(ax, axs)
    x = np.logspace(np.log10(1), np.log10(10000), num=30)
    ax.plot(x, [3*10**(-26)]*len(x), lw=lw, linestyle='--', color='gray')

    fig.savefig(f'./gamma_spectrums/madhat_plots/{outfile}', bbox_inches='tight')


def cross_section_all_waves_plots(files, outfile, nrows=2, ncols=2):
    """Plot the s and som cutoffs of cross sections."""
    # get masses for plotting
    masses = np.load('./gamma_spectrums/integrated_spectra.npy')['mDM']

    # set up figure
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    lw=2

    figsize = [2+ncols*5, nrows*5]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, figsize=figsize)

    ax = axs
    colors = dict(zip(['Tau', 'Mu', 'b',  'W'],['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), lighten_color('xkcd:purple', amount=0.8), 'xkcd:hot pink']))
    for filename in files:
        # get channel and wave name
        wave = filename.partition('_')[2].partition('_')
        channel = wave[0]
        wave = wave[2].partition('_')[2].strip('.txt')

        linestyle = '-'
        alpha=0.3
        label = channel

        if wave.partition('_')[2] is not '':
            wave = wave.partition('_')[0]
            label = None
            linestyle = ':'

        # labels & titles

        titles = {'s': r"$s$-wave", 'p': r"$p$-wave", 'd': r"$d$-wave", 'som': r"Som. enh."}
        title = titles[wave]

        # get cross section values and check for finiteness
        data =  import_madhat_output(filename)
        indices = np.isfinite(data['cs'])
        # xmasses = masses[indices]
        data = data[indices]
        masses = data['Mass']
        cs = data['cs']
        dcs_up = cs + data['+dcs']
        dcs_down = cs - data['-dcs']

        # choose subplot
        if wave == 's':
            ax = axs[0, 0]
        if wave == 'p':
            ax = axs[0, 1]
        if wave == 'd':
            ax = axs[1, 0]
        if wave == 'som':
            ax = axs[1, 1]

        # pick color of plots
        color = colors[channel]

        # ax.errorbar(*smooth_plot(xmasses, cs), yerr=[dcs_down, dcs_up], label=label, color=color, lw=lw)
        central = smooth_plot(masses, cs)
        up = smooth_plot(masses, dcs_up)
        down = smooth_plot(masses, dcs_down)
        ax.plot(*central, label=label, color=color, lw=lw, linestyle=linestyle)
        ax.plot(*up, color=lighten_color(color, amount=0.7), lw=lw*.2, linestyle=linestyle)
        ax.plot(*down, color=lighten_color(color, amount=0.7), lw=lw*.2, linestyle=linestyle)
        ax.fill_between(*central, up[-1], color=lighten_color(color, amount=0.3), alpha=alpha)
        ax.fill_between(*central, down[-1], color=lighten_color(color, amount=0.3), alpha=alpha)
        ax.legend(loc='lower right')
        ax.set_title(title)

    for ax in axs.flatten():
        style_ax(ax, axs)
    # x = np.logspace(np.log10(1), np.log10(10000), num=30)
    # ax.plot(x, [3*10**(-26)]*len(x), lw=lw, linestyle='--', color='gray')

    fig.savefig(f'./gamma_spectrums/madhat_plots/{outfile}', bbox_inches='tight')


def main():
    """Create plots using madhat data."""
    path_to_files = './madhat/jackoutput/refpaperoutput/channel_*.txt'
    cross_section_channel_plots(glob.glob(path_to_files), 'fig5.pdf')

    # path_to_files = '/Users/runburg/Desktop/MADHAT/jackoutput/channel_Tau_*.txt'
    # cross_section_wave_plots(glob.glob(path_to_files), 'tau_jack_data.pdf')

    path_to_files = './madhat/jackoutput/channel_*.txt'
    cross_section_all_waves_plots(glob.glob(path_to_files), 'jack_data.pdf')


if __name__ == '__main__':
    main()
