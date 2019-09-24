#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make bar and whisker plots of all total j-factors.

Author: Jack Runburg
Date: 26-06-2019 17:08

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def configure_plot_params(fontsize=20, spacing=0.2):
    """Make plots look the way I want."""
    mpl.rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'weight': 'light'})
    mpl.rc('text', usetex=True)
    mpl.rc('text.latex', preamble=r'\usepackage{amsmath, graphicx}')
    mpl.rcParams['mathtext.default'] = 'regular'
    # mpl.rcParams['font.serif'] = 'Computer Modern Roman'
    # plt.rcParams["font.family"] = ""
    # mpl.rcParams['font.size'] = fontsize
    # mpl.rcParams['mathtext.fontset'] = 'cm'
    # mpl.rcParams['axes.labelsize'] = 'xx-large'
    # mpl.rcParams['axes.titlesize'] = 'xx-large'

    mpl.rcParams['figure.subplot.wspace'] = spacing
    mpl.rcParams['figure.subplot.hspace'] = spacing

    # mpl.rcParams['xtick.labelsize'] = 'xx-large'
    mpl.rcParams['xtick.major.size'] = 7.5
    mpl.rcParams['xtick.major.width'] = 1
    mpl.rcParams['xtick.minor.size'] = 3.75
    mpl.rcParams['xtick.minor.width'] = 0.5

    # mpl.rcParams['ytick.labelsize'] = 'xx-large'
    mpl.rcParams['ytick.major.size'] = 7.5
    mpl.rcParams['ytick.major.width'] = 1
    mpl.rcParams['ytick.minor.size'] = 3.75
    mpl.rcParams['ytick.minor.width'] = 0.5


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.

    Input can be matplotlib color string, hex string, or RGB tuple. Make amount > 1 to darken.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)

    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except KeyError:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def sorted_j_factors(angle_jfacs):
    """Return the names of the sorted j-factors."""
    jfacs = np.load(angle_jfacs)
    jfacs[:, 0] = fix_names(jfacs[:, 0])

    js = jfacs[np.arange(0, len(jfacs), 4)]
    jp = jfacs[np.arange(1, len(jfacs), 4)]
    jd = jfacs[np.arange(2, len(jfacs), 4)]
    jsom = jfacs[np.arange(3, len(jfacs), 4)]

    sorted_names = np.empty((len(js)+1, 4), dtype=object)
    sorted_names[0] = np.array(['s', 'p', 'd', 'som'])
    for i, j in enumerate([js, jp, jd, jsom]):
        sorted_names[1:, i] = np.flip(j[j[:, 3].astype(float).argsort()][:, 0])
    np.savetxt(angle_jfacs[:-4] + '_sorted_names.txt', sorted_names, fmt='%s', delimiter=' & ', newline=' \\\\\hline\n')


def error_plot(axs, data, offset=0, fmt='D', color='r', label=None, ec=0, mec=None, ms=6, lw=2):
    """Generate a scatterplot with error bars for the given j-factor data."""
    try:
        ave = data[:, 3].astype(np.float)
        low = ave - data[:, 2].astype(np.float)
        high = data[:, 4].astype(np.float) - ave
    except IndexError:
        ave = data[3].astype(np.float)
        low = [ave - data[2].astype(np.float)]
        high = [data[4].astype(np.float) - ave]

    try:
        x = data[:, 0].astype(np.float) + offset
    except ValueError:
        print('No offset', label)
        print(data[:, 0])
        x = data[:, 0]
    except IndexError:
        x = data[0].astype(np.float) + offset

    axs.errorbar(x, ave, yerr=[low, high], fmt=fmt, ecolor=ec, mec=mec, color=color, mew=lw, capsize=lw+5, label=label, ms=ms, lw=lw)


def fix_names(name_list, get_dict=False):
    """Fix names for plotting of dwarfs."""
    name_dict = {
                'ursamajor2': 'UMaII',
                'willman1': 'Wil1',
                'horologium1_k15': 'HorI',
                'segue1': 'Seg1',
                'comaberenices': 'CBerI',
                'tucana2_des': 'TucII',
                'reticulum2_k15': 'RetII',
                'draco1': 'DraI',
                'hydrus1': 'HyiI',
                'ursaminor': 'UMiI',
                'sculptor': 'SclI',
                'carina2': 'CarII',
                'aquarius2': 'AquII',
                'bootes1': 'BooI',
                'ursamajor1': 'UMaI',
                'fornax': 'FnxI',
                'canesvenatici2': 'CVenII',
                'carina': 'CarI',
                'sextans': 'SxtI',
                'leo2': 'LeoII',
                'leo1': 'LeoI',
                'canesvenatici1': 'CVenI',
                'sagittarius2': 'SgrII',
                'hercules': 'HerI',
                'crater2': 'CraII'}

    if get_dict:
        return name_dict

    names = [name_dict[item] for item in name_list]
    return names


def set_x_value_names(data, get_dict=False):
    """For a given dwarf name, return the coordinate it corresponds to."""
    coord_dict = {item: num for (num, item) in enumerate(sorted(fix_names(None, True).values()), start=1)}
    if get_dict:
        return coord_dict

    return [coord_dict[item] for item in data]


def plot_all_dwarfs(infile):
    """Generate the plots for all the dwarfs in the given file."""
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5)]
    labels = [r"$s$-wave", r"$p$-wave", r"$d$-wave", r"Sommerfeld"]
    error_colors = [lighten_color('xkcd:black', 0.4), lighten_color('xkcd:black', 0.6), lighten_color('xkcd:black', 0.8), 'xkcd:black']
    angle_labels = [r"$0.1^\circ$", r"$0.2^\circ$", r"$0.5^\circ$", r"$10^\circ$"]
    fmts = ['*', 'o', 'h', 'D']

    if len(infile) < 4:
        error_colors = ['xkcd:slate gray', 'xkcd:black']
        angle_labels = [r"$0.5^\circ$", r"$10^\circ$"]
        fmts = ['h', 'D']

    configure_plot_params()

    fig, ax = plt.subplots(figsize=(25, 25))

    start = len(fmts)//2-len(fmts)+len(fmts) % 2
    for i, (errorcolor, angle_jfacs, a, fmt) in enumerate(zip(error_colors, infile, angle_labels, fmts), start=start):
        # for errorcolor, angle_jfacs, a in zip(['xkcd:black'], infile, [r"$10^\circ$"]):
        jfacs = np.load(angle_jfacs, allow_pickle=True)
        jfacs[:, 0] = set_x_value_names(fix_names(jfacs[:, 0]))

        js = jfacs[np.arange(0, len(jfacs), 4)]
        jp = jfacs[np.arange(1, len(jfacs), 4)]
        jd = jfacs[np.arange(2, len(jfacs), 4)]
        jsom = jfacs[np.arange(3, len(jfacs), 4)]

        for d, c, l in zip([js, jp, jd, jsom], colors, labels):
            error_plot(ax, d, offset=(i+0.5)/(1.5*len(fmts)), color=c, fmt=fmt, ec=errorcolor, mec=errorcolor, ms=12, lw=2.7)

    # add grey stripes
    for idwarf in np.arange(0, len(js)+2, 2):
        ax.axvline(idwarf, color='xkcd:light gray', lw=56, zorder=-100)

    # set ticks and labels
    ax.tick_params(axis='both', labelsize=26)
    x_labels = set_x_value_names(None, get_dict=True)
    ax.set_xticks(list(x_labels.values()))
    ax.set_xticklabels(list(x_labels.keys()), rotation=90)
    ax.set_xlim(left=0.5, right=25.5)
    ax.xaxis.set_tick_params(which='minor', bottom=False)

    ax.set_ylabel(r"$\log_{10}\,\,J/(\mathrm{GeV}^2\mathrm{cm}^{-5})$", size=28)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_tick_params(which='minor', left=True, right=True)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    # make custom legend
    legend_elements = [Line2D([0], [0], color=c, label=l, lw=12) for c, l in zip(colors, labels)]
    lw = 2.7
    ms = 12
    for ec, l, f in zip(error_colors, angle_labels, fmts):
        # make empty, uncolored legend marker for error plot
        line = Line2D([], [], ls='none', marker=f, mec=ec, color='xkcd:white', label=l, ms=ms, mew=lw, lw=lw)
        barline = LineCollection(np.empty((2, 2, 2)), label=l)
        ebar = ErrorbarContainer((line, [line], [barline]), has_yerr=True)
        ebar.set_label(l)
        legend_elements.append(ebar)

    # reorder legend to beginning
    legend_elements.insert(0, legend_elements.pop(3))
    ax.legend(handles=legend_elements, labels=labels.append(angle_labels), markerscale=1.25, prop={'size': 24}, frameon=True, loc='lower right')
    outfile = './total_j_factors/j_factors/tot_j_factor.pdf'
    fig.savefig(outfile, bbox_inches='tight')


def compare_plot(ax, data, colors, labels, fmt='D', offset=None):
    """Generate the plot for the comparison to other results."""
    configure_plot_params(fontsize=16)

    if len(fmt) != len(data):
        fmt = [fmt[0]]*len(data)

    start = len(data)//2-len(data)+len(data) % 2

    for ioffset, (d, c, l, fm) in enumerate(zip(data, colors, labels, fmt), start=start):
        try:
            d[:, 0] = set_x_value_names(d[:, 0])
        except IndexError:
            d[0] = set_x_value_names([d[0]])[0]
        except KeyError:
            print('Error with key: ', d[:, 0])
        error_plot(ax, d, offset=(ioffset+0.5)/(1.5*len(labels)), color=c, label=l, ec=lighten_color(c, amount=1.5), mec=lighten_color(c, amount=1.5), fmt=fm)

    x_labels = set_x_value_names(None, get_dict=True)
    # padding = len(max(list(x_labels), key=len))

    ax.tick_params(axis='both', labelsize=24)
    ax.set_xlim(left=0.5, right=25.5)
    ax.set_xticks(list(x_labels.values()))
    ax.set_xticklabels(list(x_labels.keys()), rotation=90)
    ax.xaxis.set_tick_params(which='minor', bottom=False)
    ax.xaxis.labelpad = 0

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_tick_params(which='minor', left=True, right=True)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    for idwarf in np.arange(0, len(data[0])+2, 2):
        ax.axvline(idwarf, color='xkcd:light grey', lw=44, zorder=-100)

    ax.set_ylabel(r"$\log_{10}\,\,\,\,{J}/(\mathrm{GeV}^2\mathrm{cm}^{-5})$", size=26)
    ax.legend(markerscale=1.25, prop={'size': 20}, frameon=True, loc='lower right')


def compare_to_others_plot():
    """Generate plots to compare to other papers."""
    mpl.use('pgf')
    fig, ax = plt.subplots(figsize=(20, 20))

    my_data = "./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_half.npy"
    jfacs = np.load(my_data, allow_pickle=True)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]
    jp = jfacs[np.arange(1, len(jfacs), 4)]

    others_data = np.load('./total_j_factors/j_factors/compare_values_j_factors.npz')

    label_dict = {'1408.0002': r'Geringer-Sameth et al. 2015 \cite{Geringer-Sameth2015ApJ...801...74G}', '1802.06811': r'Pace et al. 2019 \cite{Pace:2018tin}', '1511.06296': r'Walker et al. 2016 \cite{Walker:2016adk}', '1612.06398': r'Caldwell et al. 2016 \cite{Caldwell:2016hrl}', '1804.06430': r'Koposov et al. 2018 \cite{Koposov2018MNRAS.479.5343K}', '1712.03188': r'Bergstrom et al. 2018 \cite{Bergstrom2018}', '1804.05052': r'Petac et al. 2018 \cite{Petac:2018gue}', '1702.00408': r'Boddy et al. 2017 \cite{Boddy:2017vpe}', '1601.02181': r'Zhao et al. 2016 \cite{Zhao2016}', '1711.04696': r'Zhao et al. 2018 \cite{Zhao2017}'}

    # s wave, to 0.5 degree
    data = [js, others_data['geringer2015'], others_data['pace2018'], others_data['walker2015'], others_data['caldwell2016'], others_data['koposov2018']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), 'xkcd:purple', 'xkcd:hot pink']
    labels = ['This analysis', label_dict['1408.0002'], label_dict['1802.06811'], label_dict['1511.06296'], label_dict['1612.06398'], label_dict['1804.06430']]

    compare_plot(ax, data, colors, labels, fmt='h')
    plt.title(r"Comparison of {$J$}-factors at \mbox{$0.5^\circ$} for \mbox{$s$}-wave annihilation", size=28)
    outfile = './total_j_factors/j_factors/comparison_plot_s_half.pgf'
    fig.savefig(outfile, bbox_inches='tight')

    # # p wave, half
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))

    data = [jp, others_data['zhao2017']]
    colors = [lighten_color('xkcd:azure', amount=0.5), lighten_color('xkcd:peach', amount=1.5)]
    labels = [r'This analysis, $0.5^\circ$', label_dict['1711.04696'] + r', $0.5^\circ$']

    # compare_plot(ax, data, colors, labels, fmt='h', offset=4)

    # # p wave, total
    my_data = "./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_10.npy"
    jfacs = np.load(my_data)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]
    jsom = jfacs[np.arange(3, len(jfacs), 4)]
    jp = jfacs[np.arange(1, len(jfacs), 4)]

    data.extend([jp, others_data['zhao2016_p']])
    colors.extend(['xkcd:azure', 'xkcd:coral'])
    labels.extend([r'This analysis, $10^\circ$', label_dict['1601.02181'] + r', total'])

    compare_plot(ax, data, colors, labels, fmt=['h', 'h', 'D', 'D'])
    ax.set_title(r"Comparison of integrated $J$-factors for $p$-wave annihilation", size=28)
    outfile = './total_j_factors/j_factors/comparison_plot_p.pgf'
    fig.savefig(outfile, bbox_inches='tight')

    # s wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))

    data = [js, others_data['bergstrom2017_s'], others_data['petac2018_s'], others_data['zhao2016_s']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5)]
    labels = ['This analysis', label_dict['1712.03188'], label_dict['1804.05052'], label_dict['1601.02181']]

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of total integrated $J$-factors for $s$-wave annihilation", size=28)
    outfile = './total_j_factors/j_factors/comparison_plot_s_total.pgf'
    fig.savefig(outfile, bbox_inches='tight')

    # som wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))

    data = [jsom, others_data['bergstrom2017_som'], others_data['petac2018_som'], others_data['kumar2017_som'], others_data['zhao2016_som']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), 'xkcd:purple']
    labels = ['This analysis', label_dict['1712.03188'], label_dict['1804.05052'], label_dict['1702.00408'], label_dict['1601.02181']]

    compare_plot(ax, data, colors, labels)
    ax.set_title(r"Comparison of total integrated $J$-factors for Sommerfeld-enhanced annihilation", size=28)
    outfile = './total_j_factors/j_factors/comparison_plot_som_total.pgf'
    fig.savefig(outfile, bbox_inches='tight')



if __name__ == '__main__':
    bounds = [0.5, 10]
    waves = ['s', 'p', 'd', 'som']

    angledict = {
        0.5: 'half',
        0.1: 'tenth',
        0.2: 'twotenth',
        10: '10',
        25: '25',
        50: '50'
    }

    plot_all_dwarfs([f"./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_{angledict[ub]}.npy" for ub in bounds])
    # sorted_j_factors(["./j_factors/tot_j_factors/tot_j_fac_inclusive_" + angledict[ub] + ".npy" for ub in [0.5, 10]])
    compare_to_others_plot()
