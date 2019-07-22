#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make bar and whisker plots of all total j-factors.

Author: Jack Runburg
Date: 26-06-2019 17:08

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


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


def error_plot(axs, data, offset=0, color='r', label=None, ec=0, mec=None):
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
        x = data[:, 0]
    except IndexError:
        x = data[0].astype(np.float) + offset

    axs.errorbar(x, ave, yerr=[low, high], fmt='D', ecolor=ec, mec=mec, color=color, mew=1.5, capsize=5, elinewidth=2, label=label)


def fix_names(name_list, get_dict=False):
    """Fix names for plotting of dwarfs."""
    name_dict = {
    'ursamajor2': 'Ursa Major II',
    'willman1': 'Willman I',
    'horologium1_k15': 'Horologium I',
    'segue1': 'Segue I',
    'comaberenices': 'Coma Berenices',
    'tucana2_k15': 'Tucana II',
    'reticulum2_k15': 'Reticulum II',
    'draco1': 'Draco I',
    'hydrus1': 'Hydrus I',
    'ursaminor': 'Ursa Minor',
    'sculptor': 'Sculptor',
    'carina2': 'Carina II',
    'aquarius2': 'Aquarius II',
    'bootes1': 'Bootes I',
    'ursamajor1': 'Ursa Major I',
    'fornax': 'Fornax',
    'canesvenatici2': 'Canes Venatici II',
    'carina': 'Carina',
    'sextans': 'Sextans',
    'leo2': 'Leo II',
    'leo1': 'Leo I',
    'canesvenatici1': 'Canes Venatici I',
    'sagittarius2': 'Sagittarius II',
    'hercules': 'Hercules',
    'crater2': 'Crater II'}

    if get_dict:
        return name_dict
    else:
        names = [name_dict[item] for item in name_list]
        return names


def set_x_value_names(data, get_dict=False):
    """For a given dwarf name, return the coordinate it corresponds to."""
    coord_dict = {item: num for (num, item) in enumerate(sorted(fix_names(None, True).values()), start=1)}
    if get_dict:
        return coord_dict
    else:
        return [coord_dict[item] for item in data]


def plot_all_dwarfs(infile):
    """Generate the plots for all the dwarfs in the given file."""
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', 'xkcd:light turquoise']
    labels = [r"$s$-wave, ", r"$p$-wave, ", r"$d$-wave, ", r"Sommerfeld, "]

    fig, ax = plt.subplots(figsize=(25, 25))

    for errorcolor, angle_jfacs, a in zip(['xkcd:slate grey', 'xkcd:black'], infile, [r"$0.5^\circ$", r"$10^\circ$"]):
    # for errorcolor, angle_jfacs, a in zip(['xkcd:black'], infile, [r"$10^\circ$"]):
        jfacs = np.load(angle_jfacs)
        jfacs[:, 0] = fix_names(jfacs[:, 0])

        js = jfacs[np.arange(0, len(jfacs), 4)]
        jp = jfacs[np.arange(1, len(jfacs), 4)]
        jd = jfacs[np.arange(2, len(jfacs), 4)]
        jsom = jfacs[np.arange(3, len(jfacs), 4)]

        for d, c, l in zip([js, jp, jd, jsom], colors, labels):
            error_plot(ax, d, color=c, label=l+a, ec=errorcolor, mec=errorcolor)

    for idwarf in np.arange(0, len(js), 2):
        ax.axvline(js[idwarf, 0], color='xkcd:light gray', lw=48, zorder=-100)

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', left=True, right=True)
    ax.xaxis.set_tick_params(which='minor', bottom=False)
    plt.xticks(rotation=90, fontsize=15)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    mpl.rcParams['axes.labelsize'] = 'xx-large'
    mpl.rcParams['axes.titlesize'] = 80
    mpl.rcParams['ytick.labelsize'] = 'x-large'
    mpl.rcParams['xtick.labelsize'] = 'x-large'

    plt.ylabel(r"$\log_{10}\,\,\,\,J/(\mathrm{GeV}^2\mathrm{cm}^{-5})$", fontsize=20)
    plt.legend(markerscale=2, prop={'size': 15}, frameon=True, loc='lower left')

    outfile = infile[0][:-9] + '_working_plot.pdf'
    fig.savefig(outfile, bbox_inches='tight')


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


def compare_plot(ax, data, colors, labels):
    """Generate the plot for the comparison to other results."""
    start = len(data)//2-len(data)+len(data) % 2
    for ioffset, (d, c, l) in enumerate(zip(data, colors, labels), start=start):
        try:
            d[:, 0] = set_x_value_names(d[:, 0])
        except IndexError:
            d[0] = set_x_value_names([d[0]])[0]
        except KeyError:
            pass
        error_plot(ax, d, offset=ioffset/(2*len(data)), color=c, label=l, ec=lighten_color(c, amount=1.5), mec=lighten_color(c, amount=1.5))

    mpl.rcParams['axes.labelsize'] = 'xx-large'
    mpl.rcParams['axes.titlesize'] = 25
    mpl.rcParams['ytick.labelsize'] = 'x-large'
    mpl.rcParams['xtick.labelsize'] = 'x-large'

    x_labels = set_x_value_names(None, get_dict=True)

    ax.set_xticks(list(x_labels.values()))
    ax.set_xticklabels(list(x_labels.keys()))
    ax.xaxis.set_tick_params(which='minor', bottom=False)
    plt.xticks(rotation=90, fontsize=15)

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', left=True, right=True)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    for idwarf in np.arange(1, len(data[0])+2, 2):
        ax.axvline(idwarf, color='xkcd:light grey', lw=40, zorder=-100)

    plt.ylabel(r"$\log_{10}\,\,\,\,J/(\mathrm{GeV}^2\mathrm{cm}^{-5})$", fontsize=20)
    plt.legend(markerscale=2, prop={'size': 15}, frameon=True, loc='lower left')


def compare_to_others_plot():
    """Generate plots to compare to other papers."""
    fig, ax = plt.subplots(figsize=(20, 20))

    my_data = "./j_factors/tot_j_factors/tot_j_fac_inclusive_half.npy"
    jfacs = np.load(my_data)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]

    others_data = np.load('./j_factors/compare_values_j_factors.npz')

    # s wave, to 0.5 degree
    data = [js, others_data['geringer2015'], others_data['pace2018'], others_data['walker2015'], others_data['caldwell2016'], others_data['koposov2018']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), 'xkcd:purple', 'xkcd:hot pink']
    labels = ['Our analysis', '1408.0002', '1802.06811', '1511.06296', '1612.06398', '1804.06430']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of $J$-factors at $0.5^\circ$ for $s$-wave annihilation")
    outfile = './j_factors/comparison_plot_s_half.pdf'
    fig.savefig(outfile, bbox_inches='tight')

    # s wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))
    my_data = "./j_factors/tot_j_factors/tot_j_fac_inclusive_10.npy"
    jfacs = np.load(my_data)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]
    jsom = jfacs[np.arange(3, len(jfacs), 4)]

    data = [js, others_data['bergstrom2017_s'], others_data['petac2018_s']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach']
    labels = ['Our analysis', '1712.03188', '1804.05052']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of total integrated $J$-factors for $s$-wave annihilation")
    outfile = './j_factors/comparison_plot_s_total.pdf'
    fig.savefig(outfile, bbox_inches='tight')

    # som wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))

    data = [jsom, others_data['bergstrom2017_som'], others_data['petac2018_som'], others_data['kumar2017_som']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5)]
    labels = ['Our analysis', '1712.03188', '1804.05052', '1702.00408']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of total integrated $J$-factors for Sommerfeld-enhanced annihilation")
    outfile = './j_factors/comparison_plot_som_total.pdf'
    fig.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    bounds = [0.5, 10]
    angledict = {
        0.5: 'half',
        10: '10',
        25: '25',
        50: '50'
    }
    plot_all_dwarfs([f"/Users/runburg/github/dm_halos/j_factors/tot_j_factors/tot_j_fac_inclusive_{angledict[ub]}.npy" for ub in bounds])
    # sorted_j_factors("./j_factors/tot_j_factors/tot_j_fac_inclusive_" + angledict[10] + ".npy")
    compare_to_others_plot()
