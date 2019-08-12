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


def configure_plot_params(fontsize=20, spacing=0.2):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    mpl.rcParams['font.size'] = fontsize
    mpl.rcParams['axes.labelsize'] = 'x-large'
    mpl.rcParams['axes.titlesize'] = 'xx-large'

    mpl.rcParams['figure.subplot.wspace'] = spacing
    mpl.rcParams['figure.subplot.hspace'] = spacing

    mpl.rcParams['xtick.labelsize'] = 'x-large'
    mpl.rcParams['xtick.major.size'] = 7.5
    mpl.rcParams['xtick.major.width'] = 1
    mpl.rcParams['xtick.minor.size'] = 3.75
    mpl.rcParams['xtick.minor.width'] = 0.5

    mpl.rcParams['ytick.labelsize'] = 'x-large'
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


def error_plot(axs, data, offset=0, color='r', label=None, ec=0, mec=None, ms=6, lw=2):
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

    axs.errorbar(x, ave, yerr=[low, high], fmt='D', ecolor=ec, mec=mec, color=color, mew=lw, capsize=lw+5, label=label, ms=ms, lw=lw)


def fix_names(name_list, get_dict=False):
    """Fix names for plotting of dwarfs."""
    name_dict = {
    'ursamajor2': 'U. Major II',
    'willman1': 'Willman I',
    'horologium1_k15': 'Hor. I',
    'segue1': 'Segue I',
    'comaberenices': 'Com. Ber.',
    'tucana2_des': 'Tucana II',
    'reticulum2_k15': 'Ret. II',
    'draco1': 'Draco I',
    'hydrus1': 'Hydrus I',
    'ursaminor': 'U. Minor',
    'sculptor': 'Sculptor',
    'carina2': 'Carina II',
    'aquarius2': 'Aquarius II',
    'bootes1': 'Bootes I',
    'ursamajor1': 'U. Major I',
    'fornax': 'Fornax',
    'canesvenatici2': 'Can. Ven. II',
    'carina': 'Carina',
    'sextans': 'Sextans',
    'leo2': 'Leo II',
    'leo1': 'Leo I',
    'canesvenatici1': 'Can. Ven. I',
    'sagittarius2': 'Sag. II',
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
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5)]
    labels = [r"$s$-wave, ", r"$p$-wave, ", r"$d$-wave, ", r"Sommerfeld, "]

    configure_plot_params()

    fig, ax = plt.subplots(figsize=(25, 25))

    for errorcolor, angle_jfacs, a in zip(['xkcd:slate grey', 'xkcd:black'], infile, [r"$0.5^\circ$", r"$10^\circ$"]):
    # for errorcolor, angle_jfacs, a in zip(['xkcd:black'], infile, [r"$10^\circ$"]):
        jfacs = np.load(angle_jfacs, allow_pickle=True)
        jfacs[:, 0] = set_x_value_names(fix_names(jfacs[:, 0]))

        js = jfacs[np.arange(0, len(jfacs), 4)]
        jp = jfacs[np.arange(1, len(jfacs), 4)]
        jd = jfacs[np.arange(2, len(jfacs), 4)]
        jsom = jfacs[np.arange(3, len(jfacs), 4)]

        for d, c, l in zip([js, jp, jd, jsom], colors, labels):
            error_plot(ax, d, color=c, label=l+a, ec=errorcolor, mec=errorcolor, ms=12, lw=2.7)

    for idwarf in np.arange(0, len(js)+2, 2):
        ax.axvline(idwarf, color='xkcd:light gray', lw=48, zorder=-100)

    x_labels = set_x_value_names(None, get_dict=True)

    ax.set_xticks(list(x_labels.values()))
    ax.set_xticklabels(list(x_labels.keys()))

    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', left=True, right=True)
    ax.xaxis.set_tick_params(which='minor', bottom=False)
    plt.xticks(rotation=90)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    plt.ylabel(r"$\log_{10}\,\,J/(\mathrm{GeV}^2\mathrm{cm}^{-5})$")
    plt.legend(markerscale=1.25, prop={'size': 24}, frameon=True, loc='lower right')

    outfile = infile[0][:-9] + '_working_plot.pdf'
    fig.savefig(outfile, bbox_inches='tight')


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

    my_data = "./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_half.npy"
    jfacs = np.load(my_data, allow_pickle=True)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]

    others_data = np.load('./total_j_factors/j_factors/compare_values_j_factors.npz')

    # s wave, to 0.5 degree
    data = [js, others_data['geringer2015'], others_data['pace2018'], others_data['walker2015'], others_data['caldwell2016'], others_data['koposov2018']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5), 'xkcd:purple', 'xkcd:hot pink']
    labels = ['Our analysis', '1408.0002', '1802.06811', '1511.06296', '1612.06398', '1804.06430']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of $J$-factors at $0.5^\circ$ for $s$-wave annihilation")
    outfile = './total_j_factors/j_factors/comparison_plot_s_half.pdf'
    fig.savefig(outfile, bbox_inches='tight')

    # s wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))
    my_data = "./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_10.npy"
    jfacs = np.load(my_data)
    jfacs[:, 0] = fix_names(jfacs[:, 0])
    js = jfacs[np.arange(0, len(jfacs), 4)]
    jsom = jfacs[np.arange(3, len(jfacs), 4)]

    data = [js, others_data['bergstrom2017_s'], others_data['petac2018_s']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach']
    labels = ['Our analysis', '1712.03188', '1804.05052']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of total integrated $J$-factors for $s$-wave annihilation")
    outfile = './total_j_factors/j_factors/comparison_plot_s_total.pdf'
    fig.savefig(outfile, bbox_inches='tight')

    # som wave, total
    plt.clf()
    fig, ax = plt.subplots(figsize=(20, 20))

    data = [jsom, others_data['bergstrom2017_som'], others_data['petac2018_som'], others_data['kumar2017_som']]
    colors = ['xkcd:azure', 'xkcd:coral', 'xkcd:peach', lighten_color('xkcd:light turquoise', amount=1.5)]
    labels = ['Our analysis', '1712.03188', '1804.05052', '1702.00408']

    compare_plot(ax, data, colors, labels)
    plt.title(r"Comparison of total integrated $J$-factors for Sommerfeld-enhanced annihilation")
    outfile = './total_j_factors/j_factors/comparison_plot_som_total.pdf'
    fig.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    bounds = [0.5, 10]
    angledict = {
        0.5: 'half',
        10: '10',
        25: '25',
        50: '50'
    }
    plot_all_dwarfs([f"/Users/runburg/github/dm_halos/total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_{angledict[ub]}.npy" for ub in bounds])
    # sorted_j_factors("./j_factors/tot_j_factors/tot_j_fac_inclusive_" + angledict[10] + ".npy")
    compare_to_others_plot()
