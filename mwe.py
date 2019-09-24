#!/usr/bin/env python3

# import matplotlib.pyplot as plt
# import matplotlib as mpl
#
# mpl.use('pgf')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
#
# fig, ax = plt.subplots()
#
# ax.set_xticklabels(["foo", "bar", "quux", "longer"], rotation=90)
#
# fig.savefig('test.pgf', bbox_inches='tight')

import matplotlib.pyplot as plt
import matplotlib as mpl
# from matplotlib import font_manager
import os

# fm = font_manager.json_load(os.path.expanduser("~/.matplotlib/fontlist-v310.json"))
# fm.findfont("serif", rebuild_if_missing=True)

# font_manager._rebuild()
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

# font_dict = {'family': 'serif', 'serif': ['Computer Modern Roman']}
#
# plt.rc('text', usetex=True)
# plt.rc('font', **font_dict)
mpl.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
mpl.rc('text', usetex=True)

# prop = fm.FontProperties(fname='/usr/local/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/cmr10.ttf')

# mpl.rcParams['font.serif'] = 'DejaVu Serif'
# plt.rcParams["font.family"] = ""
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'xx-large'
mpl.rcParams['axes.titlesize'] = 'xx-large'

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

mpl.use("pgf")
# plt.xticks(range(6), ["foo", "longer", "bar", "short", "quux", "very long"], rotation=90, fontproperties=prop)
plt.xticks(range(6), ["foo", "longer", "bar", "short", "quux", "very long"], rotation=90)
plt.savefig("test.pgf")
