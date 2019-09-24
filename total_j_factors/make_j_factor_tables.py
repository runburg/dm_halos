#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create tables of J-factors for paper.

Author: Jack Runburg
Date: 15-08-2019 14:46


"""
import numpy as np


def fix_names(name_list, get_dict=False):
    """Fix names for plotting of dwarfs."""
    name_dict = {
        'ursamajor2': 'Ursa Major II',
        'willman1': 'Willman 1',
        'horologium1_k15': 'Horologium I',
        'segue1': 'Segue 1',
        'comaberenices': 'Coma Berenices I',
        'tucana2_des': 'Tucana II',
        'reticulum2_k15': 'Reticulum II',
        'draco1': 'Draco I',
        'hydrus1': 'Hydrus I',
        'ursaminor': 'Ursa Minor I',
        'sculptor': 'Sculptor I',
        'carina2': 'Carina II',
        'aquarius2': 'Aquarius II',
        'bootes1': r'Bo\"{o}tes I',
        'ursamajor1': 'Ursa Major I',
        'fornax': 'Fornax I',
        'canesvenatici2': 'Canes Venatici II',
        'carina': 'Carina I',
        'sextans': 'Sextans I',
        'leo2': 'Leo II',
        'leo1': 'Leo I',
        'canesvenatici1': 'Canes Venatici I',
        'sagittarius2': 'Sagittarius II',
        'hercules': 'Hercules I',
        'crater2': 'Crater II'}

    names = [name_dict[item] for item in name_list]
    return names


angles = ['tenth', 'twotenth', 'half', '10']
header = (r"\centering\footnotesize\hspace*{-2cm}\begin{tabular}{| l | c c c c | c c c c | c c c c | c c c c |}"
          r"\cline{2-17}\multicolumn{1}{c}{} & \multicolumn{4}{|c|}{$s$-wave} & \multicolumn{4}{|c|}{$p$-wave} & \multicolumn{4}{|c|}{$d$-wave} & \multicolumn{4}{|c|}{Sommerfeld}\\\hline"
          r" Galaxy & $0.1^\circ$ & $0.2^\circ$ & $0.5^\circ$ & $10^\circ$"
          r"& $0.1^\circ$ & $0.2^\circ$ & $0.5^\circ$ & $10^\circ$"
          r"& $0.1^\circ$ & $0.2^\circ$ & $0.5^\circ$ & $10^\circ$"
          r"& $0.1^\circ$ & $0.2^\circ$ & $0.5^\circ$ & $10^\circ$\\\hline %")
footer = r"\hline \end{tabular}%"

js, jp, jd, jsom = [], [], [], []
jfacs = np.load(f'./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_tenth.npy', allow_pickle=True)
table = np.empty((len(jfacs)//4, 17), dtype='U24')
table[:, 0] = fix_names(jfacs[np.arange(0, len(jfacs), 4)][:, 0])

for i, jfactor in enumerate([js, jp, jd, jsom]):
    for i_angle, angle in enumerate(angles, start=1):
        jfacs = np.load(f'./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_{angle}.npy', allow_pickle=True)
        j_temp = jfacs[np.arange(i, len(jfacs), 4)]
        j_temp = j_temp[j_temp[:, 0].argsort()]

        table[:, i_angle + 4*i] = [rf"${format(round(ave.astype(np.float),2), '.2f')}^{{+{format(round(high.astype(np.float)-ave.astype(np.float),2), '.2f')}}}_{{-{format(round(ave.astype(np.float)-low.astype(np.float),2), '.2f')}}}$" for low, ave, high in zip(j_temp[:, 2], j_temp[:, 3], j_temp[:, 4])]

np.savetxt('./total_j_factors/j_factors/tables/full_table.txt', table, delimiter=' & ', newline=r'\\'+'\n', header=header, footer=footer, fmt='%s', comments='')


          # for angle in angles:
          #     jfacs = np.load(f'./total_j_factors/j_factors/tot_j_factors/tot_j_fac_inclusive_{angle}.npy', allow_pickle=True)
          #
          #     for i, jfactor in enumerate([js, jp, jd, jsom]):
          #         j_temp = jfacs[np.arange(i, len(jfacs), 4)]
          #         j_temp = j_temp[j_temp[:, 0].argsort()]
          #         if angle is angles[0]:
          #             jfactor.append(fix_names(j_temp[:, 0]))
          #         jfactor.append([rf"${format(round(ave.astype(np.float),2), '.2f')}^{{+{format(round(high.astype(np.float)-ave.astype(np.float),2), '.2f')}}}_{{-{format(round(ave.astype(np.float)-low.astype(np.float),2), '.2f')}}}$" for low, ave, high in zip(j_temp[:, 2], j_temp[:, 3], j_temp[:, 4])])
          #
          # np.savetxt('./total_j_factors/j_factors/tables/J_s_table.txt', np.array(js).T, delimiter=' & ',newline=r'\\'+'\n', header=header, fmt='%s', comments='')
          # np.savetxt('./total_j_factors/j_factors/tables/J_p_table.txt', np.array(jp).T, delimiter=' & ',newline=r'\\'+'\n', header=header, fmt='%s', comments='')
          # np.savetxt('./total_j_factors/j_factors/tables/J_d_table.txt', np.array(jd).T, delimiter=' & ',newline=r'\\'+'\n', header=header, fmt='%s', comments='')
          # np.savetxt('./total_j_factors/j_factors/tables/J_som_table.txt', np.array(jsom).T, delimiter=' & ',newline=r'\\'+'\n', header=header, fmt='%s', comments='')
