#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Format spectra and j factors for MADHAT.

Author: Jack Runburg
Date: 21-07-2019 17:23

"""

import numpy as np

def fix_dwarf_names(names):
    # from MADHAT
    dwarfs = ['Bootes I',
    'Bootes II',
    'Bootes III',
    'Canes Venatici I',
    'Canes Venatici II',
    'Carina',
    'Cetus II',
    'Columba I',
    'Coma Berenices',
    'Draco',
    'Draco II',
    'Eridanus II',
    'Eridanus III',
    'Fornax',
    'Grus I',
    'Grus II',
    'Hercules',
    'Horologium I',
    'Horologium II',
    'Hydra II',
    'Indus II',
    'Kim 2',
    'Leo I',
    'Leo II',
    'Leo IV',
    'Leo T',
    'Leo V',
    'Pegasus III',
    'Phoenix II',
    'Pictor I',
    'Pisces II',
    'Reticulum II',
    'Reticulum III',
    'Sagittarius II',
    'Sculptor',
    'Segue 1',
    'Segue 2',
    'Sextans',
    'Triangulum II',
    'Tucana II',
    'Tucana III',
    'Tucana IV',
    'Tucana V',
    'Ursa Major I',
    'Ursa Major II',
    'Ursa Minor',
    'Willman I']
    # calculated j-factors
    name_dict = {
    'ursamajor2': 'Ursa Major II',
    'willman1': 'Willman I',
    'horologium1_k15': 'Horologium I',
    'segue1': 'Segue I',
    'comaberenices': 'Coma Berenices',
    'tucana2_des': 'Tucana II',
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
    dwarf_dict = {dwarf: i+1 for i, dwarf in enumerate(dwarfs)}

    return [[dwarf_dict[name_dict[name]], *names[i, 2:].astype(np.float)] for i, name in enumerate(names[:, 0]) if name_dict[name] in dwarf_dict.keys()]


spectra = np.load('/Users/runburg/github/dm_halos/gamma_spectrums/integrated_spectra.npy')
choose_channel = spectra.dtype.names

for choose in choose_channel:
    with open(f'/Users/runburg/github/dm_halos/gamma_spectrums/madhat_channel_spectra/channel_{choose}.txt', 'w') as outfile:
        outfile.write(f'# {choose} channel integrated photon spectrum\n')
        outfile.write('MASS\tEnergy Spectrum\n')
        for mass, spectrum in zip(spectra['mDM'], spectra[choose]):
            outfile.write(f'{mass}\t{spectrum}\n')

jfacs = np.load('/Users/runburg/github/dm_halos/j_factors/tot_j_factors/tot_j_fac_inclusive_10.npy')

jfacs = np.array(fix_dwarf_names(jfacs))

js = jfacs[np.arange(0, len(jfacs), 4)]
jp = jfacs[np.arange(1, len(jfacs), 4)]
jd = jfacs[np.arange(2, len(jfacs), 4)]
jsom = jfacs[np.arange(3, len(jfacs), 4)]

for j, wave in zip([js, jp, jd, jsom], ['s', 'p', 'd', 'som']):
    ave = j[:, 2].astype(np.float)
    low = ave - j[:, 1].astype(np.float)
    high = j[:, 3].astype(np.float) - ave
    with open(f'./madhat/models/j_{wave}.txt', 'w') as outfile:
        outfile.write(f'# J_{wave} factors for dwarves in madhat\n')
        outfile.write('ID\tJ\tdJ+\tdJ-\n')
        for dwarf in zip(j[:, 0].astype(int), ave, high, low):
            outfile.write('{}\t{}\t{}\t{}\n'.format(*dwarf))

jfacs = np.load('/Users/runburg/github/dm_halos/j_factors/tot_j_factors/tot_j_fac_inclusive_10.npy')

jfacs = np.array(fix_dwarf_names(jfacs))

js = jfacs[np.arange(0, len(jfacs), 4)]
jp = jfacs[np.arange(1, len(jfacs), 4)]
jd = jfacs[np.arange(2, len(jfacs), 4)]
jsom = jfacs[np.arange(3, len(jfacs), 4)]

for j, wave in zip([js, jp, jd, jsom], ['s', 'p', 'd', 'som']):
    ave = j[:, 2].astype(np.float)
    print(ave)
    indices = [10**avee > 0.15 * 10**ave.max() for avee in ave]
    low = ave - j[:, 1].astype(np.float)
    high = j[:, 3].astype(np.float) - ave

    with open(f'./madhat/models/j_{wave}_15.txt', 'w') as outfile:
        outfile.write(f'# J_{wave} factors for dwarves in madhat\n')
        outfile.write('ID\tJ\tdJ+\tdJ-\n')
        for dwarf in zip(j[:, 0][indices].astype(int), ave[indices], high[indices], low[indices]):
            outfile.write('{}\t{}\t{}\t{}\n'.format(*dwarf))
