#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 19:26:52 2019

@author: runburg
"""

import numpy as np
import mpmath as mp


def undim(file, g_n, r_s, rho_s):
    # print('start: reading data from file')

    # set working precision of decimal points
    mp.dps = 25

    '''first we will read the data from the velocity distribution file and
    scale out units'''
    # first column is ignored, second column is f(E), third column is r, fourth
    # column is psi
    file = open(file+".txt", 'r')

    # initial arrays
    fe_init = []
    r_init = []
    psi = []
    v = []

    # scale out units from input file and store column values in ind. arrays
    for line in file.readlines():
        fe_init.append(float(line.split()[1])
                       * 2/rho_s*(rho_s*g_n*r_s**2)**(3.0/2))
        r_init.append(float(line.split()[2])/r_s)
        psi.append(float(line.split()[3])/(rho_s*g_n*r_s**2))

    file.close()

    # print('finished: units scaled out')
    # print('start: change variables (r,psi)->(r,v)')

    '''now we will change variable to r,psi-->r,v'''

    r = []
    fe = []
    v = []
    # initialize counter for keeping track of for loop iterations
    count = 0
    j = 0
    # calculate values of v assuming energy of particle is binding energy psi
    for rad in r_init:
        # set the initial energy value (which will be increased iteratively)
        # to the binding energy
        i = j
        # this loop finds a velocity for increasing values of the energy for a
        # given binding energy
        while i < len(psi):
            # the second term psi[j] is the energy of the particle;
            # the first term psi[i] is the binding energy
            vel = 2 * (psi[i] - psi[j])
            vel = mp.sqrt(vel)
            # data is written as r, v, fe
            v.append(vel)
            r.append(r_init[i])
            fe.append(fe_init[j])
            if j == 0:
                count += 1
            i += 1
        j += 1

    # put everything in increasing order
    dat = np.array(list(zip(r, v, fe)))
    dat = dat[dat[:, 1].argsort()]
    dat = dat[dat[:, 0].argsort(kind='mergesort')]
    r, v, fe = dat.T

    # print('finished: change of variables completed')

    # save new array to outfile
    newfile = file+"_nounits.txt"
    with open(newfile, 'wb') as outfile:
        np.savez(outfile, r=np.array(r), v=np.array(v), fe=np.array(fe))
