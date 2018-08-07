# -*- coding: utf-8 -*-
"""
SED fitter functionality tests

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import matplotlib
matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
import numpy as np
from scipy.constants import h,k,c
import matplotlib.pyplot as plt
import os
import emcee
import corner
import time
import multiprocessing
from collections import OrderedDict
import tarfile

start_time = time.time()

#Tar file

os.chdir('/home/daedalusdata/c1625914/dustem_fitter')

tar = tarfile.open('dustem_grid.tar.gz')
print(tar.extractall())

alpha = 3.4

for alpha in np.arange(2.6,5.4,0.01):

    wavelength,\
        small_grains,\
        large_grains,\
        pyroxine,\
        olivine, = np.loadtxt(tar.extractfile('SED_%.2f.RES' % alpha),
                              unpack=True,
                              skiprows=8,
                              usecols=(0,1,2,3,4))
        
    print(wavelength)

print('Complete!')