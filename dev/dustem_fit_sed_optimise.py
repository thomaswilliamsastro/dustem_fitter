# -*- coding: utf-8 -*-
"""
Optimise the SED fitting

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

start_time = time.time()

os.chdir('/home/daedalusdata/c1625914/dustem/out/themis_grid')

obs_wavelength,obs_flux,obs_error = np.loadtxt('obs_flux.txt',
                                               unpack=True)

scaling = np.abs(np.random.normal(loc=1e40,scale=1e38))
omega_star = np.abs(np.random.normal(loc=1,scale=1e-2))
alpha = np.abs(np.random.normal(loc=5,scale=1e-2*5))
y_sCM20 = np.abs(np.random.normal(loc=1,scale=1e-2))
y_lCM20 = np.abs(np.random.normal(loc=1,scale=1e-2))
y_aSilM5 = np.abs(np.random.normal(loc=1,scale=1e-2))

#Read in the SED that corresponds to this
#alpha (rounded to nearest 0.01)

wavelength,\
    small_grains,\
    large_grains,\
    pyroxine,\
    olivine, = np.loadtxt('SED_%.2f.RES' % alpha,
                          unpack=True,
                          skiprows=8,
                          usecols=(0,1,2,3,4))
    
silicates = pyroxine+olivine

print(time.time()-start_time)
new_time = time.time()

#Convert these into units more like Jy -- divide by Hz

frequency = 3e8/(wavelength*1e-6)

small_grains /= frequency
large_grains /= frequency
silicates /= frequency

#Add in a blackbody @ 5000K to represent the stars

stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1

#Normalise based on earlier value

ratio = 6.97041442e+29

stars /= ratio

#Scale everything!

total = y_sCM20*small_grains+\
        y_lCM20*large_grains+\
        y_aSilM5*silicates+omega_star*stars
        
#Include overall scaling factor

total *= scaling
        
#Convolve this SED with the various filters we have to give
#a monochromatic flux

filter_dict = OrderedDict()

filter_folder = '/home/daedalusdata/c1625914/UsefulCode/sed_fitter/filters/'
filters = ['IRAC3_6','IRAC4_5','WISE4_6','IRAC5_8','IRAC8',
           'WISE12','WISE22','MIPS24','MIPS70','PACS100','PACS160',
           'SPIRE250','SPIRE350','SPIRE500','planck353']

for filter in filters:
    
    filter_wavelength,transmission = np.loadtxt(filter_folder+filter+'.dat',
                                     unpack=True)
    
    filter_wavelength /= 1e4
    transmission /= np.max(transmission)
    
    filter_dict[filter] = filter_wavelength,transmission
    
new_time = time.time()

filter_fluxes = []
    
for key in filter_dict:
            
    #Convolve MBB with filter
            
    filter_flux = np.interp(filter_dict[key][0],wavelength,total)
             
    filter_fluxes.append(np.trapz(filter_dict[key][1]*filter_flux,filter_dict[key][0])/
                        np.trapz(filter_dict[key][1],filter_dict[key][0]))

filter_fluxes = np.abs(np.array(filter_fluxes))

print(time.time()-new_time)
new_time = time.time()

chisq = np.sum( (filter_fluxes-obs_flux)**2 / obs_error**2)

likelihood = -0.5*chisq

print(time.time()-new_time)
new_time = time.time()

print('Complete!')