# -*- coding: utf-8 -*-
"""
Prepare data for fitting THEMIS model

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import matplotlib
matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

#Read in the M33 fluxes

os.chdir('/home/daedalusdata/c1625914/SKIRT/run')

obs_wavelength = np.loadtxt('filters.dat',usecols=1)
obs_frequency = 3e8 / (obs_wavelength*1e-6)

obs_flux = []
obs_err = []

obs_file = open('m33.dat','r')
obs_file.readline()

line = obs_file.readline()
line = line.strip()
columns = line.split()

for i in range(2,len(columns),2):
    obs_flux.append(float(columns[i]))
    obs_err.append(float(columns[i+1]))
    
obs_flux = np.array(obs_flux)
obs_err = np.array(obs_err)

#Also an additional 10% errorbar

extra_err = obs_err+0.1*obs_flux

os.chdir('/home/daedalusdata/c1625914/dustem/out/')

#Read in the SED

wavelength,\
    small_grains,\
    large_grains,\
    pyroxine,\
    olivine, = np.loadtxt('SED.RES',
                          unpack=True,
                          skiprows=8,
                          usecols=(0,1,2,3,4))
    
silicates = pyroxine+olivine
    
#Convert these into units more like Jy -- divide by Hz

frequency = 3e8/(wavelength*1e-6)

small_grains /= frequency
large_grains /= frequency
silicates /= frequency

#Add in a blackbody @ 5000K to represent the stars

h = 6.63e-34
k = 1.38e-23
c = 3e8

stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1

#Lock this at 1e-38@3.6microns (as a first guess)

themis_normalise_idx = np.where(np.abs(wavelength-3.6) == np.min(np.abs(wavelength-3.6)))
ratio = stars[themis_normalise_idx]/1e-38

stars /= ratio

total = small_grains+large_grains+silicates+stars

#Only use dust -- i.e. 3.4 microns onwards

idx = np.where(obs_wavelength >= 3.368)

obs_wavelength = obs_wavelength[idx]
obs_flux = obs_flux[idx]
extra_err = extra_err[idx]

err_ratio = extra_err/obs_flux

filters = ['WISE3_4','IRAC3_6','IRAC4_5',
           'WISE4_6','IRAC5_8','IRAC8','WISE12',
           'WISE22','MIPS24','MIPS70','PACS100',
           'PACS160','SPIRE250','SPIRE350',
           'SPIRE500','planck353',
           'SCUBA2_450']

d = {'name':'M33',
     'dist':840}

for i in range(len(filters)):
    try:
        d[filters[i]] = [obs_flux[i]]
        d[filters[i]+'_err'] = [extra_err[i]]
    except:
        d[filters[i]] = ['']
        d[filters[i]+'_err'] = ['']

df = pd.DataFrame(data=d)

df.to_csv('/home/daedalusdata/c1625914/dustem_fitter/fluxes.csv')

#Normalise at 250microns

obs_normalise_idx = np.where(obs_wavelength == 250)
themis_normalise_idx = np.where(np.abs(wavelength-250) == np.min(np.abs(wavelength-250)))

ratio = total[themis_normalise_idx]/obs_flux[obs_normalise_idx]

# obs_flux *= ratio
extra_err = obs_flux*err_ratio

#Turn into Jy

# total /= 1e-23
# total /= 4*np.pi
# total /= (840*1000*3.1e18)**2

#Write out these modified fluxes and errors

np.savetxt('themis_grid/obs_flux.txt',
           np.c_[obs_wavelength,obs_flux,extra_err])
    
plt.figure()

#Plot on the actual (normalised) data

plt.errorbar(obs_wavelength,obs_flux,
             yerr=extra_err,c='k',
             ls='none',
             marker='o')

#Plot on the THEMIS model

plt.plot(wavelength,small_grains,'b',
         label='sCM20')
plt.plot(wavelength,large_grains,'g',
         label='lCM20')
plt.plot(wavelength,silicates,'r',
         label='aSilM5')
plt.plot(wavelength,stars,'m',
         label='Stars')
plt.plot(wavelength,total,'k',
         label='Total')

plt.xscale('log')
plt.yscale('log')

plt.xlim([1,1000])
plt.ylim([1e-40,1e-35])

plt.xlabel(r'$\lambda$ ($\mu$m)')
plt.ylabel(r'erg/s/H/Hz')

plt.legend(loc='upper left',
           frameon=False)

plt.show()

print('Complete!')