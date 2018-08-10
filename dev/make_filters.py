# -*- coding: utf-8 -*-
"""
Create a Pandas dataframe for the filters

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import pandas as pd
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/home/daedalusdata/c1625914/dustem_fitter/')

#Create a dataframe with 1) central wavelength and 2) uncorrelated calibration uncertainty

d = {'Spitzer_3.6':[3.6,0.015],
     'Spitzer_4.5':[4.5,0.015],
     'Spitzer_5.8':[5.8,0.015],
     'Spitzer_8.0':[8.0,0.015],
     'Spitzer_24':[24.0,0.004],
     'Spitzer_70':[70.0,0.045],
     'Spitzer_160':[160.0,0.05],
     'WISE_3.4':[3.368,0.029],
     'WISE_4.6':[4.618,0.034],     
     'WISE_12':[12.082,0.046],
     'WISE_22':[22.194,0.056],
     'PACS_70':[70.0,0.02],
     'PACS_100':[100.0,0.02],
     'PACS_160':[160.0,0.02],
     'SPIRE_250':[250.0,0.015],
     'SPIRE_350':[350.0,0.015],
     'SPIRE_500':[500.0,0.015],
     'Planck_350':[350.0,0.064],
     'Planck_550':[550.0,0.061],
     'Planck_850':[850,0.0078],
     'SCUBA2_450':[450.0,0.12],
     'SCUBA2_850':[850,0.08],
     'IRAS_12':[12,0.2],
     'IRAS_25':[25,0.2],
     'IRAS_60':[60,0.2],
     'IRAS_100':[100,0.2]}

df = pd.DataFrame(data=d)

df.to_csv('filters.csv')

#Also pull out the Planck 350 and 550 filters and save those

planck_hdu = fits.open('dev/HFI_RIMO_R1.10.fits')[6]

planck_freq = []
planck_transmission = []

for i in range(planck_hdu.data.shape[0]):
    planck_freq.append(planck_hdu.data[i][0])
    planck_transmission.append(planck_hdu.data[i][1])
    
planck_freq = np.array(planck_freq)
planck_transmission = np.array(planck_transmission)

planck_freq *= 3e8*1e-7
planck_wavelength = (3e8*1e10)/(planck_freq*1e9)

np.savetxt('filters/Planck_550.dat',np.c_[planck_wavelength[1:],planck_transmission[1:]])

planck_hdu = fits.open('dev/HFI_RIMO_R1.10.fits')[7]

planck_freq = []
planck_transmission = []

for i in range(planck_hdu.data.shape[0]):
    planck_freq.append(planck_hdu.data[i][0])
    planck_transmission.append(planck_hdu.data[i][1])
    
planck_freq = np.array(planck_freq)
planck_transmission = np.array(planck_transmission)

planck_freq *= 3e8*1e-7
planck_wavelength = (3e8*1e10)/(planck_freq*1e9)

np.savetxt('filters/Planck_350.dat',np.c_[planck_wavelength[1:],planck_transmission[1:]])

print('Complete!')