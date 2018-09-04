# -*- coding: utf-8 -*-
"""
General useful functions

@author: Tom Williams

v0.00.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd

from scipy.constants import h,k,c

def read_sed(isrf,
             alpha,
             sCM20_df,
             lCM20_df,
             aSilM5_df,
             frequency):
    
    #Read in the SED that corresponds to this
    #combination of ISRF strength and alpha 
    #(rounded to nearest 0.01)
    
    if '%.2f' % isrf in ['-0.00']:
        isrf = np.abs(isrf)

    small_grains = sCM20_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    large_grains = lCM20_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    silicates = aSilM5_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    
    #Convert these into units more like Jy -- divide by Hz
    
    small_grains /= frequency
    large_grains /= frequency
    silicates /= frequency
    
    return small_grains,large_grains,silicates

def define_stars(flux_df,
                 gal_row,
                 frequency):
    
    #Create a blackbody of 5000K to represent the stars, and lock this
    #at the 3.6micron IRAC flux (or the 3.4 WISE flux if not available)
    
    stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1
    
    wavelength = 3e8/frequency
    wavelength *= 1e6
    
    try:
        if not np.isnan(flux_df['Spitzer_3.6'][gal_row]):
            idx = np.where(np.abs(wavelength-3.6) == np.min(np.abs(wavelength-3.6)))
            ratio = flux_df['Spitzer_3.6'][gal_row]/stars[idx]
        else:
            idx = np.where(np.abs(wavelength-3.4) == np.min(np.abs(wavelength-3.4)))
            ratio = flux_df['WISE_3.4'][gal_row]/stars[idx]
    except KeyError:
        idx = np.where(np.abs(wavelength-3.4) == np.min(np.abs(wavelength-3.4)))
        ratio = flux_df['WISE_3.4'][gal_row]/stars[idx]
    
    stars *= ratio
    
    return stars