# -*- coding: utf-8 -*-
"""
General useful functions

@author: Tom Williams

v1.00.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd

from scipy.constants import h,k,c

def convert_to_luminosity(flux,
                          distance,
                          frequency):
    
    #Get rid of the Hz^-1 term
    
    flux *= frequency
    
    #And the m^2 term -- distance is in Mpc!
    
    flux *= 4*np.pi*(distance*3.086e16*1e6)**2
    
    flux *= 1e-26
    
    #We now have a flux in W so convert to Lsun
    
    flux /= 3.828e26
    
    return flux

def read_sed(isrf,
             alpha,
             sCM20_df,
             lCM20_df,
             aSilM5_df):
    
    #Read in the SED that corresponds to this
    #combination of ISRF strength and alpha 
    #(rounded to nearest 0.01)
    
    if '%.2f' % isrf in ['-0.00']:
        isrf = np.abs(isrf)

    small_grains = sCM20_df['alpha_sCM20:%.2f,logU:%.2f' % (alpha,isrf)].values.copy()
    large_grains = lCM20_df['alpha_sCM20:%.2f,logU:%.2f' % (alpha,isrf)].values.copy()
    silicates = aSilM5_df['alpha_sCM20:%.2f,logU:%.2f' % (alpha,isrf)].values.copy()
    
    return small_grains,large_grains,silicates

def define_stars(flux_df,
                 gal_row,
                 filter_df,
                 frequency,
                 idx_key):
    
    #Create a blackbody of 5000K to represent the stars, and lock this
    #to the shortest wavelength
    
    stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1
    
    wavelength = 3e8/frequency
    wavelength *= 1e6
    
    idx = np.where(np.abs(wavelength-filter_df[idx_key][0]) == np.min(np.abs(wavelength-filter_df[idx_key][0])))
    ratio = flux_df[idx_key][gal_row]/stars[idx]
    
    stars *= ratio
    
    return stars