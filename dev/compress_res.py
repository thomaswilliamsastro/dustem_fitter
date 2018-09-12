"""
Compress each dust component

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import pandas as pd
import os
import dill
import glob
import numpy as np
from tqdm import tqdm

os.chdir(os.getcwd())

#First, create the names

alpha = np.arange(2.6,5.41,0.01)
isrf = np.arange(-1,3.51,0.01)

names = []

for i in alpha:
    for j in isrf:
        
        names.append('%.2f,%.2f' % (i,j))
        
        #And read in the wavelength
        
        if i == alpha[0] and j == isrf[0]:
            
            wavelength = np.loadtxt('grid/SED_%.2f_%.2f.RES' % (i,j),
                                    unpack=True,
                                    skiprows=8,
                                    usecols=(0))
        
#Set up a pandas DataFrame for each component (small- and large-carbons,
#and silicates)

sCM20_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))
lCM20_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))
aSilM5_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))

#Also save in the wavelength

sCM20_df['wavelength'] = wavelength

#Read in each grid and put into the dataframe

for i in tqdm(alpha,
              desc='Total'):
    for j in tqdm(isrf,
                  leave=False,
                  desc='alpha=%.2f' % i):
     
        small_grains,\
            large_grains,\
            pyroxine,\
            olivine, = np.loadtxt('grid/SED_%.2f_%.2f.RES' % (i,j),
                                  unpack=True,
                                  skiprows=8,
                                  usecols=(1,2,3,4))
     
        identifier = '%.2f,%.2f' % (i,j)
         
        sCM20_df[identifier] = small_grains
        lCM20_df[identifier] = large_grains
        aSilM5_df[identifier] = pyroxine+olivine

sCM20_df.to_hdf('../models.h5',
                mode='w',key='sCM20')
lCM20_df.to_hdf('../models.h5',
                key='lCM20')
aSilM5_df.to_hdf('../models.h5',
                key='aSilM5')

print('Complete!')