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

os.chdir('/home/daedalusdata/c1625914/dustem_fitter/dev/grid')

#Set up a pandas DataFrame for each component (small- and large-carbons,
#and silicates)

sCM20_df = pd.DataFrame()
lCM20_df = pd.DataFrame()
aSilM5_df = pd.DataFrame()

alpha = np.arange(2.6,5.41,0.01)
isrf = np.arange(-1,3.51,0.01)

for i in tqdm(alpha,
              desc='Total'):
    for j in tqdm(isrf,
                  leave=False,
                  desc='alpha=%.2f' % i):
    
        small_grains,\
            large_grains,\
            pyroxine,\
            olivine, = np.loadtxt('SED_%.2f_%.2f.RES' % (i,j),
                                  unpack=True,
                                  skiprows=8,
                                  usecols=(1,2,3,4))
            
        if i == alpha[0] and j == isrf[0]:
            sCM20_df['wavelength'] =  np.loadtxt('SED_%.2f_%.2f.RES' % (i,j),
                                                 unpack=True,
                                                 skiprows=8,
                                                 usecols=(0))
    
        identifier = '%.2f,%.2f' % (i,j)
        
        sCM20_df[identifier] = small_grains
        lCM20_df[identifier] = large_grains
        aSilM5_df[identifier] = pyroxine+olivine

sCM20_df.to_hdf('/home/daedalusdata/c1625914/dustem_fitter/dustem_res.h5',
                mode='w',key='sCM20')
lCM20_df.to_hdf('/home/daedalusdata/c1625914/dustem_fitter/dustem_res.h5',
                key='lCM20')
aSilM5_df.to_hdf('/home/daedalusdata/c1625914/dustem_fitter/dustem_res.h5',
                key='aSilM5')

print('Complete!')