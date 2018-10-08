"""
Compress each dust component

@author: Tom Williams
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import pandas as pd
import os
import glob
import numpy as np
from tqdm import tqdm

os.chdir(os.getcwd())

#First, create the names

res_names = glob.glob('grid/*.RES')

wavelength = np.loadtxt(res_names[0],
                        unpack=True,
                        skiprows=8,
                        usecols=(0))

names = []

for res_name in res_names:
    
    res_name = res_name.split('/')[1].split('.RES')[0]
    columns = res_name.split('_')
    names.append('alpha_sCM20:'+columns[1]+',logU:'+columns[2])
        
#Set up a pandas DataFrame for each component (small- and large-carbons,
#and silicates), and save wavelength separately

frequency = 3e8/(wavelength*1e-6)

wavelength_df = pd.DataFrame(columns=['wavelength'],data=wavelength)

sCM20_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))
lCM20_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))
aSilM5_df = pd.DataFrame(columns=names,data=np.zeros([len(wavelength),len(names)]))

#Read in each grid and put into the dataframe

for res_name in tqdm(res_names,
              desc='Total'):
      
    small_grains,\
        large_grains,\
        pyroxine,\
        olivine, = np.loadtxt(res_name,
                              unpack=True,
                              skiprows=8,
                              usecols=(1,2,3,4))
      
    res_name_split = res_name.split('/')[1].split('.RES')[0]
    columns = res_name_split.split('_')
    identifier = 'alpha_sCM20:'+columns[1]+',logU:'+columns[2]
          
    sCM20_df[identifier] = small_grains/frequency
    lCM20_df[identifier] = large_grains/frequency
    aSilM5_df[identifier] = (pyroxine+olivine)/frequency
 
sCM20_df.to_hdf('../core/models.h5',
                mode='w',key='sCM20')
wavelength_df.to_hdf('../core/models.h5',
                     key='wavelength')
lCM20_df.to_hdf('../core/models.h5',
                key='lCM20')
aSilM5_df.to_hdf('../core/models.h5',
                key='aSilM5')

print('Complete!')