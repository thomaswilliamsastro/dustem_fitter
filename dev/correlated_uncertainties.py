# -*- coding: utf-8 -*-
"""
Create a Pandas dataframe for the correlated uncertainties

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import pandas as pd
import numpy as np
import copy

filters = ['IRAC3_6','IRAC4_5','IRAC5_8','IRAC8',
           'MIPS24','MIPS70',
           'PACS70','PACS100','PACS160',
           'SPIRE250','SPIRE350','SPIRE500']

df = pd.DataFrame(0,index=np.arange(len(filters)),columns=filters)

df['name'] = pd.Series(filters)

#N.B all these uncertainties are doubled to account for extended
#sources!

#IRAC uncertainties. From Reach+ (2005), these are 1.8%, 1.9%,
#2.0% and 2.1% at 3.6, 4.5, 5.8 and 8 micron respectively.

idx = df[df['name'].str.contains('IRAC')].index

df['IRAC3_6'][idx] = 0.018*2
df['IRAC4_5'][idx] = 0.019*2
df['IRAC5_8'][idx] = 0.02*2
df['IRAC8'][idx] = 0.021*2

#MIPS uncertainties. From Engelbrach+ (2007) and Gordon+ (2007),
#these are 2% at 24micron and 5% at 70micron

idx = df[df['name'].str.contains('MIPS')].index

df['MIPS24'][idx] = 0.02*2
df['MIPS70'][idx] = 0.05*2

#PACS uncertainties. These are 2% (Muller+ 2011; Balog+ 2014)

idx = df[df['name'].str.contains('PACS')].index 
df['PACS70'][idx] = 0.02*2
df['PACS100'][idx] = 0.02*2
df['PACS160'][idx] = 0.02*2

#SPIRE uncertainties. These are 4% (Bendo+ 2013; Griffin+ 2013)

idx = df[df['name'].str.contains('SPIRE')].index 
df['SPIRE250'][idx] = 0.04*2
df['SPIRE350'][idx] = 0.04*2
df['SPIRE500'][idx] = 0.04*2

print(df)

df.to_csv('/home/daedalusdata/c1625914/dustem_fitter/corr_uncert.csv')

flux = {'IRAC3_6':1,
        'IRAC4_5':1,
        'IRAC5_8':1,
        'IRAC8':1}

corr_err = np.matrix(np.zeros([len(flux),len(flux)]))

keys = []

for key in flux:
    keys.append(key)
 
for i in range(len(flux)):
    for j in range(len(flux)):
         
        corr_err[i,j] = flux[keys[i]]*flux[keys[j]]* \
                            df[keys[j]][df.index[df['name'] == keys[j]][0]]* \
                            df[keys[i]][df.index[df['name'] == keys[i]][0]]