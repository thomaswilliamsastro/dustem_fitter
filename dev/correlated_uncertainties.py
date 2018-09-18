# -*- coding: utf-8 -*-
"""
Create a Pandas dataframe for the correlated uncertainties

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import pandas as pd
import numpy as np
import os

os.chdir(os.getcwd())

filters = ['Spitzer_3.6',
           'Spitzer_4.5',
           'Spitzer_5.8',
           'Spitzer_8.0',
           'Spitzer_24',
           'Spitzer_70',
           'Spitzer_160',
           'WISE_3.4',
           'WISE_4.6',     
           'WISE_12',
           'WISE_22',
           'PACS_70',
           'PACS_100',
           'PACS_160',
           'SPIRE_250',
           'SPIRE_350',
           'SPIRE_500',
           'Planck_350',
           'Planck_550',
           'Planck_850',
           'SCUBA2_450',
           'SCUBA2_850',
           'IRAS_12',
           'IRAS_25',
           'IRAS_60',
           'IRAS_100']

df = pd.DataFrame(0,index=np.arange(len(filters)),columns=filters)

df['name'] = pd.Series(filters)

#N.B all these uncertainties are doubled to account for extended
#sources!

#IRAC uncertainties. From Reach+ (2005), these are 1.8%, 1.9%,
#2.0% and 2.1% at 3.6, 4.5, 5.8 and 8 micron respectively.

idx = df[df['name'].str.contains('Spitzer_3.6|Spitzer_4.5|Spitzer_5.8|Spitzer_8.0')].index

df['Spitzer_3.6'][idx] = 0.018*2
df['Spitzer_4.5'][idx] = 0.019*2
df['Spitzer_5.8'][idx] = 0.02*2
df['Spitzer_8.0'][idx] = 0.021*2

#MIPS uncertainties. From Engelbrach+ (2007) and Gordon+ (2007),
#these are 2% at 24micron and 5% at 70micron

idx = df[df['name'].str.contains('Spitzer_24|Spitzer_70|Spitzer_160')].index
print(idx)

df['Spitzer_24'][idx] = 0.02*2
df['Spitzer_70'][idx] = 0.05*2
df['Spitzer_160'][idx] = 0.05*2

#PACS uncertainties. These are 2% (Muller+ 2011; Balog+ 2014)

idx = df[df['name'].str.contains('PACS')].index 
df['PACS_70'][idx] = 0.02*2
df['PACS_100'][idx] = 0.02*2
df['PACS_160'][idx] = 0.02*2

#SPIRE uncertainties. These are 4% (Bendo+ 2013; Griffin+ 2013)

idx = df[df['name'].str.contains('SPIRE')].index 
df['SPIRE_250'][idx] = 0.04*2
df['SPIRE_350'][idx] = 0.04*2
df['SPIRE_500'][idx] = 0.04*2

df.to_csv('../corr_uncert.csv')

print('Complete!')