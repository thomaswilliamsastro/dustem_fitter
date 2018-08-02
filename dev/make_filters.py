# -*- coding: utf-8 -*-
"""
Create a Pandas dataframe for the filters

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import pandas as pd

d = {'WISE3_4':[3.368],
     'IRAC3_6':[3.6],
     'IRAC4_5':[4.5],
     'WISE4_6':[4.618],
     'IRAC5_8':[5.8],
     'IRAC8':[8.0],
     'WISE12':[12.082],
     'WISE22':[22.194],
     'MIPS24':[24.0],
     'MIPS70':[70.0],
     'PACS100':[100.0],
     'PACS160':[160.0],
     'SPIRE250':[250.0],
     'SPIRE350':[350.0],
     'SPIRE500':[500.0],
     'planck353':[850.0],
     'SCUBA2_450':[450.0],
     'SCUBA2_850':[850]}

df = pd.DataFrame(data=d)

df.to_csv('/home/daedalusdata/c1625914/dustem_fitter/filters.csv')

print('Complete!')