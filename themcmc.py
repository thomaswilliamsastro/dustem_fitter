# -*- coding: utf-8 -*-
"""
General purpose THEMIS MCMC fitter

@author: Tom Williams

v1.00: 20180817.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd
import time

#Argument parsing
import argparse

#OS I/O stuff
import os
import sys

sys.path.append(os.getcwd())

#THEMCMC imports

import themcmc_sampler
import plotting
import general
import code_snippets

#Set up the argument parser

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='THEMCMC optional settings.')
parser.add_argument('--method',type=str,default='default',metavar='',
                    help="Method for fitting the data. Options are 'default', 'abundfree', 'ascfree'")
parser.add_argument('--plot',action='store_true',default=False,
                    help="Plot SED and corner plot.")
parser.add_argument('--skirtoutput',action='store_true',default=False,
                    help="Write out SKIRT mix code snippet. (N.B. does not work with multiu!)")
parser.add_argument('--dustemoutput',action='store_true',default=False,
                    help="Write out DustEM GRAIN.dat file (N.B. does not work with multiu!).")
parser.add_argument('--fluxes',type=str,default='fluxes',metavar='',
                    help="File containing 'fluxes'.csv to fit.")
parser.add_argument('--filters',type=str,default='filters',metavar='',
                    help="File containing 'filters'.csv to include.")

args = parser.parse_args()

#Create folders for sample pickle jars and plots, if they don't exits

if not os.path.exists('plots'):
    os.mkdir('plots')
if not os.path.exists('plots/sed'):
    os.mkdir('plots/sed')
if not os.path.exists('plots/corner'):
    os.mkdir('plots/corner')
if not os.path.exists('samples'):
    os.mkdir('samples')
if not os.path.exists('skirt_output'):
    os.mkdir('skirt_output')
if not os.path.exists('dustem_output'):
    os.mkdir('dustem_output')

if __name__ == "__main__":

    start_time = time.time()
    os.chdir(os.getcwd())
    
    #For each galaxy, read in fluxes and see if we have sufficient points to fit
    
    flux_df = pd.read_csv(args.fluxes+'.csv')
    filter_df = pd.read_csv(args.filters+'.csv')
    
    keys = []
    
    for gal_row in range(len(flux_df)):
        
        gal_name = flux_df['name'][gal_row]
        dist = flux_df['dist_best'][gal_row]
        
        for filter_name in filter_df.dtypes.index[1:]:
            
            try:
                            
                if np.isnan(flux_df[filter_name][gal_row]) == False:
                    
                    if flux_df[filter_name][gal_row] > 0:
                        
                        #Fit only the data with no flags
                        
                        try:
                        
                            if pd.isnull(flux_df[filter_name+'_flag'][gal_row]):
                                keys.append(filter_name) 
                                
                        except KeyError:
                            keys.append(filter_name)                    
                    
            except KeyError:
                pass
            
        #Only fit if we cover the entire SED adequately. This means either 3.4
        #or 3.6, Spizer 8 or WISE 12 or IRAS 12, MIPS 24 or WISE 22 or IRAS 25,
        #MIPS 70 or PACS 70 or IRAS 60, PACS 100 or IRAS 100, MIPS 160 or PACS 160,
        #SPIRE 250, SPIRE 350 or Planck 350, SPIRE 500 or Planck 550.
        
        if 'Spitzer_3.6' not in keys and 'WISE_3.4' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'Spitzer_8.0' not in keys and 'WISE_12' not in keys and 'IRAS_12' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'Spitzer_24' not in keys and 'WISE_22' not in keys and 'IRAS_25' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'Spitzer_70' not in keys and 'PACS_70' not in keys and 'IRAS_60' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'PACS_100' not in keys and 'IRAS_100' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'Spitzer_160' not in keys and 'PACS_160' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'SPIRE_250' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'SPIRE_350' not in keys and 'Planck_350' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        if 'SPIRE_500' not in keys and 'Planck_550' not in keys:
            print('Insufficient points to fit '+gal_name)
            continue
        
        #Else, we have sufficient coverage so fit!
            
        samples,filter_dict,keys = themcmc_sampler.sample(method=args.method,
                                                          flux_file=args.fluxes,
                                                          filter_file=args.filters,
                                                          gal_row=gal_row)
            
        if args.plot:
            
            print('Plotting')
            
            plotting.plot_sed(method=args.method,
                              flux_df=flux_df,
                              filter_df=filter_df,
                              gal_row=gal_row,
                              samples=samples,
                              filter_dict=filter_dict,
                              keys=keys)
            
            plotting.plot_corner(method=args.method,
                                 samples=samples,
                                 gal_name=gal_name,
                                 distance=dist)
        
        #Finally, write out code snippets for dustEM and SKIRT, if requested
        
        if args.dustemoutput:
            
            print('Writing DustEM GRAIN.dat file')
            
            code_snippets.dustemoutput(method=args.method,
                                       samples=samples,
                                       gal_name=gal_name)
                    
        if args.skirtoutput:
            
            print('Writing SKIRT code snippet')
            
            code_snippets.skirtoutput(method=args.method,
                                      samples=samples,
                                      gal_name=gal_name)
    
    print('Code complete, took %.2fm' % ( (time.time() - start_time)/60 ))