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

os.chdir(os.getcwd())

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
    
    #For each galaxy, read in fluxes
    
    flux_df = pd.read_csv(args.fluxes+'.csv')
    filter_df = pd.read_csv(args.filters+'.csv')
    
    for gal_row in range(len(flux_df)):
          
        keys = []
        
        gal_name = flux_df['name'][gal_row]
        dist = flux_df['dist'][gal_row]
        
        #Fit the data!
            
        samples,filter_dict,keys = themcmc_sampler.sample(method=args.method,
                                                          flux_file=args.fluxes,
                                                          filter_file=args.filters,
                                                          gal_row=gal_row)
            
        if args.plot:
            
            if not os.path.isfile('plots/sed/'+gal_name+'_'+args.method+'.png'):
                
                print('Plotting '+gal_name)
            
                plotting.plot_sed(method=args.method,
                                  flux_df=flux_df,
                                  filter_df=filter_df,
                                  gal_row=gal_row,
                                  samples=samples,
                                  filter_dict=filter_dict)
                
                plotting.plot_corner(method=args.method,
                                     samples=samples,
                                     gal_name=gal_name,
                                     distance=dist)
        
        #Finally, write out code snippets for dustEM and SKIRT, if requested
        
        if args.dustemoutput:
            
            if not os.path.isfile('dustem_output/GRAIN_'+gal_name+'_'+args.method+'.dat'):
            
                print('Writing DustEM GRAIN.dat file for '+gal_name)
                
                code_snippets.dustemoutput(method=args.method,
                                           samples=samples,
                                           gal_name=gal_name)
                    
        if args.skirtoutput:
            
            if not os.path.isfile('skirt_output/template_'+gal_name+'_'+arg.smethod+'.ski'):
            
                print('Writing SKIRT code snippet '+gal_name)
                
                code_snippets.skirtoutput(method=args.method,
                                          samples=samples,
                                          gal_name=gal_name)
    
    print('Code complete, took %.2fm' % ( (time.time() - start_time)/60 ))