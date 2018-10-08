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

os.chdir(os.getcwd())
sys.path.append(os.getcwd())

#THEMCMC imports

import sampler_themcmc
import plotting
import general
import code_snippets

try:
    from schwimmbad import MPIPool
except:
    pass

#Set up the argument parser

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='THEMCMC optional settings.')
parser.add_argument('--method',type=str,default='default',metavar='',
                    help="Method for fitting the data. Options are 'default', 'abundfree', 'ascfree'")
parser.add_argument('--plot',action='store_true',default=False,
                    help="Plot SED and corner plot.")
parser.add_argument('--units',type=str,default='flux',
                    help="Units for the plot. Options are flux (Jy) or luminosity (Lsun)")
parser.add_argument('--skirtoutput',action='store_true',default=False,
                    help="Write out SKIRT mix code snippet. (N.B. does not work with multiu!)")
parser.add_argument('--dustemoutput',action='store_true',default=False,
                    help="Write out DustEM GRAIN.dat file (N.B. does not work with multiu!).")
parser.add_argument('--fluxes',type=str,default='fluxes',metavar='',
                    help="File containing 'fluxes' to fit.")
parser.add_argument('--mpi',action='store_true',default=False,
                    help="Run with MPI (requires Schwimmbad, and --bind-to none).")

args = parser.parse_args()

#Create folders for sample pickle jars and plots, if they don't exits

os.chdir(os.getcwd())

if not os.path.exists('../plots') and args.plot:
    os.mkdir('../plots')
if not os.path.exists('../plots/sed') and args.plot:
    os.mkdir('../plots/sed')
if not os.path.exists('../plots/corner') and args.plot:
    os.mkdir('../plots/corner')
if not os.path.exists('../samples'):
    os.mkdir('../samples')
if not os.path.exists('../skirt_output') and args.skirtoutput:
    os.mkdir('../skirt_output')
if not os.path.exists('../dustem_output') and args.dustemoutput:
    os.mkdir('../dustem_output')
    
def main(gal_row):
        
    gal_name = flux_df['name'][gal_row]
    
    try:
    
        dist = flux_df['dist'][gal_row]
        
    except KeyError:
                      
        raise Exception('No distance found!')
        
    samples_df,filter_dict = sampler_themcmc.sample(method=args.method,
                                                         flux_file=args.fluxes,
                                                         filter_file='filters.csv',
                                                         gal_row=gal_row,
                                                         pandas_dfs=pandas_dfs,
                                                         mpi=args.mpi)
        
    if args.plot:
        
        if not os.path.isfile('../plots/sed/'+gal_name+'_'+args.method+'.png'):
            
            print('Plotting SED')
        
            plotting.plot_sed(method=args.method,
                              flux_df=flux_df,
                              filter_df=filter_df,
                              pandas_dfs=pandas_dfs,
                              gal_row=gal_row,
                              samples_df=samples_df,
                              filter_dict=filter_dict,
                              units=args.units,
                              distance=dist)
            
        if not os.path.isfile('../plots/corner/'+gal_name+'_'+args.method+'.png'):
            
            print('Plotting corner')
            
            plotting.plot_corner(method=args.method,
                                 samples_df=samples_df,
                                 gal_name=gal_name,
                                 distance=dist)
    
    #Finally, write out code snippets for dustEM and SKIRT, if requested
    
    if args.dustemoutput:
        
        if not os.path.isfile('../dustem_output/GRAIN_'+gal_name+'_'+args.method+'.dat'):
        
            print('Writing DustEM GRAIN.dat file')
            
            code_snippets.dustemoutput(method=args.method,
                                       samples_df=samples_df,
                                       gal_name=gal_name)
                
    if args.skirtoutput:
        
        if not os.path.isfile('../skirt_output/template_'+gal_name+'_'+args.method+'.ski'):
        
            print('Writing SKIRT code snippet')
            
            code_snippets.skirtoutput(method=args.method,
                                      samples_df=samples_df,
                                      gal_name=gal_name)

if __name__ == "__main__":

    start_time = time.time()
    
    #For each galaxy, read in fluxes
    
    flux_df = pd.read_csv('../'+args.fluxes)
    filter_df = pd.read_csv('../filters.csv')
    
    #Read in and zip up dataframes
    
    sCM20_df = pd.read_hdf('models.h5','sCM20')
    lCM20_df = pd.read_hdf('models.h5','lCM20')
    aSilM5_df = pd.read_hdf('models.h5','aSilM5')
    wavelength_df = pd.read_hdf('models.h5','wavelength')
    
    pandas_dfs = [sCM20_df,lCM20_df,
                  aSilM5_df,wavelength_df]
    
    if args.mpi:
        
        mpi_pool = MPIPool()
        
        if not mpi_pool.is_master():
            mpi_pool.wait()
            sys.exit(0)
        
        mpi_pool.map( main,
                      range(len(flux_df)) )
        
        mpi_pool.close()
        
    else:
        
        for gal_row in range(len(flux_df)):
             
            main(gal_row)
    
    print('Code complete, took %.2fm' % ( (time.time() - start_time)/60 ))
