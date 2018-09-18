# -*- coding: utf-8 -*-
"""
THEMCMC sampler

@author: Tom Williams

v1.00.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd
import dill
from tqdm import tqdm
from scipy.constants import h,k,c
import os
from psutil import virtual_memory
import sys
from collections import OrderedDict

#emcee-related imports

import emcee
from multiprocessing import Pool

#For calculating redshift, given distance

import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

#THEMCMC imports

import general

#MAIN SAMPLING FUNCTION

def sample(method,
           flux_file,
           filter_file,
           gal_row):
    
    #Read in models
    
    global sCM20_df,lCM20_df,aSilM5_df
    
    #Read in the DustEM grid
    
    sCM20_df = pd.read_hdf('models.h5','sCM20')
    lCM20_df = pd.read_hdf('models.h5','lCM20')
    aSilM5_df = pd.read_hdf('models.h5','aSilM5')
    
    #Read in the useful Pandas dataframes
    
    global flux_df,filter_df,corr_uncert_df
    
    flux_df = pd.read_csv(flux_file+'.csv')
    filter_df = pd.read_csv(filter_file+'.csv')
    corr_uncert_df = pd.read_csv('corr_uncert.csv')
    
    #Define the wavelength grid (given by the dustEM output)
    
    global wavelength
    wavelength = sCM20_df['wavelength'].values.copy()
    
    global frequency
    frequency = 3e8/(wavelength*1e-6)
    
    default_total = sCM20_df['5.00,0.00'].values.copy()+\
                    lCM20_df['5.00,0.00'].values.copy()+\
                    aSilM5_df['5.00,0.00'].values.copy()
                    
    #If we're not fitting redshift, use the distance to calculate redshift
    
    if method != 'fit-z':
        
        global z
        z = z_at_value(Planck15.luminosity_distance,flux_df['dist'][gal_row]*u.Mpc)
    
    #Create a dictionary of the filters
    
    global filter_dict
    filter_dict = {}
    
    for filter_name in filter_df.dtypes.index[1:]:
        
        filter_wavelength,transmission = np.loadtxt('filters/'+filter_name+'.dat',
                                         unpack=True)
        
        filter_wavelength /= 1e4
        transmission /= np.max(transmission)
        
        filter_dict[filter_name] = filter_wavelength,transmission
        
    gal_name = flux_df['name'][gal_row]
    
    #Pull out fluxes and flux errors
    
    obs_flux = []
    obs_error = []
    
    global keys
    keys = []
    
    for key in filter_dict:
        
        try:
                        
            if np.isnan(flux_df[key][gal_row]) == False and \
                np.isnan(flux_df[key+'_err'][gal_row]) == False:
                
                if flux_df[key][gal_row] > 0:
                    
                    #Fit only the data with no flags
                    
                    try:
                    
                        if pd.isnull(flux_df[key+'_flag'][gal_row]):
                    
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])
                            keys.append(key)  
                            
                    except:
                        
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])
                            keys.append(key)                                         
                
        except KeyError:
            pass
        
    obs_flux = np.array(obs_flux)
    obs_error = np.array(obs_error)
    
    stars = general.define_stars(flux_df,
                                 gal_row,
                                 frequency)
    
    #Set an initial guess for the scaling variable. Lock this to the SPIRE250 flux, which
    #all galaxies must have
    
    idx = np.where(np.abs(wavelength-filter_df['SPIRE_250'][0]) == np.min(np.abs(wavelength-filter_df['SPIRE_250'][0])))
    
    initial_dust_scaling = flux_df['SPIRE_250'][gal_row]/default_total[idx[0]] * (3e8/(250*1e-6))
    
    #Read in the pickle jar if it exists, else do the fitting

    if os.path.exists('samples/'+gal_name+'_'+method+'.h5'):
 
        print('Reading in '+gal_name+' samples')
     
        samples_df = pd.read_hdf('samples/'+gal_name+'_'+method+'.h5',
                                 'samples')
            
    else:  
        
        #Make sure the program doesn't run into swap

        ram_footprint = sys.getsizeof(sCM20_df)*4 #approx footprint (generous!)
        mem = virtual_memory().available #free available RAM
        
        processes = int(np.floor(mem/ram_footprint))
        
        print('Will fit using '+str(processes)+' processes')
        
        pos = []
        nwalkers = 500
        
        ####DEFAULT THEMIS MIX####
        
        if method == 'default':
            
            #Set up the MCMC. We have 3 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #and the overall scaling factor for the dust grains.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default.
                     
            ndim = 3
             
            for i in range(nwalkers):
                
                isrf_var = np.random.normal(loc=0,scale=1e-2)
                omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                dust_scaling_var = np.abs(np.random.normal(loc=initial_dust_scaling,
                                                           scale=initial_dust_scaling*1e-2))
            
                pos.append([isrf_var,
                            omega_star_var,
                            dust_scaling_var])  
                
        ####ALLOWING VARYING ABUNDANCES####
        
        if method == 'abundfree':
            
            #Set up the MCMC. We have 6 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #deviation from default abundance of the small- and large-carbon grains
            #and silicates,
            #and the overall scaling factor for the dust grains.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default.
                     
            ndim = 6
             
            for i in range(nwalkers):
                
                isrf_var = np.random.normal(loc=0,scale=1e-2)
                omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                y_sCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                y_lCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                y_aSilM5_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                dust_scaling_var = np.abs(np.random.normal(loc=initial_dust_scaling,
                                                           scale=initial_dust_scaling*1e-2))
            
                pos.append([isrf_var,
                            omega_star_var,
                            y_sCM20_var,
                            y_lCM20_var,
                            y_aSilM5_var,
                            dust_scaling_var])     
    
        ####VARYING SMALL CARBON GRAIN SIZE DISTRIBUTION####
        
        if method == 'ascfree':
            
            #Set up the MCMC. We have 7 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #the overall scaling factor for the dust grains.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default.
                     
            ndim = 7
             
            for i in range(nwalkers):
                
                isrf_var = np.random.normal(loc=0,scale=1e-2)
                omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                alpha_var = np.abs(np.random.normal(loc=5,scale=1e-2*5))
                y_sCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                y_lCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                y_aSilM5_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                dust_scaling_var = np.abs(np.random.normal(loc=initial_dust_scaling,
                                                           scale=initial_dust_scaling*1e-2))
            
                pos.append([isrf_var,
                            omega_star_var,
                            alpha_var,
                            y_sCM20_var,
                            y_lCM20_var,
                            y_aSilM5_var,
                            dust_scaling_var])
                
        ####FIT REDSHIFT####
        
        if method == 'fit-z':
            
            #Set up the MCMC. We have 4 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #the overall scaling factor for the dust grains,
            #and the redshift of the object.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default. Set the initial redshift guess to 1.
                     
            ndim = 4
             
            for i in range(nwalkers):
                
                isrf_var = np.random.normal(loc=0,scale=1e-2)
                omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                dust_scaling_var = np.abs(np.random.normal(loc=initial_dust_scaling,
                                                           scale=initial_dust_scaling*1e-2))
                z_var = np.abs(np.random.normal(loc=1,
                                                scale=1e-2))
            
                pos.append([isrf_var,
                            omega_star_var,
                            dust_scaling_var,
                            z_var])
                
        #Run this MCMC. Since emcee pickles any arguments passed to it, use as few
        #as possible and rely on global variables instead!
        
        pool = Pool(processes)
        
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        lnprob, 
                                        args=(method,
                                              obs_flux,
                                              obs_error,
                                              stars),
                                        pool=pool)
         
        #500 steps for the 500 walkers, but throw away
        #the first 250 as a burn-in
        
        nsteps = 500
        for i,result in tqdm(enumerate(sampler.sample(pos,
                                                  iterations=nsteps)),
                             total=nsteps,
                             desc='Fitting '+gal_name):
            pos,probability,state = result
            
        pool.close()
            
        samples = sampler.chain[:, 250:, :].reshape((-1, ndim))
        
        # Convert samples to pandas dataframe and save out
        
        if method == 'default':
        
            samples_df = pd.DataFrame( OrderedDict( (("log$_{10}$ U",samples[:,0]),
                                                     ("$\Omega_\\ast$",samples[:,1]),
                                                     ("log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)",samples[:,2]) )))
            
        if method == 'abundfree':
            
            samples_df = pd.DataFrame( OrderedDict( (("log$_{10}$ U",samples[:,0]),
                                                     ("$\Omega_\\ast$",samples[:,1]),
                                                     ("log$_{10}$ M$_\mathregular{sCM20}$",samples[:,2]),
                                                     ("log$_{10}$ M$_\mathregular{lCM20}$",samples[:,3]),
                                                     ("log$_{10}$ M$_\mathregular{aSilM5}$",samples[:,4]),
                                                     ("log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)",samples[:,5]) )))
            
        if method == 'ascfree':
            
            samples_df = pd.DataFrame( OrderedDict( (("log$_{10}$ U",samples[:,0]),
                                                     ("$\Omega_\\ast$",samples[:,1]),
                                                     ("$\alpha_\mathregular{sCM20}$",samples[:,2]),
                                                     ("log$_{10}$ M$_\mathregular{sCM20}$",samples[:,3]),
                                                     ("log$_{10}$ M$_\mathregular{lCM20}$",samples[:,4]),
                                                     ("log$_{10}$ M$_\mathregular{aSilM5}$",samples[:,5]),
                                                     ("log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)",samples[:,6]) )))
            
        if method == 'fit-z':
        
            samples_df = pd.DataFrame( OrderedDict( (("log$_{10}$ U",samples[:,0]),
                                                     ("$\Omega_\\ast$",samples[:,1]),
                                                     ("log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)",samples[:,2]),
                                                     ("$z$",samples[:,3]) ))) 
        
        samples_df.to_hdf('samples/'+gal_name+'_'+method+'.h5',
                          'samples',mode='w')
            
    return samples_df,filter_dict,keys
        
#EMCEE-RELATED FUNCTIONS

def lnlike(theta,
           method,
           obs_flux,
           obs_error,
           stars):


    global z
    
    if method == 'default':
        
        isrf,\
            omega_star,\
            dust_scaling = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1    
        
    if method == 'abundfree':
        
        isrf,\
            omega_star,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta     
            
        alpha = 5  
    
    if method == 'ascfree':
    
        isrf,\
            omega_star,\
            alpha,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta
            
    if method == 'fit-z':
        
        isrf,\
            omega_star,\
            dust_scaling,\
            z = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1
        
    small_grains,\
        large_grains,\
        silicates = general.read_sed(isrf,
                                     alpha,
                                     sCM20_df,
                                     lCM20_df,
                                     aSilM5_df,
                                     frequency)    
    
    #Scale everything accordingly
    
    total = dust_scaling*y_sCM20*small_grains+\
            dust_scaling*y_lCM20*large_grains+\
            dust_scaling*y_aSilM5*silicates+omega_star*stars
            
    filter_fluxes = filter_convolve(total,
                                    z)
    
    #Build up a matrix for the various uncertainties
    
    rms_err = np.matrix(np.zeros([len(obs_flux),len(obs_flux)]))
    
    for i in range(len(obs_error)):
        rms_err[i,i] = obs_error[i]**2
        
    #Uncorrelated calibration errors
        
    uncorr_err = np.matrix(np.zeros([len(obs_flux),len(obs_flux)]))
    
    i = 0
    
    for key in keys:
        
        uncorr_err[i,i] = (filter_df[key][1]*obs_flux[i])**2
        
        i += 1
        
    #And finally, correlated calibration errors
    
    corr_err = np.matrix(np.zeros([len(obs_flux),len(obs_flux)]))
    
    for i in range(len(obs_flux)):
        for j in range(len(obs_flux)):
            
            corr_err[i,j] = obs_flux[i]*obs_flux[j]* \
                            corr_uncert_df[keys[j]][corr_uncert_df.index[corr_uncert_df['name'] == keys[j]][0]]* \
                            corr_uncert_df[keys[i]][corr_uncert_df.index[corr_uncert_df['name'] == keys[i]][0]]
    
    total_err = rms_err+uncorr_err+corr_err
    
    flux_diff = (filter_fluxes-obs_flux)[np.newaxis]
    
    chisq = flux_diff*total_err.I*flux_diff.T
    
    likelihood = -0.5*chisq
    
    return likelihood

def lnprob(theta,
           method,
           obs_flux,
           obs_error,
           stars):
    
    lp = priors(theta,
                method)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,
                       method,
                       obs_flux,
                       obs_error,
                       stars)

def priors(theta,
           method):
    
    #log_ISRF must be between -1 and 3.5.
    #alpha_sCM20 must be between 2.6 and 5.4.
    #Multiplicative factors must be greater than 0.
    #Redshift must be between 0 and 15

    global z
    
    if method == 'default':
        
        isrf,\
            omega_star,\
            dust_scaling = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1        
        
    if method == 'abundfree':
        
        isrf,\
            omega_star,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta    
            
        alpha = 5      
    
    if method == 'ascfree':
    
        isrf,\
            omega_star,\
            alpha,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta
            
    if method == 'fit-z':
        
        isrf,\
            omega_star,\
            dust_scaling,\
            z = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1
    
    if -1 <= isrf <= 3.5 and \
        omega_star > 0 and \
        2.6<=alpha<=5.4 and \
        y_sCM20 > 0 and \
        y_lCM20 > 0 and \
        y_aSilM5 > 0 and \
        dust_scaling > 0 and \
        0<=z<=15:
        return 0.0
    else:
        return -np.inf

def filter_convolve(flux,
                    z):
    
    #Convolve this SED with the various filters we have to give
    #a monochromatic flux
 
    filter_fluxes = []
    
    #Account for redshift
    
    wavelength_redshifted = wavelength/(1+z)
    
    for key in keys:
                 
        #Convolve MBB with filter
        
        filter_flux = np.interp(filter_dict[key][0],wavelength_redshifted,flux)
                 
        filter_fluxes.append( np.abs( (np.trapz(filter_dict[key][1]*filter_flux,filter_dict[key][0])/
                                      np.trapz(filter_dict[key][1],filter_dict[key][0])) ) )
    
    return np.array(filter_fluxes)
