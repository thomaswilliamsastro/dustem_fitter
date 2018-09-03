# -*- coding: utf-8 -*-
"""
THEMCMC sampler

@author: Tom Williams

v0.00.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import numpy as np
import pandas as pd

#emcee-related imports

import emcee
from multiprocessing import Pool

from tqdm import tqdm
from scipy.constants import h,k,c
import os

#THEMCMC imports

import general

#MAIN SAMPLING FUNCTION

def sample(method,
           flux_file,
           filter_file,
           gal_row):
    
    #Read in the DustEM grid
    
    global sCM20_df
    sCM20_df = pd.read_hdf('models.h5','sCM20')
    global lCM20_df
    lCM20_df = pd.read_hdf('models.h5','lCM20')
    global aSilM5_df
    aSilM5_df = pd.read_hdf('models.h5','aSilM5')
    
    #Read in the useful Pandas dataframes
    
    global flux_df
    flux_df = pd.read_csv(flux_file+'.csv')
    global filter_df
    filter_df = pd.read_csv(filter_file+'.csv')
    global corr_uncert_df
    corr_uncert_df = pd.read_csv('corr_uncert.csv')
    
    #Define the wavelength grid (given by the dustEM output)
    
    global wavelength
    wavelength = sCM20_df['wavelength'].values.copy()
    
    global frequency
    frequency = 3e8/(wavelength*1e-6)
    
    default_total = sCM20_df['5.00,0.00'].values.copy()+\
                    lCM20_df['5.00,0.00'].values.copy()+\
                    aSilM5_df['5.00,0.00'].values.copy()
    
    #Create a dictionary of the filters
    
    global filter_dict
    filter_dict = {}
    
    for filter_name in filter_df.dtypes.index[1:]:
        
        filter_wavelength,transmission = np.loadtxt('filters/'+filter_name+'.dat',
                                         unpack=True)
        
        filter_wavelength /= 1e4
        transmission /= np.max(transmission)
        
        filter_dict[filter_name] = filter_wavelength,transmission
        
    distance = flux_df['dist_best'][gal_row]
    gal_name = flux_df['name'][gal_row]
    
    #Pull out fluxes and flux errors
    
    obs_flux = []
    obs_error = []
    
    global keys
    keys = []
    
    for key in filter_dict:
        
        try:
                        
            if np.isnan(flux_df[key][gal_row]) == False:
                
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
    
    #Create a blackbody of 5000K to represent the stars, and lock this
    #at the 3.6micron IRAC flux (or the 3.4 WISE flux if not available)
    
    stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1
    
    try:
        if not np.isnan(flux_df['Spitzer_3.6'][gal_row]):
            idx = np.where(np.abs(wavelength-3.6) == np.min(np.abs(wavelength-3.6)))
            ratio = flux_df['Spitzer_3.6'][gal_row]/stars[idx]
        else:
            idx = np.where(np.abs(wavelength-3.4) == np.min(np.abs(wavelength-3.4)))
            ratio = flux_df['WISE_3.4'][gal_row]/stars[idx]
    except KeyError:
        idx = np.where(np.abs(wavelength-3.4) == np.min(np.abs(wavelength-3.4)))
        ratio = flux_df['WISE_3.4'][gal_row]/stars[idx]
    
    stars *= ratio
    
    #Set an initial guess for the scaling variable. Lock this by a default THEMIS SED
    #to the 250 micron flux
    
    idx = np.where(np.abs(wavelength-250) == np.min(np.abs(wavelength-250)))
    initial_dust_scaling = flux_df['SPIRE_250'][gal_row]/default_total[idx] * (3e8/250e-6)
    
    #Read in the pickle jar if it exists, else do the fitting
    
    if os.path.exists('samples/'+gal_name+'_samples_scfree.hkl'):
        
        print('Reading in '+gal_name+' pickle jar')
        
        with open('samples/'+gal_name+'_samples_'+args.method+'.hkl', 'rb') as samples_dj:
            samples = dill.load(samples_dj)    
            
    else:
        
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
                     
            ndim,nwalkers = (3,500)
             
            pos = []
             
            for i in range(nwalkers):
                
                isrf_var = np.random.normal(loc=0,scale=1e-2)
                omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
                dust_scaling_var = np.abs(np.random.normal(loc=initial_dust_scaling,
                                                           scale=initial_dust_scaling*1e-2))
            
                pos.append([isrf_var,
                            omega_star_var,
                            dust_scaling_var])           

        ####VARYING SMALL CARBON GRAIN SIZE DISTRIBUTION####
        
        if method == 'ascfree':
            
            #Set up the MCMC. We have 7 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #Power-law slope for small carbon grains,
            #deviation from default abundance of the small- and large-carbon grains
            #and silicates,
            #and the overall scaling factor for the dust grains.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default.
                     
            ndim,nwalkers = (7,500)
             
            pos = []
             
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
                
        #Run this MCMC. Since emcee pickles any arguments passed to it, use as few
        #as possible and rely on global variables instead!
        
        pool = Pool()
        
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
        
        # Save samples to dill pickle jar
        with open('samples/'+gal_name+'_samples_'+method+'.hkl', 'wb') as samples_dj:
            dill.dump(samples, samples_dj)
            
        return samples
        
#EMCEE-RELATED FUNCTIONS

def lnlike(theta,
           method,
           obs_flux,
           obs_error,
           stars):
    
    if method == 'default':
        
        isrf,\
            omega_star,\
            dust_scaling = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1      
    
    if method == 'ascfree':
    
        isrf,\
            omega_star,\
            alpha,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta
        
    small_grains,\
        large_grains,\
        silicates = general.read_sed(isrf,
                                     alpha)    
    
    #Scale everything accordingly
    
    total = dust_scaling*y_sCM20*small_grains+\
            dust_scaling*y_lCM20*large_grains+\
            dust_scaling*y_aSilM5*silicates+omega_star*stars
            
    filter_fluxes = general.filter_convolve(total)
    
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
    #Multiplicative factors must be greater than 0
    
    if method == 'default':
        
        isrf,\
            omega_star,\
            dust_scaling = theta
            
        alpha = 5
        y_sCM20 = 1
        y_lCM20 = 1
        y_aSilM5 = 1          
    
    if method == 'ascfree':
    
        isrf,\
            omega_star,\
            alpha,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = theta
    
    if -1 <= isrf <= 3.5 and \
        omega_star > 0 and \
        2.6<=alpha<=5.4 and \
        y_sCM20 > 0 and \
        y_lCM20 > 0 and \
        y_aSilM5 > 0 and \
        dust_scaling > 0:
        return 0.0
    else:
        return -np.inf