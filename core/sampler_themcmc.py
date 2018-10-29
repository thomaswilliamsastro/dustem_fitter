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
from psutil import virtual_memory
import os
import sys
from collections import OrderedDict

#emcee-related imports

import emcee
from multiprocessing import Pool, cpu_count   

#For calculating redshift, given distance

import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

#THEMCMC imports

import general
from fortran_funcs import covariance_matrix,trapz

#MAIN SAMPLING FUNCTION

def sample(method,
           components,
           flux_file,
           filter_file,
           gal_row,
           pandas_dfs,
           mpi,
           overwrite):
    
    #Read in models
    
    global sCM20_df,lCM20_df,aSilM5_df
    
    #Read in the DustEM grid
    
    sCM20_df,\
        lCM20_df,\
        aSilM5_df,\
        wavelength_df = pandas_dfs
    
    #Read in the useful Pandas dataframes
    
    global flux_df,filter_df,corr_uncert_df
    
    flux_df = pd.read_csv('../'+flux_file)
    filter_df = pd.read_csv('../'+filter_file)
    corr_uncert_df = pd.read_csv('corr_uncert.csv')
    
    #Define the wavelength grid (given by the dustEM output)
    
    global wavelength
    wavelength = wavelength_df['wavelength'].values.copy()
    
    global frequency
    frequency = 3e8/(wavelength*1e-6)
    
    default_total = sCM20_df['alpha_sCM20:5.00,logU:0.00'].values.copy()+\
                    lCM20_df['alpha_sCM20:5.00,logU:0.00'].values.copy()+\
                    aSilM5_df['alpha_sCM20:5.00,logU:0.00'].values.copy()
        
    global z
    z = z_at_value(Planck15.luminosity_distance,flux_df['dist'][gal_row]*u.Mpc)
    
    #Create a dictionary of the filters
    
    global filter_dict
    filter_dict = {}
    
    for filter_name in filter_df.dtypes.index[1:]:
        
        filter_wavelength,transmission = np.loadtxt('../filters/'+filter_name+'.dat',
                                         unpack=True)
        
        filter_wavelength /= 1e4
        transmission /= np.max(transmission)
        
        filter_dict[filter_name] = filter_wavelength,transmission
        
    gal_name = flux_df['name'][gal_row]
    
    #Pull out fluxes and flux errors
    
    obs_flux = []
    obs_error = []
    obs_wavelengths = []
    
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
                            obs_wavelengths.append(filter_df[key][0])
                            
                    except:
                        
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])
                            keys.append(key)        
                            obs_wavelengths.append (filter_df[key][0])                                
                
        except KeyError:
            pass
        
    obs_flux = np.array(obs_flux)
    obs_error = np.array(obs_error)
    obs_wavelengths = np.array(obs_wavelengths)
    
    #Define stars, locking the initial scaling guess to shortest wavelength
    
    idx = np.where( obs_wavelengths == np.min(obs_wavelengths) )
    idx_key = keys[idx[0][0]]
    
    stars = general.define_stars(flux_df,
                                 gal_row,
                                 filter_df,
                                 frequency,
                                 idx_key)
    
    #Read in the pickle jar if it exists, else do the fitting

    if os.path.exists('../samples/'+gal_name+'_'+method+'_'+str(components)+'comp.h5') and not overwrite:
 
        print('Reading in '+gal_name+' samples')
     
        samples_df = pd.read_hdf('../samples/'+gal_name+'_'+method+'_'+str(components)+'comp.h5',
                                 'samples')
            
    else:
        
        #Make sure the program doesn't run into swap

        ram_footprint = sys.getsizeof(sCM20_df)*4 #approx footprint (generous!)
        mem = virtual_memory().available #free available RAM
        
        procs = cpu_count()
        
        #Run with the minimum processors that will either (a) not quite run into
        #swap, (b) maxes out the machine or (c) 4 processes (since it doesn't scale well
        #beyond that)
        
        processes = np.min([int(np.floor(mem/ram_footprint)),procs,4])
    
        print('Fitting '+gal_name+' using '+str(processes)+' processes')
        
        #Build up the matrices for the errors
        
        rms_err = np.array(np.zeros([len(obs_flux),len(obs_flux)]))
        
        for i in range(len(obs_error)):
            rms_err[i,i] = obs_error[i]**2
            
        #Uncorrelated calibration errors
            
        uncorr_err = np.array(np.zeros([len(obs_flux),len(obs_flux)]))
        
        i = 0
        
        for key in keys:
            
            uncorr_err[i,i] = (filter_df[key][1]*obs_flux[i])**2
            
            i += 1
            
        #And finally, correlated calibration errors
        
        corr_err = np.array(np.zeros([len(obs_flux),len(obs_flux)]))
        
        for i in range(len(obs_flux)):
            for j in range(i+1):
                
                corr_err[i,j] = obs_flux[i]*obs_flux[j]* \
                                corr_uncert_df[keys[j]][corr_uncert_df.index[corr_uncert_df['name'] == keys[j]][0]]* \
                                corr_uncert_df[keys[i]][corr_uncert_df.index[corr_uncert_df['name'] == keys[i]][0]]
                                
                corr_err[j,i] = corr_err[i,j]
        
        global total_err    
        total_err = rms_err+uncorr_err+corr_err
        
        pos = []
        nwalkers = 500
        
        ####DEFAULT THEMIS MIX####
        
        #Set up logU. 0 for one component, 0 and 3 for two and evenly spaced
        #between those for more
        
        log_u_selection = np.linspace(0,3,components)
        
        #Find an initial dust scaling for each component. Lock to flux closest to
        #particular wavelength peak
        
        dust_scaling = []
        
        for log_u in log_u_selection:
            
            total = sCM20_df['alpha_sCM20:5.00,logU:%.2f' % log_u].values.copy()+\
                lCM20_df['alpha_sCM20:5.00,logU:%.2f' % log_u].values.copy()+\
                aSilM5_df['alpha_sCM20:5.00,logU:%.2f' % log_u].values.copy()
                
            idx_max = np.where(total == np.max(total))[0][0]
            
            idx = np.where(np.abs(obs_wavelengths-wavelength[idx_max]) == np.min(np.abs(obs_wavelengths-wavelength[idx_max])))
    
            idx_key = keys[idx[0][0]]
        
            dust_scaling.append(flux_df[idx_key][gal_row]/total[idx_max])
        
        if method == 'default':
            
            #Set up the MCMC. We have 1+(2*components) free parameters.
            #Stellar scaling,
            #ISRF strength, 
            #overall scaling factor for the dust grains for each component.
            
            #Initial guesses: 
            
            #Since we've already normalised the stellar parameter, set this to 1. 
            
            #Log ISFR is selected from the setup above.
            
            #The scaling is based on the peak of each component's SED, from above.
                     
            ndim = 1+2*components
             
            for i in range(nwalkers):
                
                values_var = []
                
                values_var.append( np.abs(np.random.normal(loc=1,scale=1e-2)) )
                
                for component in range(components):
                    
                    values_var.append(np.random.normal(loc=log_u_selection[component],
                                                       scale=1e-2))
                    values_var.append(np.abs(np.random.normal(loc=dust_scaling[component],
                                                              scale=1e-2*dust_scaling[component])))
            
                pos.append(values_var)  
                
        ####ALLOWING VARYING ABUNDANCES####
        
        if method == 'abundfree':
            
            #Set up the MCMC. We have 1+(5*components) free parameters.
            #Stellar scaling,
            #ISRF strength, 
            #overall scaling factor for the dust grains for each component,
            #and deviation from default abundance for each grain type.
            
            #Initial guesses: 
            
            #Since we've already normalised the stellar parameter, set this to 1. 
            
            #Log ISFR is selected from the setup above.
            
            #The scaling is based on the peak of each component's SED, from above.
             
            #The deviations all default to 1
                     
            ndim = 1+5*components
             
            for i in range(nwalkers):
                
                values_var = []
                
                values_var.append( np.abs(np.random.normal(loc=1,scale=1e-2)) )
                
                for component in range(components):
                    
                    values_var.append(np.random.normal(loc=log_u_selection[component],
                                                       scale=1e-2))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=dust_scaling[component],
                                                              scale=1e-2*dust_scaling[component])))
            
                pos.append(values_var)  
    
        ####VARYING SMALL CARBON GRAIN SIZE DISTRIBUTION####
        
        if method == 'ascfree':
            
            #Set up the MCMC. We have 2+(5*components) free parameters.
            #Stellar scaling,
            #alpha_sCM20
            #ISRF strength, 
            #overall scaling factor for the dust grains for each component,
            #and deviation from default abundance for each grain type.
            
            #Initial guesses: 
            
            #Since we've already normalised the stellar parameter, set this to 1. 
            
            #Log ISFR is selected from the setup above.
            
            #The scaling is based on the peak of each component's SED, from above.
             
            #The deviations all default to 1
            
            #Set up the MCMC. We have 7 free parameters.
            #ISRF strength,
            #Stellar scaling,
            #the overall scaling factor for the dust grains.
            
            #Set the initial guesses for the slope and abundances at the default 
            #THEMIS parameters. The overall scaling is given by the ratio to 250 micron
            #earlier, and since we've already normalised the stellar parameter set this
            #to 1. The ISRF is 10^0, i.e. MW default.
                     
            ndim = 2+5*components
            
            for i in range(nwalkers):
                
                values_var = []
                
                values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                values_var.append(np.abs(np.random.normal(loc=5,scale=1e-2*5)))
                
                for component in range(components):
                    
                    values_var.append(np.random.normal(loc=log_u_selection[component],
                                                       scale=1e-2))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=1,scale=1e-2)))
                    values_var.append(np.abs(np.random.normal(loc=dust_scaling[component],
                                                              scale=1e-2*dust_scaling[component])))
            
                pos.append(values_var)
                
        #Run this MCMC. Since emcee pickles any arguments passed to it, use as few
        #as possible and rely on global variables instead!
        
        pool = Pool(processes)
        
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        lnprob, 
                                        args=(method,
                                              components,
                                              obs_flux,
                                              stars),
                                        pool=pool)
         
        #Set a number of steps for the walkers, and throw away
        #the first half as burn-in
        
        #If using MPI this gets very messy so don't use
        #tqdm
        
        nsteps = 500
        
        if mpi:
            
            for i,result in enumerate(sampler.sample(pos,
                                                     iterations=nsteps)):
                pos,probability,state = result
            
        else:
        
            for i,result in tqdm(enumerate(sampler.sample(pos,
                                                      iterations=nsteps)),
                                 total=nsteps,
                                 desc='Fitting '+gal_name):
                pos,probability,state = result
            
        pool.close()
            
        samples = sampler.chain[:, int(np.floor(nsteps/2)):, :].reshape((-1, ndim))
        
        # Convert samples to pandas dataframe and save out
        
        samples_dict = OrderedDict()
        samples_dict["$\Omega_\\ast$"] = samples[:,0]
        
        if method == 'ascfree':
            
            samples_dict["$\\alpha_\mathregular{sCM20}$"] = samples[:,1]
              
        for component in range(components):
        
            if method == 'default':
            
                samples_dict["log$_{10}$ U$_"+str(component+1)+"$"] = samples[:,2*component+1]
                samples_dict["log$_{10}$ M$_\mathregular{dust,"+str(component+1)+"}$ (M$_\odot$)"] = samples[:,2*component+2]
            
            if method == 'abundfree':
                
                samples_dict["log$_{10}$ U$_"+str(component+1)+"$"] = samples[:,5*component+1]
                samples_dict["log$_{10}$ M$_\mathregular{sCM20,"+str(component+1)+"}$"] = samples[:,5*component+2]
                samples_dict["log$_{10}$ M$_\mathregular{lCM20,"+str(component+1)+"}$"] = samples[:,5*component+3]
                samples_dict["log$_{10}$ M$_\mathregular{aSilM5,"+str(component+1)+"}$"] = samples[:,5*component+4]
                samples_dict["log$_{10}$ M$_\mathregular{dust,"+str(component+1)+"}$ (M$_\odot$)"] = samples[:,5*component+5]
            
            if method == 'ascfree':
                
                samples_dict["log$_{10}$ U$_"+str(component+1)+"$"] = samples[:,5*component+2]
                samples_dict["log$_{10}$ M$_\mathregular{sCM20,"+str(component+1)+"}$"] = samples[:,5*component+3]
                samples_dict["log$_{10}$ M$_\mathregular{lCM20,"+str(component+1)+"}$"] = samples[:,5*component+4]
                samples_dict["log$_{10}$ M$_\mathregular{aSilM5,"+str(component+1)+"}$"] = samples[:,5*component+5]
                samples_dict["log$_{10}$ M$_\mathregular{dust,"+str(component+1)+"}$ (M$_\odot$)"] = samples[:,5*component+6]
        
        samples_df = pd.DataFrame(samples_dict)
        
        samples_df.to_hdf('../samples/'+gal_name+'_'+method+'_'+str(components)+'comp.h5',
                          'samples',mode='w')
            
    return samples_df,filter_dict
        
#EMCEE-RELATED FUNCTIONS

def lnprob(theta,
           method,
           components,
           obs_flux,
           stars):
    
    lp = priors(theta,
                method,
                components)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,
                       method,
                       components,
                       obs_flux,
                       stars)
    
def priors(theta,
           method,
           components):
    
    #log_ISRF must be between -2 and 7.
    #alpha_sCM20 must be between 2 and 7.
    #Multiplicative factors must be greater than 0.
    #Redshift must be between 0 and 15.
    
    #Also insist that each component increases in ISRF
    #strength

    global z
    
    if not 0<=z<=15:
        return -np.inf
    
    omega_star = theta[0]
    
    if not omega_star>=0:
        return -np.inf
    
    if method == 'ascfree':
        
        alpha = theta[1]
        
        if not 2<=alpha<=7:
            return -np.inf
    
    for component in range(components):
    
        if method == 'default':
            
            isrf = theta[2*component+1]
            dust_scaling = theta[2*component+2]
            
            if not -2<=isrf<=7 or not dust_scaling>=0:
                return -np.inf
            
            if component > 0:
            
                if not isrf>=theta[2*(component-1)+1]:
                    return -np.inf
            
        if method == 'abundfree':
            
            isrf = theta[5*component+1]
            y_sCM20 = theta[5*component+2]
            y_lCM20 = theta[5*component+3]
            y_aSilM5 = theta[5*component+4]
            dust_scaling = theta[5*component+5]   
            
            if not -2<=isrf<=7 or not dust_scaling>=0\
                or not y_sCM20>=0 or not y_lCM20>=0\
                or not y_aSilM5>=0:
                return -np.inf 
            
            if component > 0:
            
                if not isrf>=theta[5*(component-1)+1]:
                    return -np.inf
        
        if method == 'ascfree':
        
            isrf = theta[5*component+2]
            y_sCM20 = theta[5*component+3]
            y_lCM20 = theta[5*component+4]
            y_aSilM5 = theta[5*component+5]
            dust_scaling = theta[5*component+6]   
            
            if not -2<=isrf<=7 or not dust_scaling>=0\
                or not y_sCM20>=0 or not y_lCM20>=0\
                or not y_aSilM5>=0:
                return -np.inf 
            
            if component > 0:
            
                if not isrf>=theta[5*(component-1)+2]:
                    return -np.inf
                
    return 0.0

def lnlike(theta,
           method,
           components,
           obs_flux,
           stars):

    global z
    
    omega_star = theta[0]
    
    if method == 'ascfree':
        
        alpha = theta[1]
    
    for component in range(components):
    
        if method == 'default':
            
            isrf = theta[2*component+1]
            dust_scaling = theta[2*component+2]
                
            y_sCM20 = 1
            y_lCM20 = 1
            y_aSilM5 = 1    
            alpha = 5
            
        if method == 'abundfree':
            
            isrf = theta[5*component+1]
            y_sCM20 = theta[5*component+2]
            y_lCM20 = theta[5*component+3]
            y_aSilM5 = theta[5*component+4]
            dust_scaling = theta[5*component+5]
                
            alpha = 5
        
        if method == 'ascfree':
        
            isrf = theta[5*component+2]
            y_sCM20 = theta[5*component+3]
            y_lCM20 = theta[5*component+4]
            y_aSilM5 = theta[5*component+5]
            dust_scaling = theta[5*component+6]
            
        small_grains,\
            large_grains,\
            silicates = general.read_sed(isrf,
                                         alpha,
                                         sCM20_df,
                                         lCM20_df,
                                         aSilM5_df) 
            
        if component == 0:
            
            total = dust_scaling*y_sCM20*small_grains+\
                        dust_scaling*y_lCM20*large_grains+\
                        dust_scaling*y_aSilM5*silicates
                    
        else:
            
            total += dust_scaling*y_sCM20*small_grains+\
                        dust_scaling*y_lCM20*large_grains+\
                        dust_scaling*y_aSilM5*silicates 
    
    #Include stars
    
    total += omega_star*stars
            
    filter_fluxes = filter_convolve(total,
                                    z)
    
    flux_diff = (filter_fluxes-obs_flux)[np.newaxis]
    
    chisq = covariance_matrix(flux_diff,total_err,flux_diff.T)
    
    likelihood = -0.5*chisq
    
    return likelihood

def filter_convolve(flux,
                    z):
    
    #Convolve this SED with the various filters we have to give
    #a monochromatic flux
    
    #Account for redshift
    
    wavelength_redshifted = wavelength * (1+z)
    
    filter_fluxes = [np.abs( (trapz(filter_dict[key][1]*np.interp(filter_dict[key][0],
                                                                  wavelength_redshifted,
                                                                  flux),\
                                    filter_dict[key][0])/
                                    trapz(filter_dict[key][1],filter_dict[key][0])) ) for key in keys]
    
    return np.array(filter_fluxes)
