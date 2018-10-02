# -*- coding: utf-8 -*-
"""
Plot commands for THEMCMC

@author: Tom Williams

v1.00.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

#For plotting
import matplotlib
matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'cm'
import corner

import numpy as np
import pandas as pd

#For calculating redshift, given distance

import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

#THEMCMC imports

import general
from fortran_funcs import trapz

def plot_sed(method,
             flux_df,
             filter_df,
             pandas_dfs,
             gal_row,
             samples_df,
             filter_dict,
             units,
             distance):
    
    gal_name = flux_df['name'][gal_row]
    
    #Read in the DustEM grid
    
    sCM20_df,\
        lCM20_df,\
        aSilM5_df,\
        wavelength_df = pandas_dfs
    
    wavelength = wavelength_df['wavelength'].values.copy()
    
    #Take redshift into account
    
    if method != 'fit-z':
    
        z = z_at_value(Planck15.luminosity_distance,flux_df['dist'][gal_row]*u.Mpc)    
        wavelength_redshifted = wavelength * (1+z)
        
    frequency = 3e8/(wavelength*1e-6)
    
    #Convert the samples dataframe back into an array for plotting
    
    samples = np.zeros(samples_df.shape)
    
    i = 0
    
    for col_name in samples_df.dtypes.index:
        
        col_values = samples_df[col_name].tolist()
        
        samples[:,i] = col_values
        i += 1
    
    #Pull out fluxes
    
    obs_flux = []
    obs_error = []
    obs_wavelength = []
    obs_flag = []    
    keys = []
    
    for key in filter_dict:
        
        try:
                        
            if np.isnan(flux_df[key][gal_row]) == False:
                
                if flux_df[key][gal_row] > 0:
                    
                    #Fit only the data with no flags
                    
                    try:
                            
                        if pd.isnull(flux_df[key+'_flag'][gal_row]):       
                            obs_wavelength.append(filter_df[key][0])
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])                     
                            obs_flag.append(0)  
                                                    
                        else:       
                            obs_wavelength.append(filter_df[key][0])
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])                     
                            obs_flag.append(1)
                            
                        keys.append(key) 
                            
                    except:
                            obs_wavelength.append(filter_df[key][0])
                            obs_flux.append(flux_df[key][gal_row])
                            obs_error.append(flux_df[key+'_err'][gal_row])    
                            obs_flag.append(0)    
                            keys.append(key)                              
                
        except KeyError:
            pass
        
    obs_flux = np.array(obs_flux)
    obs_error = np.array(obs_error)
    obs_wavelength = np.array(obs_wavelength)
    obs_flag = np.array(obs_flag)
    
    idx = np.where( obs_wavelength[obs_flag == 0] == np.min(obs_wavelength[obs_flag == 0]) )
    
    filtered_keys = []
    
    for i in range(len(obs_flag)):
        
        if not obs_flag[i]:
            filtered_keys.append(keys[i])
    
    idx_key = filtered_keys[idx[0][0]]
    
    #Generate stars
    
    stars = general.define_stars(flux_df,
                                 gal_row,
                                 filter_df,
                                 frequency,
                                 idx_key)

    #For errorbars
    
    y_to_percentile_stars = []
    y_to_percentile_small = []
    y_to_percentile_large = []
    y_to_percentile_silicates = []
    y_to_percentile_total = []
    
    for i in range(150):
        
        if method == 'default':
            
            isrf,\
                omega_star,\
                dust_scaling = samples[np.random.randint(len(samples)),:]
             
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
                dust_scaling = samples[np.random.randint(len(samples)),:]    
                
            alpha = 5       
                
        if method == 'ascfree':
            
            isrf,\
                omega_star,\
                alpha,\
                y_sCM20,\
                y_lCM20,\
                y_aSilM5,\
                dust_scaling = samples[np.random.randint(len(samples)),:]
                
        if method == 'fit-z':
        
            isrf,\
                omega_star,\
                dust_scaling,\
                z = samples[np.random.randint(len(samples)),:]
                
            alpha = 5
            y_sCM20 = 1
            y_lCM20 = 1
            y_aSilM5 = 1
            
            wavelength_redshifted = wavelength * (1+z)
               
        small_grains,\
            large_grains,\
            silicates = general.read_sed(isrf,
                                         alpha,
                                         sCM20_df,
                                         lCM20_df,
                                         aSilM5_df)
        
        y = omega_star*stars   
        y_to_percentile_stars.append(y)
        
        y = y_sCM20*small_grains*dust_scaling   
        y_to_percentile_small.append(y)
        
        y = y_lCM20*large_grains*dust_scaling 
        y_to_percentile_large.append(y)
        
        y = y_aSilM5*silicates*dust_scaling 
        y_to_percentile_silicates.append(y)
        
        y = dust_scaling*(y_sCM20*small_grains+\
                           y_lCM20*large_grains+\
                           y_aSilM5*silicates)+\
                           omega_star*stars
        y_to_percentile_total.append(y)
        
        
    y_upper_stars = np.percentile(y_to_percentile_stars,84,
                                  axis=0)
    y_lower_stars = np.percentile(y_to_percentile_stars,16,
                                  axis=0)
    y_upper_small = np.percentile(y_to_percentile_small,84,
                                  axis=0)
    y_lower_small = np.percentile(y_to_percentile_small,16,
                                  axis=0)
    y_upper_large = np.percentile(y_to_percentile_large,84,
                                  axis=0)
    y_lower_large = np.percentile(y_to_percentile_large,16,
                                  axis=0)
    y_upper_silicates = np.percentile(y_to_percentile_silicates,84,
                                  axis=0)
    y_lower_silicates = np.percentile(y_to_percentile_silicates,16,
                                  axis=0)
    y_upper = np.percentile(y_to_percentile_total,84,
                            axis=0)
    y_lower = np.percentile(y_to_percentile_total,16,
                            axis=0)
    
    #And calculate the median lines
    
    y_median_stars = np.percentile(y_to_percentile_stars,50,
                                   axis=0)
    y_median_small = np.percentile(y_to_percentile_small,50,
                                   axis=0)
    y_median_large = np.percentile(y_to_percentile_large,50,
                                   axis=0)
    y_median_silicates = np.percentile(y_to_percentile_silicates,50,
                                       axis=0)
    y_median_total = np.percentile(y_to_percentile_total,50,
                                   axis=0)
    
    #If outputting luminosity, convert all these fluxes accordingly
    
    if units in ['luminosity']:
        
        y_upper_stars = general.convert_to_luminosity(y_upper_stars,
                                                      distance,
                                                      frequency)
        y_lower_stars = general.convert_to_luminosity(y_lower_stars,
                                                      distance,
                                                      frequency)
        y_upper_small = general.convert_to_luminosity(y_upper_small,
                                                      distance,
                                                      frequency)
        y_lower_small = general.convert_to_luminosity(y_lower_small,
                                                      distance,
                                                      frequency)
        y_upper_large = general.convert_to_luminosity(y_upper_large,
                                                      distance,
                                                      frequency)
        y_lower_large = general.convert_to_luminosity(y_lower_large,
                                                      distance,
                                                      frequency)
        y_upper_silicates = general.convert_to_luminosity(y_upper_silicates,
                                                          distance,
                                                          frequency)
        y_lower_silicates = general.convert_to_luminosity(y_lower_silicates,
                                                          distance,
                                                          frequency)
        y_upper = general.convert_to_luminosity(y_upper,
                                                distance,
                                                frequency)
        y_lower = general.convert_to_luminosity(y_upper,
                                                distance,
                                                frequency)
    
        y_median_stars = general.convert_to_luminosity(y_median_stars,
                                                       distance,
                                                       frequency)
        y_median_small = general.convert_to_luminosity(y_median_small,
                                                       distance,
                                                       frequency)
        y_median_large = general.convert_to_luminosity(y_median_large,
                                                       distance,
                                                       frequency)
        y_median_silicates = general.convert_to_luminosity(y_median_silicates,
                                                           distance,
                                                           frequency)
        y_median_total = general.convert_to_luminosity(y_median_total,
                                                       distance,
                                                       frequency)
        
        #And the actual fluxes!
        
        obs_flux = general.convert_to_luminosity(obs_flux,
                                                 distance,
                                                 3e8/(obs_wavelength*1e-6))
        obs_error = general.convert_to_luminosity(obs_error,
                                                  distance,
                                                  3e8/(obs_wavelength*1e-6))
                
    #Calculate residuals
                   
    flux_model = filter_convolve(y_median_total,
                                 wavelength,
                                 filter_dict,
                                 keys)
    
    residuals = (obs_flux-flux_model)/obs_flux
    residual_err = obs_error/obs_flux
        
    residuals = np.array(residuals)*100
    residual_err = np.array(residual_err)*100
    
    residual_upper = (y_upper-y_median_total)*100/y_median_total
    residual_lower = (y_lower-y_median_total)*100/y_median_total
    
    fig1 = plt.figure(figsize=(10,6))
    frame1 = fig1.add_axes((.1,.3,.8,.6))
    
    #Plot the best fit and errorbars. The flux errors here are only
    #RMS, so include overall calibration from Clark et al. (2018).
    
    for i in range(len(keys)):
        
        calib_uncert = {'Spitzer_3.6':0.03,
                        'Spitzer_4.5':0.03,
                        'Spitzer_5.8':0.03,
                        'Spitzer_8.0':0.03,
                        'Spitzer_24':0.05,
                        'Spitzer_70':0.1,
                        'Spitzer_160':0.12,
                        'WISE_3.4':0.029,
                        'WISE_4.6':0.034,     
                        'WISE_12':0.046,
                        'WISE_22':0.056,
                        'PACS_70':0.07,
                        'PACS_100':0.07,
                        'PACS_160':0.07,
                        'SPIRE_250':0.055,
                        'SPIRE_350':0.055,
                        'SPIRE_500':0.055,
                        'Planck_350':0.064,
                        'Planck_550':0.061,
                        'Planck_850':0.0078,
                        'SCUBA2_450':0.12,
                        'SCUBA2_850':0.08,
                        'IRAS_12':0.2,
                        'IRAS_25':0.2,
                        'IRAS_60':0.2,
                        'IRAS_100':0.2}[keys[i]]
        
        obs_error[i] += calib_uncert*obs_flux[i]
    
    #THEMIS models
    
    plt.fill_between(wavelength,y_lower_stars,y_upper_stars,
                     facecolor='m', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.3)
    plt.fill_between(wavelength,y_lower_small,y_upper_small,
                     facecolor='b', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.3)
    plt.fill_between(wavelength,y_lower_large,y_upper_large,
                     facecolor='g', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.3)
    plt.fill_between(wavelength,y_lower_silicates,y_upper_silicates,
                     facecolor='r', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.3)
    plt.fill_between(wavelength,y_lower,y_upper,
                     facecolor='k', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.4)
    
    plt.plot(wavelength,y_median_stars,
             c='m',
             ls='--',
             label='Stars')
    plt.plot(wavelength,y_median_small,
             c='b',
             ls='-.',
             label='sCM20')
    plt.plot(wavelength,y_median_large,
             c='g',
             dashes=[2,2,2,2],
             label='lCM20')
    plt.plot(wavelength,y_median_silicates,
             c='r',
             dashes=[5,2,10,2],
             label='aSilM5')
    plt.plot(wavelength,y_median_total,
             c='k',
             label='Total')
    
    #Observed fluxes
    
    plt.errorbar(obs_wavelength[obs_flag == 0],
                 obs_flux[obs_flag == 0],
                 yerr=obs_error[obs_flag == 0],
                 c='r',
                 marker='o',
                 markersize=4,
                 ls='none',
                 zorder=99)
    
    plt.errorbar(obs_wavelength[obs_flag == 1],
                 obs_flux[obs_flag == 1],
                 yerr=obs_error[obs_flag == 1],
                 c='r',
                 mfc='white',
                 marker='o',
                 markersize=4,
                 ls='none',
                 zorder=98)
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlim([1,1000])
    plt.ylim([0.5*10**np.floor(np.log10(np.min(obs_flux[obs_flag == 0]))-1),
              10**np.ceil(np.log10(np.max(obs_flux[obs_flag == 0]))+1)])
    
    #Move the legend outside of the plot so it doesn't overlap with anything
    
    plt.subplots_adjust(left=0.1,right = 0.75)
    
    plt.legend(loc=2,
               fontsize=14,
               frameon=False,
               bbox_to_anchor=(1.01, 0.5))
    
    if units in ['flux']:
    
        plt.ylabel(r'$F_\nu$ (Jy)',
                   fontsize=14)
        
    elif units in ['luminosity']:
        
        plt.ylabel(r'$\lambda L_\lambda$ ($L_\odot$)',
                   fontsize=14)
        
    else:
        
        print('Unknown unit type specified! Defaulting to Jy')
        plt.ylabel(r'$F_\nu$ (Jy)',
                   fontsize=14)
    
    plt.yticks(fontsize=14)
    
    plt.tick_params(labelbottom=False)
    
    #Add in residuals
    
    frame2=fig1.add_axes((.1,.1,.8,.2))
    
    #Include the one-sigma errors in the residuals
    
    plt.fill_between(wavelength,residual_lower,residual_upper,
                     facecolor='k', interpolate=True,lw=0.5,
                     edgecolor='none', alpha=0.4)
    
    plt.errorbar(obs_wavelength[obs_flag == 0],
                 residuals[obs_flag == 0],
                 yerr=residual_err[obs_flag == 0],
                 c='r',
                 marker='o',
                 markersize=4,
                 ls='none',
                 zorder=99)
    
    plt.errorbar(obs_wavelength[obs_flag == 1],
                 residuals[obs_flag == 1],
                 yerr=residual_err[obs_flag == 1],
                 c='r',
                 mfc='white',
                 marker='o',
                 markersize=4,
                 ls='none',
                 zorder=98)
    
    plt.axhline(0,ls='--',c='k')
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.xscale('log')
    
    plt.xlabel(r'$\lambda$ ($\mu$m)',
               fontsize=14)
    plt.ylabel('Residual (%)',
               fontsize=14)
    
    plt.xlim([1,1000])
    plt.ylim([-100,100])
    
    plt.savefig('../plots/sed/'+gal_name+'_'+method+'.png',
                bbox_inches='tight',
                dpi=150)
    plt.savefig('../plots/sed/'+gal_name+'_'+method+'.pdf',
                bbox_inches='tight')
        
def plot_corner(method,
                samples_df,
                gal_name,
                distance):
    
    #Relative mass abundances for each component
    
    dgr_sCM20 = 0.17e-2
    dgr_lCM20 = 0.63e-3
    dgr_aSilM5 = 0.255e-2*2
    
    #Convert the samples dataframe back into an array for plotting,
    #also append the column names for corner labels later.
    
    samples = np.zeros(samples_df.shape)
    
    i = 0
    corner_labels = []
    
    for col_name in samples_df.dtypes.index:
        
        col_values = samples_df[col_name].tolist()
        
        corner_labels.append(col_name)
        
        samples[:,i] = col_values
        i += 1
    
    #Start by converting the scaling factor to a hydrogen mass
    
    scaling_factor_idx = corner_labels.index("log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)")
    
    samples[:,scaling_factor_idx] *= 1e-23
    samples[:,scaling_factor_idx] *= (distance*1e6*3.0857e18)**2
    samples[:,scaling_factor_idx] *= 1.67e-27
    samples[:,scaling_factor_idx] /= 2e30
    samples[:,scaling_factor_idx] *= 4*np.pi
        
    #In the cases where we vary abundances, turn this into a more meaningful
    #dust mass for each component
    
    if method in ['ascfree','abundfree']:
        
        idx_sCM20 = corner_labels.index("log$_{10}$ M$_\mathregular{sCM20}$")
        idx_lCM20 = corner_labels.index("log$_{10}$ M$_\mathregular{lCM20}$")
        idx_aSilM5 = corner_labels.index("log$_{10}$ M$_\mathregular{aSilM5}$")
        
        corner_labels.append(r'log$_{10}$ M$_\mathregular{sCM20}$')
        corner_labels.append(r'log$_{10}$ M$_\mathregular{lCM20}$')
        corner_labels.append(r'log$_{10}$ M$_\mathregular{aSilM5}$')
        
        samples[:,idx_sCM20] = dgr_sCM20*samples[:,scaling_factor_idx]*samples[:,idx_sCM20]
        samples[:,idx_sCM20] = np.log10(samples[:,idx_sCM20])
        samples[:,idx_lCM20] = dgr_lCM20*samples[:,scaling_factor_idx]*samples[:,idx_lCM20]
        samples[:,idx_lCM20] = np.log10(samples[:,idx_lCM20])
        samples[:,idx_aSilM5] = dgr_aSilM5*samples[:,scaling_factor_idx]*samples[:,idx_aSilM5]
        samples[:,idx_aSilM5] = np.log10(samples[:,idx_aSilM5])
        
    #Finally, for all cases the final column is total dust mass
    
    if method == 'default':
        
        samples[:,scaling_factor_idx] = samples[:,scaling_factor_idx]*dgr_sCM20 + \
                                        samples[:,scaling_factor_idx]*dgr_lCM20 + \
                                        samples[:,scaling_factor_idx]*dgr_aSilM5
        
    if method in ['ascfree','abundfree']:
    
        samples[:,scaling_factor_idx] = 10**samples[:,idx_sCM20] + \
                                        10**samples[:,idx_lCM20] + \
                                        10**samples[:,idx_aSilM5]
    
    samples[:,scaling_factor_idx] = np.log10(samples[:,scaling_factor_idx])
    
    #Calculate median values to plot as 'truths' on the corner plot
    
    medians = np.median(samples,axis=0)
    samples_plot_range = [0.995]*samples.shape[1]
    
    corner.corner(samples, labels=corner_labels,
                  quantiles=[0.16,0.84],
                  show_titles=True,
                  truths=medians,
                  range=samples_plot_range,
                  truth_color='k')
    
    plt.savefig('../plots/corner/'+gal_name+'_'+method+'.png',
                bbox_inches='tight',
                dpi=150)
    plt.savefig('../plots/corner/'+gal_name+'_'+method+'.pdf',
                bbox_inches='tight')

def filter_convolve(flux,
                    wavelength,
                    filter_dict,
                    keys):
    
    #Convolve this SED with the various filters we have to give
    #a monochromatic flux
 
    filter_fluxes = []
    
    for key in keys:
                 
        #Convolve MBB with filter
        
        filter_flux = np.interp(filter_dict[key][0],wavelength,flux)
                 
        filter_fluxes.append( np.abs( (np.trapz(filter_dict[key][1]*filter_flux,filter_dict[key][0])/
                                      np.trapz(filter_dict[key][1],filter_dict[key][0])) ) )
    
    return np.array(filter_fluxes)
