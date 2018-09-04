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

#THEMCMC imports

import general

def plot_sed(method,
             flux_df,
             filter_df,
             gal_row,
             samples,
             filter_dict):
    
    gal_name = flux_df['name'][gal_row]
    
    #Read in the DustEM grid
    
    sCM20_df = pd.read_hdf('models.h5','sCM20')
    lCM20_df = pd.read_hdf('models.h5','lCM20')
    aSilM5_df = pd.read_hdf('models.h5','aSilM5')
    
    wavelength = sCM20_df['wavelength'].values.copy()
    frequency = 3e8/(wavelength*1e-6)
    
    #Generate stars
    
    stars = general.define_stars(flux_df,
                                 gal_row,
                                 frequency)
    
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
               
        small_grains,\
            large_grains,\
            silicates = general.read_sed(isrf,
                                         alpha,
                                         sCM20_df,
                                         lCM20_df,
                                         aSilM5_df,
                                         frequency)
        
        x,y = plt.plot(wavelength,omega_star*stars)[0].get_data()    
        y_to_percentile_stars.append(y)
        
        x,y = plt.plot(wavelength,y_sCM20*small_grains*dust_scaling)[0].get_data()    
        y_to_percentile_small.append(y)
        
        x,y = plt.plot(wavelength,y_lCM20*large_grains*dust_scaling)[0].get_data()    
        y_to_percentile_large.append(y)
        
        x,y = plt.plot(wavelength,y_aSilM5*silicates*dust_scaling)[0].get_data()    
        y_to_percentile_silicates.append(y)
        
        x,y = plt.plot(wavelength,dust_scaling*(y_sCM20*small_grains+\
                                    y_lCM20*large_grains+\
                                    y_aSilM5*silicates)+\
                                    omega_star*stars)[0].get_data()
        y_to_percentile_total.append(y)
        plt.close()
        
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
    #RMS, so include overall calibration from Chris' paper
    
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
              0.8*10**np.ceil(np.log10(np.max(obs_flux[obs_flag == 0]))+1)])
    
    plt.legend(loc='upper left',
               fontsize=14,
               frameon=False)
    
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
    
    plt.savefig('plots/sed/'+gal_name+'_'+method+'.png',
                bbox_inches='tight',
                dpi=150)
    plt.savefig('plots/sed/'+gal_name+'_'+method+'.pdf',
                bbox_inches='tight')
        
def plot_corner(method,
                samples,
                gal_name,
                distance):
    
    #Relative mass abundances for each component
    
    dgr_sCM20 = 0.17e-2
    dgr_lCM20 = 0.63e-3
    dgr_aSilM5 = 0.255e-2*2
    
    #Start by converting the scaling factor to a hydrogen mass
    
    samples[:,-1] *= 1e-23
    samples[:,-1] *= (distance*1000*3.0857e18)**2
    samples[:,-1] *= 1.67e-27
    samples[:,-1] /= 2e30
    samples[:,-1] *= 4*np.pi
    
    #In all cases, the first 2 parts of samples are the ISRF and
    #the stellar scaling factor, which we don't have to do anything
    #to
    
    corner_labels = [r'log$_{10}$ U',r'$\Omega_\ast$']
    
    #In the case of varying alpha_sC20, the third column is alpha which
    #doesn't require any conversions
    
    if method == 'ascfree':
    
        corner_labels.append(r'$\alpha_\mathregular{sCM20}$')
        
    #In the cases where we vary abundances, turn this into a more meaningful
    #dust mass for each component
    
    if method in ['ascfree','abundfree']:
        
        corner_labels.append(r'log$_{10}$ M$_\mathregular{sCM20}$')
        corner_labels.append(r'log$_{10}$ M$_\mathregular{lCM20}$')
        corner_labels.append(r'log$_{10}$ M$_\mathregular{aSilM5}$')
        
        samples[:,-4] = dgr_sCM20*samples[:,-1]*samples[:,-4]
        samples[:,-4] = np.log10(samples[:,-4])
        samples[:,-3] = dgr_lCM20*samples[:,-1]*samples[:,-3]
        samples[:,-3] = np.log10(samples[:,-3])
        samples[:,-2] = dgr_aSilM5*samples[:,-1]*samples[:,-2]
        samples[:,-2] = np.log10(samples[:,-2])
        
    #Finally, for all cases the final column is total dust mass
    
    corner_labels.append(r'log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)')
    
    if method == 'default':
        
        samples[:,-1] = samples[:,-1]*dgr_sCM20 + \
                        samples[:,-1]*dgr_lCM20 + \
                        samples[:,-1]*dgr_aSilM5
        
    if method in ['ascfree','abundfree']:
    
        samples[:,-1] = 10**samples[:,-4] + \
                        10**samples[:,-3] + \
                        10**samples[:,-2]
    
    samples[:,-1] = np.log10(samples[:,-1])
    
    #Calculate median values to plot as 'truths' on the corner plot
    
    medians = np.median(samples,axis=0)
    range = [0.995]*samples.shape[1]
    
    corner.corner(samples, labels=corner_labels,
                  quantiles=[0.16,0.84],
                  show_titles=True,
                  truths=medians,
                  range=range,
                  truth_color='k')
    
    plt.savefig('plots/corner/'+gal_name+'_'+method+'.png',
                bbox_inches='tight',
                dpi=150)
    plt.savefig('plots/corner/'+gal_name+'_'+method+'.pdf',
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