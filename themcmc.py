# -*- coding: utf-8 -*-
"""
General purpose THEMIS MCMC fitter

@author: Tom Williams

v1.00: 20180817.
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

#Argument parsing
import argparse

#OS I/O stuff
import os
import sys

sys.path.append(os.getcwd())

import numpy as np
import pandas as pd
import time

#THEMCMC imports

import themcmc_sampler
import plotting
import general

# #For plotting
# import matplotlib
# matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
# import matplotlib.pyplot as plt
# plt.rcParams['mathtext.fontset'] = 'cm'
# import corner
# 
# import dill
# import tarfile
# from multiprocessing import Pool

# 
# import numpy as np
# 
# import emcee

# import datetime
# from tqdm import tqdm

#Set up the argument parser

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Options for the THEMCMC fitter.')
parser.add_argument('--method',type=str,default='default',metavar='',
                    help="Method for fitting the data. Options are 'default', 'ascfree'")
parser.add_argument('--plot',type=bool,default=True,metavar='',
                    help="Plot SED and corner plot.")
parser.add_argument('--skirtoutput',type=bool,default=False,metavar='',
                    help="Write out SKIRT mix code snippet. (N.B. does not work with multiu!)")
parser.add_argument('--dustemoutput',type=bool,default=False,metavar='',
                    help="Write out DustEM GRAIN.dat file (N.B. does not work with multiu!).")
parser.add_argument('--fluxes',type=str,default='fluxes',metavar='',
                    help="File containing fluxes to fit.")
parser.add_argument('--filters',type=str,default='filters',metavar='',
                    help="File containing filters to include.")

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
        
        else:
            
            samples = themcmc_sampler.sample(method=args.method,
                                             flux_file=args.fluxes,
                                             filter_file=args.filters,
                                             gal_row=gal_row)
        

                
        if args.plot:
        
            #For error plotting later
            
            y_to_percentile_stars = []
            y_to_percentile_small = []
            y_to_percentile_large = []
            y_to_percentile_silicates = []
            y_to_percentile_total = []
            
            for isrf,\
                omega_star,\
                alpha,\
                y_sCM20,\
                y_lCM20,\
                y_aSilM5,\
                dust_scaling in samples[np.random.randint(len(samples), size=150)]:
                       
                small_grains,\
                    large_grains,\
                    silicates = read_sed(isrf,
                                         alpha)
                
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
            
        #Calculate the percentiles
        
        percentiles =  zip(*np.percentile(samples, [16, 50, 84],axis=0))
         
        isrf,\
            omega_star,\
            alpha,\
            y_sCM20,\
            y_lCM20,\
            y_aSilM5,\
            dust_scaling = map(lambda percentiles: (percentiles[1],
                           percentiles[2]-percentiles[1],
                           percentiles[1]-percentiles[0]),
                           percentiles)
            
        if plot:
                    
            #Calculate residuals
                       
            flux_model = filter_convolve(y_median_total)
            
            residuals = (obs_flux-flux_model)/obs_flux
            residual_err = obs_error/obs_flux
                
            residuals = np.array(residuals)*100
            residual_err = np.array(residual_err)*100
            
            #And residuals for the points not fitted
            
            flux_model_not_fit = filter_convolve_not_fit(y_median_total,
                                                         keys_not_fit)
            
            residuals_not_fit = (obs_flux_not_fit-flux_model_not_fit)/obs_flux_not_fit
            residual_err_not_fit = obs_error_not_fit/obs_flux_not_fit
                
            residuals_not_fit = np.array(residuals_not_fit)*100
            residual_err_not_fit = np.array(residual_err_not_fit)*100
            
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
                
            for i in range(len(keys_not_fit)):
                
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
                                'IRAS_100':0.2}[keys_not_fit[i]]
                
                obs_error_not_fit[i] += calib_uncert*obs_flux_not_fit[i]
            
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
            
            plt.errorbar(obs_wavelength,obs_flux,
                         yerr=obs_error,
                         c='k',
                         ls='none',
                         marker='.',
                         zorder=99)
            
            plt.errorbar(obs_wavelength_not_fit,
                         obs_flux_not_fit,
                         yerr=obs_error_not_fit,
                         c='r',
                         ls='none',
                         marker='.',
                         zorder=98)
            
            plt.xscale('log')
            plt.yscale('log')
            
            plt.xlim([1,1000])
            plt.ylim([0.5*10**np.floor(np.log10(np.min(obs_flux))-1),
                      10**np.ceil(np.log10(np.max(obs_flux))+1)])
            
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
            
            plt.errorbar(obs_wavelength,residuals,yerr=residual_err,c='k',marker='.',
                         ls='none',
                         zorder=99)
            plt.errorbar(obs_wavelength_not_fit,residuals_not_fit,
                         yerr=residual_err_not_fit,c='r',marker='.',
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
            
            plt.savefig('plots/png/sed/'+gal_name+'_scfree.png',
                        bbox_inches='tight',
                        dpi=150)
            plt.savefig('plots/pdf/sed/'+gal_name+'_scfree.pdf',
                        bbox_inches='tight')
        
        #Corner plot -- convert the scaling factor to a hydrogen mass
        
        samples[:,-1] *= 1e-23
        samples[:,-1] *= (distance*1000*3.0857e18)**2
        samples[:,-1] *= 1.67e-27
        samples[:,-1] /= 2e30
        samples[:,-1] *= 4*np.pi
        
        mass_h = samples[:,-1].copy()
        
        #We now have a hydrogen mass in solar masses! Use the
        #relative abundances to turn this into a dust mass
        
        dgr_sCM20 = 0.17e-2
        dgr_lCM20 = 0.63e-3
        dgr_aSilM5 = 0.255e-2*2
        
        samples[:,-1] = samples[:,-1]*dgr_sCM20*samples[:,3] + \
                        samples[:,-1]*dgr_lCM20*samples[:,4] + \
                        samples[:,-1]*dgr_aSilM5*samples[:,5]
        
        samples[:,-1] = np.log10(samples[:,-1])
        
        #Because of varying grain sizes the abundances aren't really
        #physically meaningful. Turn into a dust mass for each component
        
        samples[:,3] = dgr_sCM20*mass_h*samples[:,3]
        samples[:,3] = np.log10(samples[:,3])
        samples[:,4] = dgr_lCM20*mass_h*samples[:,4]
        samples[:,4] = np.log10(samples[:,4])
        samples[:,5] = dgr_aSilM5*mass_h*samples[:,5]
        samples[:,5] = np.log10(samples[:,5])
        
        #Recalculate percentiles
        
        percentiles =  zip(*np.percentile(samples[:,3:], [16, 50, 84],axis=0))
         
        m_sCM20,\
            m_lCM20,\
            m_aSilM5,\
            m_dust = map(lambda percentiles: (percentiles[1],
                           percentiles[2]-percentiles[1],
                           percentiles[1]-percentiles[0]),
                           percentiles)
            
        if plot:
        
            corner.corner(samples, labels=[r'log$_{10}$ U',
                                           r'$\Omega_\ast$',
                                           r'$\alpha_\mathregular{sCM20}$',
                                           r'log$_{10}$ M$_\mathregular{sCM20}$',
                                           r'log$_{10}$ M$_\mathregular{lCM20}$',
                                           r'log$_{10}$ M$_\mathregular{aSilM5}$',
                                           r'log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)'],
                          quantiles=[0.16,0.84],
                          show_titles=True,
                          truths=[isrf[0],
                                  omega_star[0],
                                  alpha[0],
                                  m_sCM20[0],
                                  m_lCM20[0],
                                  m_aSilM5[0],
                                  m_dust[0]],
                          range=[0.995,0.995,0.995,0.995,0.995,0.995,0.995],
                          truth_color='k')
            
            plt.savefig('plots/png/corner/'+gal_name+'_scfree.png',
                        bbox_inches='tight',
                        dpi=150)
            plt.savefig('plots/pdf/corner/'+gal_name+'_scfree.pdf',
                        bbox_inches='tight')
        
        #Finally, write out code snippets for dustEM and SKIRT, if requested
        
        if dustem_output:
            
            with open('GRAIN_orig.DAT', 'r') as file :
                filedata = file.read()
                
                filedata = filedata.replace('[1.000000]','%.6f' % 10**isrf[0])
                filedata = filedata.replace('[0.170E-02]', '%.3E' % (0.17e-2*y_sCM20[0]))
                filedata = filedata.replace('[-5.00E-00]', '%.2E' % (alpha[0]))
                filedata = filedata.replace('[0.630E-03]', '%.3E' % (0.63e-3*y_lCM20[0]))
                filedata = filedata.replace('[0.255E-02]', '%.3E' % (0.255e-2*y_aSilM5[0]))
                
                # Write the converted file out
                with open('dustem_output/GRAIN_'+gal_name+'.dat', 'w') as file:
                    file.write(filedata)
                    
        if skirt_output:
            
            #Default, unchanging parameters
            
            #sCM20
            
            m_h = 1.6737236e-27
            amin = 0.0004e-6
            amax = 4.9e-6
            a_t = 0.e-6
            a_c = 0.05e-6
            a_u = 0.0001e-6
            gamma = 1
            zeta = 0
            eta = 0
            C = 1.717e-43
            rho_sCM20 = 1.6
            rho_sCM20 *= 100**3/1000
            
            a = np.linspace(amin,amax,100000)
            
            #lCM20
            
            a0_lCM20 = 0.007e-6
            rho_lCM20 = 1.57
            rho_lCM20 *= 100**3/1000
            
            #aSilM5
            
            a0_aSilM5 = 0.008e-6
            rho_aSilM5 = 2.19
            rho_aSilM5 *= 100**3/1000
            
            #Calculate new proportionality factors
                    
            idx = np.where(a<=a_t)        
            f_ed = np.zeros(len(a))       
            f_ed = np.exp(-((a-a_t)/a_c)**gamma)        
            f_ed[idx] = 1        
            f_cv = (1+np.abs(zeta) * (a/a_u)**eta )**np.sign(zeta)       
            omega_a_sCM20 = a**-alpha[0] * f_ed * f_cv        
            m_d_n_h = 4/3*np.pi*rho_sCM20*np.trapz(a**3*omega_a_sCM20,a)
            
            prop_factor_sCM20 = y_sCM20[0]*0.17e-2/(m_d_n_h/m_h)
            
            
            omega_a_lCM20 = a**-1 * np.exp(-0.5*np.log(a/a0_lCM20)**2)        
            m_d_n_h = 4/3*np.pi*rho_lCM20*np.trapz(a**3*omega_a_lCM20,a)
            
            prop_factor_lCM20 = y_lCM20[0]*0.63e-3/(m_d_n_h/m_h)
            
            
            omega_a_aSilM5 = a**-1 * np.exp(-0.5*np.log(a/a0_aSilM5)**2)
            m_d_n_h = 4/3*np.pi*rho_aSilM5*np.trapz(a**3*omega_a_aSilM5,a)
    
            prop_factor_aSilM5 = y_aSilM5[0]*0.255e-2/(m_d_n_h/m_h)
            
            #Write these new values out
            
            with open('template.ski', 'r') as file :
                filedata = file.read()
                
                filedata = filedata.replace('[alpha]','-%.1f' % (alpha[0]))
                filedata = filedata.replace('[y_sCM20]','%.11e' % (prop_factor_sCM20))
                filedata = filedata.replace('[y_lCM20]','%.11e' % (prop_factor_lCM20))
                filedata = filedata.replace('[y_aSilM5]','%.11e' % (prop_factor_aSilM5))
                
                # Write the converted file out
                with open('skirt_output/template_'+gal_name+'.ski', 'w') as file:
                    file.write(filedata)
    
    print('Code complete, took %.2fm' % ( (time.time() - start_time)/60 ))
