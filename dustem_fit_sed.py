# -*- coding: utf-8 -*-
"""
Fit the M33 dust SED

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import matplotlib
matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'cm'
import numpy as np
from scipy.constants import h,k,c
import os
import emcee
import corner
import time
import datetime
from collections import OrderedDict
from multiprocessing import Pool

def read_sed(alpha):
    
    #Read in the SED that corresponds to this
    #alpha (rounded to nearest 0.01)
    
    wavelength,\
        small_grains,\
        large_grains,\
        pyroxine,\
        olivine, = np.loadtxt('SED_%.2f.RES' % alpha,
                              unpack=True,
                              skiprows=8,
                              usecols=(0,1,2,3,4))
        
    silicates = pyroxine+olivine
    
    #Convert these into units more like Jy -- divide by Hz
    
    frequency = 3e8/(wavelength*1e-6)
    
    small_grains /= frequency
    large_grains /= frequency
    silicates /= frequency
    
    #Add in a blackbody @ 5000K to represent the stars
    
    stars = 2*h*frequency**3/c**2 * (np.exp( (h*frequency)/(k*5000) ) -1)**-1
    
    #Normalise based on earlier value
    
    ratio = 6.97041442e+29
    
    stars /= ratio
    
    return wavelength,small_grains,large_grains,silicates,stars

def filter_convolve(wavelength,flux,
                    filter_dict):
    
    #Convolve this SED with the various filters we have to give
    #a monochromatic flux
 
    filter_fluxes = []
    
    for key in filter_dict:
                 
        #Convolve MBB with filter
        
        filter_flux = np.interp(filter_dict[key][0],wavelength,flux)
                 
        filter_fluxes.append(np.trapz(filter_dict[key][1]*filter_flux,filter_dict[key][0])/
                            np.trapz(filter_dict[key][1],filter_dict[key][0]))
        
    filter_fluxes = np.abs(np.array(filter_fluxes))
    
    return filter_fluxes

def lnprob(theta,
           obs_flux,obs_err,
           filter_dict):
    
    lp = priors(theta)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,
                       obs_flux,obs_err,
                       filter_dict)

def priors(theta):
    
    #alpha must be between 2.6 and 5.4,
    #insist that all the multiplicative
    #factors are greater than 0
    
    num_h,\
        omega_star,\
        alpha,\
        y_sCM20,\
        y_lCM20,\
        y_aSilM5 = theta
    
    if 2.6<=alpha<=5.4 and num_h > 0 and \
        omega_star > 0 and y_sCM20 > 0 and \
        y_lCM20 > 0 and y_aSilM5 > 0:
        return 0.0
    else:
        return -np.inf
    
def lnlike(theta,
           obs_flux,obs_err,
           filter_dict):
    
    num_h,\
        omega_star,\
        alpha,\
        y_sCM20,\
        y_lCM20,\
        y_aSilM5 = theta
        
    wavelength,\
        small_grains,\
        large_grains,\
        silicates,\
        stars = read_sed(alpha)    
    
    #Scale everything!
    
    total = y_sCM20*small_grains+\
            y_lCM20*large_grains+\
            y_aSilM5*silicates+omega_star*stars
            
    #Include number of hydrogen atoms
    
    total *= num_h
            
    filter_fluxes = filter_convolve(wavelength,total,
                                    filter_dict)
    
    chisq = np.sum( (filter_fluxes-obs_flux)**2 / obs_err**2)
    
    likelihood = -0.5*chisq
    
    return likelihood

start_time = time.time()

m33_distance = 840 #kpc

os.chdir('/home/daedalusdata/c1625914/dustem/out/themis_grid')

#Read in the M33 wavelengths, fluxes and flux errors

obs_wavelength,obs_flux,obs_error = np.loadtxt('obs_flux.txt',
                                               unpack=True)

wavelength = np.loadtxt('SED_5.40.RES',
                        unpack=True,
                        skiprows=8,
                        usecols=(0))

filter_dict = OrderedDict()

filter_folder = '/home/daedalusdata/c1625914/UsefulCode/sed_fitter/filters/'
filters = ['IRAC3_6','IRAC4_5','WISE4_6','IRAC5_8','IRAC8',
           'WISE12','WISE22','MIPS24','MIPS70','PACS100','PACS160',
           'SPIRE250','SPIRE350','SPIRE500','planck353']

for filter in filters:
    
    filter_wavelength,transmission = np.loadtxt(filter_folder+filter+'.dat',
                                     unpack=True)
    
    filter_wavelength /= 1e4
    transmission /= np.max(transmission)
    
    filter_dict[filter] = filter_wavelength,transmission

#Set up the MCMC. We have 6 free parameters, overall scaling factor,
#stellar scaling, power-law slope, and the ratio of the small- and 
#large-carbon grains and silicates. Set the initial guesses at the default
#THEMIS parameters
         
ndim,nwalkers = (6,500)
 
pos = []
 
for i in range(nwalkers):
    
    scaling_var = np.abs(np.random.normal(loc=1e40,scale=1e38))
    omega_star_var = np.abs(np.random.normal(loc=1,scale=1e-2))
    alpha_var = np.abs(np.random.normal(loc=5,scale=1e-2*5))
    y_sCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
    y_lCM20_var = np.abs(np.random.normal(loc=1,scale=1e-2))
    y_aSilM5_var = np.abs(np.random.normal(loc=1,scale=1e-2))

    pos.append([scaling_var,
                omega_star_var,
                alpha_var,
                y_sCM20_var,
                y_lCM20_var,
                y_aSilM5_var])
    
#Run this MCMC.

pool = Pool()

sampler = emcee.EnsembleSampler(nwalkers, 
                                ndim, 
                                lnprob, 
                                args=(obs_flux,obs_error,
                                      filter_dict),
                                pool=pool)
 
#500 steps for the 500 walkers, but throw away
#the first 250 as a burn-in

nsteps = 500
for i,result in enumerate(sampler.sample(pos,
                                         iterations=nsteps)):
    pos,probability,state = result
    print('MCMC run {:.1f}%'.format(float(i)*100/nsteps)+' complete',
          end='\r')
    
pool.close()
    
samples = sampler.chain[:, 250:, :].reshape((-1, ndim))
# samples = sampler.chain.reshape((-1, ndim))

print(str(datetime.datetime.now())+': MCMC Complete! Took %.2fs' % (time.time() - start_time))

#For error plotting later

y_to_percentile_stars = []
y_to_percentile_small = []
y_to_percentile_large = []
y_to_percentile_silicates = []
y_to_percentile_total = []

for scaling,\
    omega_star,\
    alpha,\
    y_sCM20,\
    y_lCM20,\
    y_aSilM5 in samples[np.random.randint(len(samples), size=250)]:
    
    wavelength,\
        small_grains,\
        large_grains,\
        silicates,\
        stars = read_sed(alpha)
    
    x,y = plt.plot(wavelength,omega_star*stars*scaling)[0].get_data()    
    y_to_percentile_stars.append(y)
    
    x,y = plt.plot(wavelength,y_sCM20*small_grains*scaling)[0].get_data()    
    y_to_percentile_small.append(y)
    
    x,y = plt.plot(wavelength,y_lCM20*large_grains*scaling)[0].get_data()    
    y_to_percentile_large.append(y)
    
    x,y = plt.plot(wavelength,y_aSilM5*silicates*scaling)[0].get_data()    
    y_to_percentile_silicates.append(y)
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
y_upper = y_upper_stars+y_upper_small+y_upper_large+y_upper_silicates
y_lower = y_lower_stars+y_lower_small+y_lower_large+y_lower_silicates

#Calculate the percentiles

percentiles =  zip(*np.percentile(samples, [16, 50, 84],axis=0))
 
scaling,\
    omega_star,\
    alpha,\
    y_sCM20,\
    y_lCM20,\
    y_aSilM5 = map(lambda percentiles: (percentiles[1],
                   percentiles[2]-percentiles[1],
                   percentiles[1]-percentiles[0]),
                   percentiles)
    
#From this, get the median fit lines

wavelength,\
    small_grains,\
    large_grains,\
    silicates,\
    stars = read_sed(alpha[0]) 

#Now scale everything!

total = y_sCM20[0]*small_grains+\
        y_lCM20[0]*large_grains+\
        y_aSilM5[0]*silicates+omega_star[0]*stars
        
total *= scaling[0]
        
#Calculate residuals
           
flux_model = filter_convolve(wavelength,total,
                             filter_dict)
    
residuals = (obs_flux-flux_model)/obs_flux
residuals *= 100
residual_err = obs_error/obs_flux
residual_err *= 100

fig1 = plt.figure(figsize=(10,6))
frame1 = fig1.add_axes((.1,.3,.8,.6))

#Plot the best fit and errorbars

#Observed fluxes

plt.errorbar(obs_wavelength,obs_flux,
             yerr=obs_error,
             c='k',
             ls='none',
             marker='.')

#THEMIS models

plt.fill_between(wavelength,y_lower_stars,y_upper_stars,
                 facecolor='m', interpolate=True,lw=0.5,
                 edgecolor='m', alpha=0.2)
plt.fill_between(wavelength,y_lower_small,y_upper_small,
                 facecolor='b', interpolate=True,lw=0.5,
                 edgecolor='b', alpha=0.2)
plt.fill_between(wavelength,y_lower_large,y_upper_large,
                 facecolor='g', interpolate=True,lw=0.5,
                 edgecolor='g', alpha=0.2)
plt.fill_between(wavelength,y_lower_silicates,y_upper_silicates,
                 facecolor='r', interpolate=True,lw=0.5,
                 edgecolor='r', alpha=0.2)
plt.fill_between(wavelength,y_lower,y_upper,
                 facecolor='k', interpolate=True,lw=0.5,
                 edgecolor='k', alpha=0.3)

plt.plot(wavelength,omega_star[0]*stars*scaling[0],
         c='m',
         label='Stars')
plt.plot(wavelength,y_sCM20[0]*small_grains*scaling[0],
         c='b',
         label='sCM20')
plt.plot(wavelength,y_lCM20[0]*large_grains*scaling[0],
         c='g',
         label='lCM20')
plt.plot(wavelength,y_aSilM5[0]*silicates*scaling[0],
         c='r',
         label='aSilM5')
plt.plot(wavelength,total,
         c='k',
         label='Total')

plt.xscale('log')
plt.yscale('log')

plt.xlim([1,1000])
plt.ylim([0.5*10**np.floor(np.log10(np.min(obs_flux))-1),
          10**np.ceil(np.log10(np.max(obs_flux))+1)])

plt.legend(loc='upper left',
           frameon=False)

plt.ylabel(r'$F_\nu$ (Jy)',
           fontsize=14)

plt.yticks(fontsize=14)

plt.tick_params(labelbottom=False)

#Add in residuals

frame2=fig1.add_axes((.1,.1,.8,.2))

plt.errorbar(obs_wavelength,residuals,yerr=residual_err,c='k',marker='.',
             ls='none')

plt.axhline(0,ls='--',c='k')

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.xscale('log')

plt.xlabel(r'$\lambda$ ($\mu$m)',
           fontsize=14)
plt.ylabel('Residual (%)')

plt.xlim([1,1000])
plt.ylim([-100,100])

plt.savefig('/home/daedalusdata/c1625914/Images/201806/THEMIS_sed.png',
            bbox_inches='tight',
            dpi=150)
plt.savefig('/home/daedalusdata/c1625914/Images/201806/THEMIS_sed.pdf',
            bbox_inches='tight',
            dpi=150)

#Corner plot -- convert the scaling factor to a hydrogen mass

samples[:,0] *= 1e-23
samples[:,0] *= (m33_distance*1000*3.0857e18)**2
samples[:,0] *= 1.67e-27
samples[:,0] /= 2e30
samples[:,0] *= 4*np.pi

mass_h = samples[:,0].copy()

#We now have a hydrogen mass in solar masses! Use the
#relative abundances to turn this into a dust mass

dgr_sCM20 = 0.17e-2
dgr_lCM20 = 0.63e-3
dgr_aSilM5 = 0.255e-2*2

samples[:,0] = samples[:,0]*dgr_sCM20*samples[:,3] + \
                samples[:,0]*dgr_lCM20*samples[:,4] + \
                samples[:,0]*dgr_aSilM5*samples[:,5]

samples[:,0] = np.log10(samples[:,0])

#Because of varying grain sizes the abundances aren't really
#physically meaningful. Turn into a dust mass for each component

samples[:,3] = dgr_sCM20*mass_h*samples[:,3]
samples[:,3] = np.log10(samples[:,3])
samples[:,4] = dgr_lCM20*mass_h*samples[:,4]
samples[:,4] = np.log10(samples[:,4])
samples[:,5] = dgr_aSilM5*mass_h*samples[:,5]
samples[:,5] = np.log10(samples[:,5])

#And log the stellar scaling parameter

samples[:,1] = np.log10(samples[:,1])

#Recalculate percentiles

percentiles =  zip(*np.percentile(samples, [16, 50, 84],axis=0))
 
m_dust,\
    omega_star,\
    alpha,\
    y_sCM20,\
    y_lCM20,\
    y_aSilM5 = map(lambda percentiles: (percentiles[1],
                   percentiles[2]-percentiles[1],
                   percentiles[1]-percentiles[0]),
                   percentiles)

corner.corner(samples, labels=[r'log$_{10}$ M$_\mathregular{dust}$ (M$_\odot$)',
                               r'log$_{10}$ $\Omega_\ast$',
                               r'$\alpha_\mathregular{sCM20}$',
                               r'log$_{10}$ M$_\mathregular{sCM20}$',
                               r'log$_{10}$ M$_\mathregular{lCM20}$',
                               r'log$_{10}$ M$_\mathregular{aSilM5}$'],
              quantiles=[0.16,0.84],
              show_titles=True,
              truths=[m_dust[0],
                      omega_star[0],
                      alpha[0],
                      y_sCM20[0],
                      y_lCM20[0],
                      y_aSilM5[0]],
              range=[0.995,0.995,0.995,0.995,0.995,0.995],
              truth_color='k')

# plt.show()

plt.savefig('/home/daedalusdata/c1625914/Images/201806/THEMIS_corner.png',
            bbox_inches='tight',
            dpi=150)
plt.savefig('/home/daedalusdata/c1625914/Images/201806/THEMIS_corner.pdf',
            bbox_inches='tight',
            dpi=150)

#Finally, write out any values we might want

percentiles_to_write = np.vstack((m_dust,
                                  omega_star,
                                  alpha,
                                  y_sCM20,
                                  y_lCM20,
                                  y_aSilM5))

np.savetxt('percentiles.txt',
           percentiles_to_write,
           header='log10(mdust)\nomega_star\nalpha_sCM20\ny_sCM20\ny_lCM20\ny_aSilM5')

print('Code Complete! Took %.2fs' % (time.time() - start_time))