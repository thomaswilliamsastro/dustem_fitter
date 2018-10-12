# -*- coding: utf-8 -*-
"""
Parameters for THEMCMC

@author: Tom Williams
"""

###Input Parameters###

fluxes = 'fluxes_m33_px.csv' #Path to pandas dataframe containing fluxes
filters = ['Spitzer_3.6','Spitzer_4.5',
           'Spitzer_5.8','Spitzer_8.0',
           'Spitzer_24','Spitzer_70',
           'WISE_3.4','WISE_4.6',
           'WISE_12','WISE_22',
           'PACS_100','PACS_160',
           'SPIRE_250',
           'SCUBA2_450','SCUBA2_850'] #List of filters to include

###Dust Model Parameters###

method = 'ascfree' #Options are default, abundfree and ascfree
components = 2 #Number of dust components to fit
overwrite_samples = False #Rerun the MCMC if a samples file already exists

###Output Parameters###

plot_sed = False #Produce SED and corner plots
units = 'flux' #Units for SED plot. Either flux (Jy) or luminosity (Lsun)
overwrite_sed_plot = True #If SED already exists, overwrite

plot_corner = False
overwrite_corner_plot = True #If corner plot already exists, overwrite

skirt_output = False #Produce a SKIRT code snippet
dustem_output = False #Produce a DustEM code snippet

overwrite_plots = True #Overwrite any already created plots

###MPI Settings###

mpi = True #Run w/MPI or not
mpi_processes = 5 #If running MPI, use this many processes
                  #(should be number of cores/4)
