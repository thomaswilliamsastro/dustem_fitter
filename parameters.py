# -*- coding: utf-8 -*-
"""
Parameters for THEMCMC

@author: Tom Williams

v1.00: 20181002.
"""

###Input Parameters###

fluxes = 'fluxes.csv' #Path to pandas dataframe containing fluxes
filters = ['Spitzer_3.6','Spitzer_4.5',
           'Spitzer_5.8','Spitzer_8.0',
           'Spitzer_24','Spitzer_70',
           'WISE_3.4','WISE_4.6',
           'WISE_12','WISE_22',
           'PACS_100','PACS_160',
           'SPIRE_250',
           'SCUBA2_450','SCUBA2_850'] #List of filters to include

method = 'ascfree' #Options are default, abundfree, ascfree and fit-z

###Output Parameters###

plot = True #Produce SED and corner plots
units = 'flux' #Units for SED plot. Either flux (Jy) or luminosity (Lsun)
skirt_output = False #Produce a SKIRT code snippet
dustem_output = False #Produce a DustEM code snippet

###MPI Settings###

mpi = False #Run w/MPI or not
mpi_processes = 5 #If running MPI, use this many processes
                  #(should be number of cores/4)