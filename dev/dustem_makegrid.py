# -*- coding: utf-8 -*-
"""
Create a DustEM grid

@author: Tom Williams
"""

#Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import numpy as np
import multiprocessing
import time

def run_dustem(alpha_isrf):
    
    i,j = alpha_isrf
    
    with open('data/GRAIN_orig.DAT', 'r') as grain_file:
        filedata = grain_file.read()
     
        # Replace the target string
        filedata = filedata.replace('5.00','%.2f' % i)
        filedata = filedata.replace('1.000000','%.6f' % 10**j)
         
        # Write the file out again
        with open('data/GRAIN_'+str(i)+'_'+str(j)+'.DAT', 'w') as grain_file:
            grain_file.write(filedata)
     
    #Run the DustEM code
     
    os.system('./src/dustem GRAIN_'+str(i)+'_'+str(j)+'.DAT '+str(i)+'_'+str(j)+'.RES')
     
    #Remove the GRAIN.DAT files
     
    os.remove('data/GRAIN_'+str(i)+'_'+str(j)+'.DAT')
    
if __name__ == '__main__':
    
    start = time.time()

    os.chdir(os.getcwd())
    
    #Move into the DustEM directory
    
    os.chdir('../dustem')
    
    #Check if DustEM has already been compiled
    
    if not os.path.exists('src/dustem'):
    
        #Edit DM constants.f90
        
        with open('src/DM_constants.f90', 'r') as constants_file:
            filedata = constants_file.read()
        
        # Replace the target string
        filedata = filedata.replace('/Users/lverstra/Desktop/dustem4.2', str(os.getcwd()))
        
        # Write the file out again
        with open('src/DM_constants.f90', 'w') as constants_file:
            constants_file.write(filedata)
            
        #And compile DustEM
        
        os.chdir('src')
        os.system('make clean')
        os.system('make')
        
        os.chdir('../')
        
    #Set up values for alpha_sCM20 and the ISRF strength
    
    alpha = [2.6,5.4]
    isrf = [-1, 3.5]
    
#     alpha = np.arange(2.6,5.41,0.01)
#     isrf = np.arange(-1,3.51,0.01)
    
    alpha_isrf = [(x,y) for x in alpha for y in isrf]
    
    #Set up multiprocessing pool with number of processors
    
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
     
    pool.map(run_dustem,alpha_isrf)
    
    print('Complete! Took %.2fs' % (time.time()-start))