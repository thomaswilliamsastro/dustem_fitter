# -*- coding: utf-8 -*-
"""
Create the SED grid for the varying
small carbon grain sizes

@author: Tom Williams
"""
from __future__ import absolute_import, print_function, division
import os
import numpy as np
import subprocess
import stat

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

os.chdir('/home/daedalusdata/c1625914/dustem/data')

alpha = np.arange(2.6,5.41,0.01)
isrf = np.arange(-1,3.51,0.01)

for i in alpha:
    for j in isrf:
        with open('GRAIN_orig.DAT', 'r') as file :
            filedata = file.read()
        
            # Replace the target string
            filedata = filedata.replace('5.00','%.2f' % i)
            filedata = filedata.replace('1.000000','%.6f' % 10**j)
            
            # Write the file out again
            with open('GRAIN.DAT', 'w') as file:
                file.write(filedata)
        
        #Write a script to run the DustEM code
        
        script_file = open('run.sh','w+')
        script_file.write('#!/bin/csh\n')
        script_file.write('#\n')   
        script_file.write('../src/dustem\n')    
        script_file.write('mv /home/daedalusdata/c1625914/dustem/out/SED.RES /home/daedalusdata/c1625914/dustem_fitter/dev/grid/SED_%.2f_%.2f.RES' %(i,j))  
        
        script_file.close()
        
        st = os.stat('./run.sh')
        os.chmod('./run.sh', st.st_mode | stat.S_IEXEC)
          
        for path in execute(["./run.sh"]):
            print(path, end="")
    
      
print('Complete!')