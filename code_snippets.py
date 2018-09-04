# -*- coding: utf-8 -*-
"""
Code snippets for DustEM and SKIRT

@author: Tom Williams

v0.00.
"""

def dustemoutput(method,
                 samples):
    
    with open('GRAIN_orig.DAT', 'r') as file:
        
        #Read in the template file
        filedata = file.read()
        
        #Calculate the sample medians
        medians = np.median(samples,axis=0)
        
        #In all cases, we need to replace the ISRF strength
        filedata = filedata.replace('[1.000000]','%.6f' % 10**medians[0])
        
        #In the cases where we vary abundances, update as appropriate
        
        if method == 'ascfree':
            
            filedata = filedata.replace('[0.170E-02]', '%.3E' % (0.17e-2*medians[-4]))       
            filedata = filedata.replace('[0.630E-03]', '%.3E' % (0.63e-3*medians[-3]))
            filedata = filedata.replace('[0.255E-02]', '%.3E' % (0.255e-2*medians[-2]))
            
        #Where we vary alpha_sCM20, update this
        
        if method == 'ascfree':
            
            filedata = filedata.replace('[-5.00E-00]', '%.2E' % (medians[2]))
    
    # Write the converted file out
    with open('dustem_output/GRAIN_'+gal_name+'.dat', 'w') as file:
        file.write(filedata)
                    
def skirtoutput(method,
                samples):
    
    #REMEMBER WE CAN ONLY TRUST ALPHA TO 2 DP!
    
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