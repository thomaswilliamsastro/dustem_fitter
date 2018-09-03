# -*- coding: utf-8 -*-
"""
General useful functions

@author: Tom Williams

v0.00.
"""

def read_sed(isrf,
             alpha):
    
    #Read in the SED that corresponds to this
    #combination of ISRF strength and alpha 
    #(rounded to nearest 0.01)
    
    if '%.2f' % isrf in ['-0.00']:
        isrf = np.abs(isrf)

    small_grains = sCM20_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    large_grains = lCM20_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    silicates = aSilM5_df['%.2f,%.2f' % (alpha,isrf)].values.copy()
    
    #Convert these into units more like Jy -- divide by Hz
    
    small_grains /= frequency
    large_grains /= frequency
    silicates /= frequency
    
    return small_grains,large_grains,silicates

def filter_convolve(flux):
    
    #Convolve this SED with the various filters we have to give
    #a monochromatic flux
 
    filter_fluxes = []
    
    for key in keys:
                 
        #Convolve MBB with filter
        
        filter_flux = np.interp(filter_dict[key][0],wavelength,flux)
                 
        filter_fluxes.append( np.abs( (np.trapz(filter_dict[key][1]*filter_flux,filter_dict[key][0])/
                                      np.trapz(filter_dict[key][1],filter_dict[key][0])) ) )
    
    return np.array(filter_fluxes)