# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:12:54 2021

@author: oryan

This algorithm handles all of the imports that APySPAM will require. It takes no inputs and returns two arrays: 
one is an array of the filter values, and the other is an array of SSPs. 

Doing it this way allows us to only load once, and to check that the correct folders exist. 
"""

import os
import numpy as np
import glob
from SEDs import SED

class Imports():
    def Filters():
        dir_path = os.getcwd()
        filter_folder = dir_path+'\\Filters\\'
        filters = glob.glob(filter_folder+'*.*')
        
        filter_data = []
        
        for i in range(len(filters)):
            temp = np.loadtxt(filters[i])
            filter_data.append(temp)
            
        return filter_data
    
    def SSPs(metallicity):
        dir_path = os.getcwd()
        metals_folder = dir_path+'\\Spectral_Data\\'
        
        Spectral_Density_Array_1 = SED.Metal_Interpolation(metallicity[0],metals_folder)
        Spectral_Density_Array_2 = SED.Metal_Interpolation(metallicity[1],metals_folder)
        
        Wavelength = Spectral_Density_Array_1[1:,0]
        
        return Spectral_Density_Array_1, Spectral_Density_Array_2, Wavelength
    
    def Export(Vectors,Fluxes,Formed_Stellar_Masses,SFRs,directory_result):
        pass
    
        