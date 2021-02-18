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
    
    def SSPs(metals):
        # This ones a little bit more complicated than the previous import. Here, we will require the user inputted metallicities as well
        # as the way to call the function which interpolates them. This shouldn't be too complicated, but I'm doing this just before JC
        # so will leave this as an exercise for later.
        pass
    
        