# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:06:45 2020

@author: oryan
"""

import pandas as pd
import numpy as np

class Setup_Parameters:
  def Setup_Parameters(self,Input_Counter):
      Inputs = pd.read_csv(r'C:\Users\oryan\Documents\PySPAM\All_Gal_Inputs.csv')
      Inputs_Data = Inputs.values
      
      Galaxy_Inputs = Inputs_Data[Input_Counter,:]
      
      self.params.n = 2000
      
      self.params.phi1 = Galaxy_Inputs[10]      #NOTE: Original Value 5.0
      self.params.theta1 = Galaxy_Inputs[12] #NOTE: Original Value = 5.0
      self.params.rscale1 = [0.0,0.0,1.0] #NOTE: Original Value = [1.0, 1.0, 1.0]    #These parameters are only important in the soft point mass runs.
      self.params.rout1 = Galaxy_Inputs[8]    #NOTE: Original Value = 1.0
      self.params.mass1 = Galaxy_Inputs[6]    #NOTE: Original Value = 1.0
      self.params.epsilon1 = 0.3   #NOTE: Original Value = 0.3
      self.params.eps1 = self.params.epsilon1*self.params.epsilon1   
      self.params.heat1 = 0.0  #NOTE: Original Value 0.0
      self.params.opt1 = 1    #NOTE: Original Value 1
      self.params.Galactic_Age_1 = 10           #Specify this in Gyrs
      self.params.Gas_Fraction_1 = 0.17
      self.params.metal_1 = 0.005    #
      
      self.params.phi2 = Galaxy_Inputs[11]    #NOTE: Original Value = 0.0
      self.params.theta2 = Galaxy_Inputs[13]    #NOTE: Original Value = 0.0
      self.params.rscale2 = [0.0,0.0,0.0] #NOTE: Original Value = [0.3,0.3,0.3]
      self.params.rout2 = Galaxy_Inputs[9]  #NOTE: Original Value = 0.50
      self.params.mass2 = Galaxy_Inputs[7]   #NOTE:Have changed this from 0.50
      self.params.epsilon2 = 0.3   #NOTE: Original Value = 0.3
      self.params.eps2 = self.params.epsilon2*self.params.epsilon2
      self.params.heat2 = 0.0   #NOTE: Original Value = 0.0
      self.params.opt2 = 1      #NOTE: Original Value = 1
      self.params.Galactic_Age_2 = 8
      self.params.Gas_Fraction_2 = 0.14
      self.params.metal_2 = 0.05   #Set the metalicity of the particles in galaxy 2
      
      self.params.n1 = int(self.params.n*(self.params.mass1/(self.params.mass1 + self.params.mass2)))    #NOTE: Original Value 1000
      self.params.n2 = int(self.params.n*(self.params.mass2/(self.params.mass1 + self.params.mass2)))    #NOTE: Have changed this from 500
      
      while (self.params.n1 + self.params.n2) < self.params.n:
          random = np.random.random()
          if random <= 0.5:
              self.params.n1 += 1
          elif random > 0.5:
              self.params.n2 += 1
      
      r_x = Galaxy_Inputs[0]
      r_y = Galaxy_Inputs[1]
      r_z = Galaxy_Inputs[2]

      v_x = Galaxy_Inputs[3]
      v_y = Galaxy_Inputs[4]
      v_z = Galaxy_Inputs[5]
      
      self.params.sec_vec = [r_x, r_y, r_z, v_x, v_y, v_z]
      self.params.use_sec_vec = True

      self.params.inclination_degree = 20# -1.38 #-22.47999  #NOTE: Changed from 20
      self.params.omega_degree = 0 #0.41816  #NOTE: Original Value = 0.0
      self.params.rmin = 0.90            #NOTE:Changed from 0.90
      self.params.velocity_factor = 1#-1.38  #NOTE: Original Value = 0.90      

      self.params.h = 0.01
      self.params.tstart = Galaxy_Inputs[34]  #NOTE: Original Value = -5
      self.params.tend = 0#Galaxy_Inputs[35]#3.6642833366547523
      self.params.time = self.params.tstart - self.params.tend#-5.307844817363995  #NOTE: Original Value = -5
      self.params.tIsSet = False  #NOTE: Original Value = True
      
      self.params.redshift = Galaxy_Inputs[38]    #This varible will change the redshift of our galaxies. This will have to be looked up for object.
      
      self.params.Galaxy_Name = str(Galaxy_Inputs[37])          #This specifies folder name for the program. Have a folder: '\dir...\Colour_Evolution_Movie'+galaxyname+'\...
      
      self.params.display_scale = 2*Galaxy_Inputs[39] + 2  # In Galaxy units (x15kpc). This will change the scale of the overall interaction plot at the end.
                                        # WARNING!!! Making this too small may crash the code. Start big, and work down to better scales.