# -*- coding: utf-8 -*-
"""
This algorithm calculates the star formation rate and mass formed in a timestep using 
the Kennicutt-Schmidt model.
@author: oryan
"""

import numpy as np
import sys

class KS_Model:
    def SFR_KS(n1,n2,Coords,r1,r2,DU,Masses,Sec_Position,TU,timestep):
        n = n1 + n2
        t = TU*timestep*1e9
        Coords_temp = np.zeros(Coords.shape)
        Coords_temp += Coords
        Coords_temp *= DU
        Coords_temp = Coords_temp[:n,:]
        Formed_Mass = np.zeros(n)
        Mass_Used = np.zeros(n)
        SFR = np.zeros(n)
        # Total_Spanned_Volume = n1*((4/3)*np.pi*1e3**3)
        Prim_Area = np.pi*(DU*r1*1e3)**2
        Sec_Area = np.pi*(DU*r2*1e3)**2
        Particle_Area_1 = Prim_Area/n1
        Particle_Area_2 = Sec_Area/n2
        
        for p in range(n):
                                    
            Bin_Mass = Masses[p]
            # if p <= n1:
            #     Particle_Area = Particle_Area_1
            # elif p > n1:
            #     Particle_Area = Particle_Area_2
            
            Particle_Area = KS_Model.Criteria_Calculator(p,Coords_temp)
            
            Avg_Sigma_Density = (Bin_Mass)/Particle_Area
            # Avg_Sigma_Density = 0.01*Avg_Sigma_Density
            Avg_Sigma_SFR = 8.33e-4*(Avg_Sigma_Density)**1.4         # This is found in Solar Masses per kpc squared per year.
            Avg_SFR = 0.03*Avg_Sigma_SFR*(Particle_Area/1e6)                       # Get rid of the z direction
            # if p == 0:
            #     print(' ')
            #     print(Avg_Sigma_Density)
            #     print(Avg_SFR)
            Total_Mass_Formed = Avg_SFR*t
                        
            Formed_Mass[p] = Total_Mass_Formed
            Mass_Used[p] = Total_Mass_Formed
            
            SFR[p] = Avg_SFR
                    
        M1 = np.sum(Mass_Used[:n1])
        M2 = np.sum(Mass_Used[n1:])
        
        # print(' ')
        # print('The Total SFR in the Primary Galaxy is: ', np.sum(SFR[:n1]))
        # print('The Total SFR in the Secondary Galaxy is: ', np.sum(SFR[n1:]))
        # print(' ')
        # print('The Mean SFR in the Primary Galaxy is: ', np.mean(SFR[:n1]))
        # print('The Mean SFR in the Secondary Galaxy is: ', np.mean(SFR[n1:]))
        # sys.exit()
        
        return SFR,Formed_Mass, M1, M2, Mass_Used
    
    
    def Criteria_Calculator(p,Coords):
      criteria_r_index = np.zeros(len(Coords),dtype='int')
      Particle_Frame_Coords = Coords[:,:] - Coords[p,:]
      x = Particle_Frame_Coords[:,0]
      y = Particle_Frame_Coords[:,1]
      z = Particle_Frame_Coords[:,2]
    
      Radii = np.sqrt(x*x + y*y + z*z)

      criteria_r_index = [x[0] for x in enumerate(Radii) if x[1] <= 2.0]
      Nearest_Neighbours = Radii[criteria_r_index]
      
      if np.min(Nearest_Neighbours) == 0 and np.max(Nearest_Neighbours) == 0:
          Particle_Size = 1
      else:
          Particle_Size = np.min(Nearest_Neighbours[np.nonzero(Nearest_Neighbours)])/2
          
      Particle_Area = np.pi*(Particle_Size*1e3)**2
      
      return Particle_Area