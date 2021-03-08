# -*- coding: utf-8 -*-
"""
This algorithm calculates the star formation rate and mass formed in a timestep using 
the Kennicutt-Schmidt model of star formation. For a further breakdown see [url].

Inputs: n1/n2 - [int] - The number of particles in the primary and secondary galaxies, respectively, at initialization.
        Coords - [numpy array] - A nx6 array containing the x,y,z,vx,vy,vz position of each particle.
        r1/r2 - [float] - The radii of the primary and secondary galaxy respectively.
        DU/TU - [float] - Conversion from galaxy units of time and distance to SI units.
        Masses - [2x1 numpy array] - Contains the primary and secondary masses of the galaxies.
        Sec_Position - [3x1 numpy array] - The position of the secondary galaxy centre at the current timestep.
        timestep - [float] - The length of the timestep, defined by self.params.h in Setup_Parameters.py.
        
Outputs: SFR - [nx1 numpy array] - The SFR at each particle.
         Formed Mass - [nx1 numpy array] - The mass formed due to any star formation enhancement in the particle.
         M1/M2 - [float] - Total mass formed at each particle.
         Mass_Used - The total gas mass used in the primary and secondary.

@author: oryan
"""

import numpy as np

class KS_Model:
    def SFR_KS(n1,n2,Coords,r1,r2,DU,Masses,Sec_Position,TU,timestep,Gal_m1,Gal_m2,Ages,e_times,Weights):
        n = n1 + n2
        t = TU*timestep*1e9
        Coords_temp = np.zeros(Coords.shape)
        Coords_temp += Coords
        Coords_temp *= DU
        Coords_temp = Coords_temp[:n,:]
        Formed_Mass = np.zeros(n)
        Mass_Used = np.zeros(n)
        SFR_Array = np.zeros(n)
        SFR_test = 0
        
        Baryonic_Fraction = (1 + 0.3333)/7.1333
        
        for p in range(n):
                                    
            Bin_Mass = Masses[p]
            
            Particle_Area = KS_Model.Criteria_Calculator(p,Coords_temp)
            
            if p <= n1:
                SFR_acc = Weights[p]*((1/(e_times[0]**2))*Ages[0]*np.exp(-(Ages[0])/e_times[0]))*((Baryonic_Fraction*Gal_m1*1e11)/1e9)  # 1e9 here transforms SFR from M_0/Gyr to M_0/yr
            elif p > n1:
                SFR_acc = Weights[p]*((1/(e_times[1]**2))*Ages[1]*np.exp(-(Ages[1])/e_times[1]))*((Baryonic_Fraction*Gal_m2*1e11)/1e9)
            
            Avg_Sigma_Density = (Bin_Mass)/Particle_Area
            Avg_Sigma_SFR = 8.33e-4*(Avg_Sigma_Density/10)**1.4         # This is found in Solar Masses per kpc squared per year. 
            KS_SFR = 0.03*Avg_Sigma_SFR*(Particle_Area/1e6)
            SFR = KS_SFR - SFR_acc
            if SFR < 0:
                SFR = 0
            Total_Mass_Formed = SFR*t
                        
            Formed_Mass[p] = Total_Mass_Formed
            Mass_Used[p] = Total_Mass_Formed
            
            SFR_Array[p] = KS_SFR
                    
        M1 = np.sum(Mass_Used[:n1])
        M2 = np.sum(Mass_Used[n1:])
        
        return SFR_Array,Formed_Mass, M1, M2, Mass_Used
    
    
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