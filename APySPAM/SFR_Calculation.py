# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 13:30:07 2021

@author: oryan

Calculates the SFR at different points of the interaction, and calculates the flux to add onto the tagging particles.

"""
import numpy as np

class SFR_Calculations:
    def SFR(Gas,mass1,mass2,r1,r2,Sep,h,time,Weights,n1,n,init_ages):
        """
        Calculates the total mass formed via star formation throughout the simulation time. This assumes an exponentially declining star formation history, and enhances it 
        using the beta parameter described in O'Ryan et al. paper. 

        Note, in the implementation in the code, we just calculate the star formation from enhancement. All other star formation will already be accounted for in the underlying
        GALAXEV stellar populations we have used.

        If you are using your own stellar population models, you will need to edit this file with the star formation prescription you are using and enhance it by the beta parameter.

        Parameters
        -----------
        Gas[Tuple]:
            Gas mass of both the primary and secondary galaxies.
        mass1[float64]:
            Mass of the primary galaxy in Simulation Units (x10^11 Solar Masses).
        mass2[float64]:
            Mass of the secondary galaxy in simulation units (x10^11 Solar Masses)
        r1[float64]:
            The radius of the primary galaxy in simulation units (x15kpc)
        r2[float64]:
            The radius of the secondary galaxy in simulation units (x15kpc)
        Sep[List]:
            A list of the 3D distance between the centres of the primary and secondary galaxy in simulation units. Note, upon entering this function they are the
            wrong way around as the seperations have been found using backward integration.
        h[float64]:
            The timestep in simulation units.
        time[float64]:
            The total time that the simulation took in simulation units (x564Myrs). User specified.
        Weights[List]:
            The weights calculated in the Gas distribution function. These weights are based upon the total gas mass assigned to each particle upon initialisation.
        n1[flat64]:
            The number of particles in the primary galaxy.
        n[float64]:
            The total number of particles in the whole simulation.
        init_ages[Tuple]:
            Tuple containing the ages of the primary and secondary galaxy's at the start of the simulation run.

        Return
        --------
        SFRs[List]:
            A list of the star formation rates of each particle at the final timestep.
        Population_Mass[Numpy Array]:
            The total mass formed at each aprticle at each timestep. To then be used in colour function to combine with further SEDs.
        """
        # First, initialise the array of the different population masses and SFRs.
        Population_Mass = np.zeros([n,int(time/h)])
        Age = np.zeros(int(time/h))
        SFRs = np.zeros(n)

        Sep = np.flip(Sep)
                
        # Initialise Constants
        DU = 15
        TU = 87e6   # One time unit is 87Myrs according to ArXiv: 1511.05041
        e_1 = 1.5
        e_2 = 1.5
        counter = 0
        Baryonic_Fraction = (1 + 0.3333)/7.1333
                
        Mass_Ratio_1 = mass2/mass1
        Mass_Ratio_2 = mass1/mass2
        
        # Now, we can calculate the SFR at each particle through each timestep.
        for i in range(Population_Mass.shape[1]):
            Distance_Ratio_1 = r1/Sep[i]
            Distance_Ratio_2 = r2/Sep[i]
            
            Age_1 = init_ages[0] + counter*(h*TU/1e9)
            Age_2 = init_ages[1] + counter*(h*TU/1e9)
            
            Starburst_Enhancement_prim = 0.25*(Mass_Ratio_1)*(Distance_Ratio_1**2)
            Starburst_Enhancement_sec =  0.25*(Mass_Ratio_2)*(Distance_Ratio_2**2)
           
            # Note, this is only the enhanced part of the star formation. To find total at any time step, would need to uncomment lines below and add on the SFR_N_Base.
            SFRs_1 = Starburst_Enhancement_prim*((1/(e_1**2))*Age_1*np.exp(-(Age_1)/e_1))*((Baryonic_Fraction*mass1*1e11)/1e9)  # 1e9 here transforms SFR from M_0/Gyr to M_0/yr
            SFRs_2 = Starburst_Enhancement_sec*((1/(e_2**2))*Age_2*np.exp(-(Age_2)/e_2))*((Baryonic_Fraction*mass2*1e11)/1e9)
            
            # SFR_1_Base = ((1/(e_1**2))*Age_1*np.exp(-(Age_1)/e_1))*((Baryonic_Fraction*mass1*1e11)/1e9)
            # SFR_2_Base = ((1/(e_2**2))*Age_2*np.exp(-(Age_2)/e_2))*((Baryonic_Fraction*mass2*1e11)/1e9)
                        
            if SFRs_1 > 0:
                Population_Mass[:n1,i] += Weights[:n1]*SFRs_1*h*TU
            else:
                pass
            
            if SFRs_2 > 0:
                Population_Mass[n1:,i] += Weights[n1:]*SFRs_2*h*TU
            else:
                pass
            
            counter += 1
                                
                
        SFRs[:n1] = Weights[:n1]*SFRs_1
        SFRs[n1:] = Weights[n1:]*SFRs_2
        
        return SFRs, Population_Mass
        
        