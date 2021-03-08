"""
This algorithm calculates the star formation rate and mass formed in a timestep using 
a Delayed-Tau model. For more info see [url to MkDocs].

Inputs: Gas_Masses
        Gal_m1/m2 - The total masses of the primary and secondary galaxies.
        r1/r2 - The radii of the primary and secondary galaxies
        Galaxy_Seperation - The seperation between the two galaxy centres.
        DU/TU - Conversion from galaxy units to kpc/yr.
        Ages - Galactic Ages.
        e_times - e-folding times
        timestep - Timestep defined in simulation.
        Weights - Particle weights defined in setup.
        n1 - The number of particles in the primary galaxy. Used as the cut-off index in the arrays to define which particle is primary or secondary. 
        
Outputs: SFRs - The total star formation rate at the particle position. This is both from the underlying SSP (assuming a delayed-tau model). 
         Particle_Mass_Formed - The mass formed due to star formation enhancement. 
         M1 - Total mass formed at the particle position (sum of underlying SSP model and enhancement) for Primary Galaxy.
         M2 - As M1 but for secondary.

"""

import numpy as np
import sys

class Delayed_Tau:
    def SFR_DT(Gas_Masses,Gal_m1,Gal_m2,r1,r2,Galaxy_Seperation,DU,TU,Ages,e_times,timestep,i,
               Weights,n1):
        
        Age_1, Age_2 = Ages
        e_1, e_2 = e_times
        Gas_m1, Gas_m2 = Gas_Masses
        
        if i <= 100:
            delay = i
        else:
            delay = 100
            
        Seperation_Coords = Galaxy_Seperation[i-delay]*DU
        Seperation = np.sqrt(Seperation_Coords[0]*Seperation_Coords[0] + Seperation_Coords[1]*Seperation_Coords[1] + Seperation_Coords[2]*Seperation_Coords[2])
        Distance_Ratio_1 = r1*DU/(Seperation)
        Distance_Ratio_2 = r2*DU/(Seperation)
        Mass_Ratio_prim = Gal_m2/Gal_m1       # This was originally the wrong way around. Remember, we want the effect that the secondary has ON the primary.
        Mass_Ratio_sec = Gal_m1/Gal_m2
        Baryonic_Fraction = (1 + 0.3333)/7.1333
            
        Starburst_Enhancement_prim = 1 + (Mass_Ratio_prim)*(Distance_Ratio_1**2)
        Starburst_Enhancement_sec =  1 + (Mass_Ratio_sec)*(Distance_Ratio_2**2)
           
        SFR_1_Total = Starburst_Enhancement_prim*((1/(e_1**2))*Age_1*np.exp(-(Age_1)/e_1))*((Baryonic_Fraction*Gal_m1*1e11)/1e9)  # 1e9 here transforms SFR from M_0/Gyr to M_0/yr
        SFR_2_Total = Starburst_Enhancement_sec*((1/(e_2**2))*Age_2*np.exp(-(Age_2)/e_2))*((Baryonic_Fraction*Gal_m2*1e11)/1e9)
        
        SFR_1_Base = ((1/(e_1**2))*Age_1*np.exp(-(Age_1)/e_1))*((Baryonic_Fraction*Gal_m1*1e11)/1e9)
        SFR_2_Base = ((1/(e_2**2))*Age_2*np.exp(-(Age_2)/e_2))*((Baryonic_Fraction*Gal_m2*1e11)/1e9)
        SFR_Bases = [SFR_1_Base, SFR_2_Base]
        
        SFR_Enhanced = np.asarray([SFR_1_Total - SFR_1_Base,SFR_2_Total - SFR_2_Base])
        
        print(SFR_1_Total)
        
        SFR_Enhanced[SFR_Enhanced <= 0] = 0
        
        SFRs,Particle_Mass_Formed,M1,M2 = Delayed_Tau.Mass_Formed(SFR_Enhanced,SFR_Bases,timestep,TU,Weights,n1)
        
        yr_conv = 1e9
        M1 += (SFR_1_Base*timestep*TU*yr_conv)
        M2 += (SFR_2_Base*timestep*TU*yr_conv)
                
        return SFRs,Particle_Mass_Formed,M1,M2      # Note that in this algorithm, we are finding the TOTAL STAR FORMATION RATE of at the position
                                                    # of the particle, but only finding the extra mass formed from enhancement. The underlying SSP
                                                    # deals with the un-enhanced SFR in the population.
    def Mass_Formed(SFR,SFR_Bases,timestep,TU,Weights,n1):
        SFR_1, SFR_2 = SFR
        SFR_1_Base, SFR_2_Base = SFR_Bases
        yr_conv = 1e9
        
        M1 = Weights[:n1]*SFR_1*timestep*TU*yr_conv
        M2 = Weights[n1:]*SFR_2*timestep*TU*yr_conv
        
        Particle_SFR_1 = (Weights[:n1]*SFR_1) + (Weights[:n1]*SFR_1_Base)
        Particle_SFR_2 = (Weights[n1:]*SFR_2) + (Weights[n1:]*SFR_2_Base)
        
        SFRs = np.concatenate([Particle_SFR_1,Particle_SFR_2])
        MT = np.concatenate([M1,M2])
        
        sum_m1 = np.sum(M1)
        sum_m2 = np.sum(M2)
        
        return SFRs,MT,sum_m1,sum_m2
    
    