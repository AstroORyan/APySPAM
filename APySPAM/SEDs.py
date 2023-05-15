# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 17:54:48 2021

@author: oryan
"""
import numpy as np
import datetime
import sys

class SED:
    def getSED(Spectral_Density_1,Spectral_Density_2,Ages,n1,n2,time,Weights,Part_Mass,SFR_Mass,h):
        '''
        Calculates the SED of each particle from both the initial underlying stellar population and any populations formed in star formation at each timestep.

        Parameters
        ----------
        Spectral_Density_1[np.array]:
            Array containing the SED information of the Bruzual and Charlot Model found with the relevent metallicity of the primary galaxy. This array contains flux vs wavength at different ages of simple stellar population.
            Should be about [200, 6900] dim.
        Spectral_Density_2[np.array]:
            As above, but for the secondary galaxy.
        Ages[array]:
            N_particles x 1 array containing the intial ages of each particle. Note, this is assumed to be the user given age of the primary and secondary galaxies respectively.
        n1[int]:
            Number of particles in the primary galaxy.
        n2[int]:
            Number of particles in the secondary galaxy.
        time[float]:
            The total backwards integration time of the interaction. In simulation units.
        Weights[array]:
            The weights assigned to each particle based on the gas distribution.
        Part_Mass[array]:
            The total stellar mass assigned to each particle. Based upon the weights assigned upon initialisation and the gas fraction.
        SFR_Mass[array]:
            The total stellar mass formed in star formation at every timestep.
        h[float]:
            The length of the timestep in simulation units.

        Returns
        -------
        Tuple:
            Tuple containing the final combined SED of each particle from the initial underlying stellar population and any new stellar populations formed in star formation. Also returns the wavelengths of the SEDs.
        '''
        Initial_Flux_Dist = np.zeros([n1+n2,6900])
        Additional_Flux = np.zeros([n1+n2,6900])
        
        #Load in original SED.
        Initial_Flux_Dist[:n1,:] = SED.initSED(Spectral_Density_1,Ages[0],n1)
        Initial_Flux_Dist[n1:,:] = SED.initSED(Spectral_Density_2,Ages[1],n2)

        Weights = Weights[:,np.newaxis]
        
        # Multiply by the mass each particle has.
        Initial_Flux_Dist[:n1,:] = Initial_Flux_Dist[:n1,:]*(Weights[:n1,:])*(Part_Mass[0]*1e11)
        Initial_Flux_Dist[n1:,:] = Initial_Flux_Dist[n1:,:]*(Weights[n1:,:])*(Part_Mass[1]*1e11)
        
        # Now, the really rough part... Finding the contribution from all the mass formed. 
        Additional_Flux[:n1,:] = SED.SFR_Flux(Spectral_Density_1[0,:],Spectral_Density_1[1:,1:],Ages[0],n1,time,SFR_Mass[:n1,:],h)
        Additional_Flux[n1:,:] = SED.SFR_Flux(Spectral_Density_2[0,:],Spectral_Density_2[1:,1:],Ages[1],n2,time,SFR_Mass[n1:,:],h)

        Flux_Dist = Initial_Flux_Dist + Additional_Flux
        wavelength = Spectral_Density_1[1:,0]
        
        return [Flux_Dist,wavelength]
                
        
    def initSED(SEDs,Age,n):
        '''
        Initialise the Spectral Energy Distribution of each particle. Finds the total flux of the initial underlying stellar population BEFORE accounting for any stellar mass found in star formation. 

        Parameters
        -----------
        SEDs[array]
            The SED provided by the Bruzual and Charlot Model. Is a [200, 6900] array of wavelength vs flux at different timesteps. 
        Age[float]:
            The age of the initial stellar population. Assumed to be the same age as the galaxy and inputted by the user.
        n[int]:
            Number of particles in either the primary or secondary galaxy.
        
        Returns
        --------
        Spectral_Density[array]:
            The SED of a normalised stellar population at the same age as the galaxy. Extracted from the Bruzual and Charlot model.
        '''
                
        # Wavelengths = SEDs[1:,0]
        Ages = SEDs[0,1:]/1e9
        Fluxes = SEDs[1:,1:]
        
        Age_Index = np.where(Age >= Ages)[0][-1]
        
        Spectral_Density = Fluxes[:,Age_Index - 1].reshape(Fluxes.shape[0],1)*np.ones(n)
        Spectral_Density = Spectral_Density.transpose()
        
        return Spectral_Density
    
    def SFR_Flux(Spec_Ages,Spectra,Age,n,time,Mass,h):
        '''
        Calculate the total flux of a population formed by star formation at a specific age at some point in the simulation. New population is assumed to be a simple stellar population of the mass of that formed at each timestep. Age is t into the simulation. 

        Parameters
        -----------
        Spec_Ages[array]:
            Ages of model SEDs as extracted from Bruzual and Charlot. Used in comparison to the age of the newly formed population due to star formation. 
        Spectra[array]:
            Extracted SEDs from Bruzual and Charlot.
        Age[array]:
            Ages of each new population formed from star formation. 
        time[float]:
            Total time in the simulation.
        Mass[array]:
            The mass (in Solar Units) formed at each timestep at each particle in the simulation.
        h[float]:
            Timestep length.

        Returns
        --------
        SFR_Flux[array]:
            The modified SED due to star formation. This is all new flux due to star formation and will be added onto the SED of the initial simple stellar population. 

        '''
        # The New
        new_pop_age = Age + (time * h)
        age_index = np.where(new_pop_age <= Spec_Ages)[0][0]

        sfr_spectra = Spectra[:, int(age_index)]
        sfr_spectra = sfr_spectra[np.newaxis, :]

        total_mass_formed = np.sum(Mass, axis = 1)
        total_mass_formed = total_mass_formed[:, np.newaxis]
            
        sfr_flux = sfr_spectra * total_mass_formed
                
        return sfr_flux