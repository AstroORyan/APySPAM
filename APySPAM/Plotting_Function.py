# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 15:22:15 2021

@author: oryan

Quickly creates an image of the system using imshow and a meshgrid.

"""
import matplotlib.pyplot as plt
import numpy as np
import uuid

class Plotting_Function:
    def plotting(Coordinates,part_flux,SFRs,n_filters):
        '''
        Creates a white light image of the interacting galaxy system. With movie_flag activated, this will be run every 10 time steps to create a movie of interaction. With it false, this will just be called at the end.

        Parameters
        ------------
        Coordinates:
            N_particles x 6 array of particle coordinates and velocities. In simulation units.
        part_flux:
            N_particles x N_filters array of integrated flux for each particle.
        SFRs:
            N_particles x 1 array of SFRs at of each particle. Can be used to create a SFR map of the interacting galaxies.
        n_filters:
            Integer, number of filters input into simulation.

        Returns
        --------
        None
        '''

        # First, extract the wanted dimensions. 
        x = Coordinates[:,0]
        y = Coordinates[:,1]
        
        total_flux = np.zeros(len(x) - 1)
        
        for i in range(n_filters):
            total_flux += part_flux[i]
        
        
        # Define Constants
        DU = 15
        Resolution = 1.5/DU    # in kpc
        
        # Define size of observed space
        Image = np.zeros([50,50])
        x_min = (-Image.shape[0]/2)*Resolution
        x_max = (Image.shape[0]/2)*Resolution
        
        y_min = (-Image.shape[1]/2)*Resolution
        y_max = (Image.shape[1]/2)*Resolution
        
        x_pixel_value = np.linspace(x_min,x_max,Image.shape[0])
        y_pixel_value = np.linspace(y_min,y_max,Image.shape[1])
        
        for i in range(len(x) - 1):
            if x[i] > x_max or x[i] < x_min:
                continue
            elif y[i] > y_max or y[i] < y_min:
                continue
            else:
                p = np.where(x[i] >= x_pixel_value)[0][-1]
                q = np.where(y[i] >= y_pixel_value)[0][-1]

                Image[p,q] += total_flux[i]
            
        plt.figure()
        plt.imshow(Image.T, origin='lower')
        plt.title('White Image')
        
        plt.figure()
        plt.imshow(-2.5*np.log10(Image.T) - 48.6, origin='lower')
        plt.title('AB Magnitude Map')