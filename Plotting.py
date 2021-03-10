# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:19:25 2020

@author: oryan
"""
import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import make_lupton_rgb
from astropy.convolution import Gaussian2DKernel, convolve
import os


class Plotting:
    def Plotting_Function(Population_Flux_Array,position_integer, Particle_Coords, image_scale, Name,
                        n1, n2, SFRs, Time, TU, DU, Sec_Conversion_Frame,directory_results,N_filters):
        
        Results_dir = directory_results+Name
        if os.path.isdir(Results_dir):
            pass
        else:
            os.mkdir(Results_dir)
        folders_make = ['\\SFR_Density','\\Colour_Images','\\Lupton_Images']
        for q in range(len(folders_make)):
            if os.path.isdir(Results_dir+folders_make[q]):
                pass
            else:
                os.mkdir(Results_dir+folders_make[q])
        
        fluxes, maps = Plotting.Image_Array(Population_Flux_Array,Particle_Coords,image_scale,DU,SFRs)
                                  
        maggies = Plotting.Colour_Conversion(fluxes)
        
        SFR_Map = maps[1]

        SFR_Map = np.transpose(SFR_Map)
                
        SFR_Map[SFR_Map == 0] = np.nan
        
        maggies[0] = np.transpose(maggies[0])
        maggies[1] = np.transpose(maggies[1])
        maggies[2] = np.transpose(maggies[2])
        maggies[3] = np.transpose(maggies[3])
        maggies[4] = np.transpose(maggies[4])
                
        Plotting.Colour_Plotter(maggies, Results_dir, Name, Name, Time, TU,N_filters)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pos = ax.imshow(SFR_Map,cmap=plt.cm.plasma,interpolation='nearest',origin='lower',vmin = 0, vmax = np.nanmax(SFR_Map)/2)
        ax.set_xlabel('X from Primary Galactic Center (MPc)')
        ax.set_ylabel('Y from Primary Galactic Center(MPc)')
        plt.title('SFR Density in '+str(Name)+' at t = '+ '{:.3f}'.format(round(13.7 + Time*TU,3)))
        cbar =  fig.colorbar(pos)
        cbar.set_label('SFR M$_{\odot}$yr$^{-1}$',rotation=270,labelpad=15)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.savefig(Results_dir+r'\SFR_Density\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig) 
        
        Plotting.Lupton_Plot(fluxes,Name,Time,TU,directory_results)

    def Image_Array(Population_Flux_Array,Particle_Coords,image_scale,DU,SFRs):
        u_Flux = Population_Flux_Array[:,0]
        g_Flux = Population_Flux_Array[:,1]
        r_Flux = Population_Flux_Array[:,2]
        i_Flux = Population_Flux_Array[:,3]
        z_Flux = Population_Flux_Array[:,4]
    
    #    Below lines create intermediate images for making movies of the simulation
        x = Particle_Coords[:,0]
        y = Particle_Coords[:,1]
        
    
        x_space = image_scale    #We've set this up in galaxy units!
        y_space = image_scale
        #We then set up a 1kpc resolution
        Resolution = np.ceil(image_scale/10)    #In kpc
        x_bins = int(x_space*DU/Resolution)
        y_bins = int(y_space*DU/Resolution)
        
        x_step = x_space/x_bins
        y_step = y_space/y_bins
        x_bin_width = x_space/x_bins
        y_bin_width = y_space/y_bins
        x_bin_step = x_step/x_bin_width
        y_bin_step = x_step/y_bin_width
    
        magnitude_array_u = np.zeros([x_bins,y_bins])    
        magnitude_array_g = np.zeros([x_bins,y_bins])    
        magnitude_array_r = np.zeros([x_bins,y_bins])    
        magnitude_array_i = np.zeros([x_bins,y_bins])    
        magnitude_array_z = np.zeros([x_bins,y_bins])
        Gas_Map = np.zeros([x_bins,y_bins])
        SFR_Map = np.zeros([x_bins,y_bins])
        Weight_Dist = np.zeros([x_bins,y_bins])
        
        for i in range(len(x) - 1):
            
            if (x[i] > (x_space/2)) or (x[i] < -(x_space/2)):
                continue
            if (y[i] > (y_space/2)) or (y[i] < -(y_space/2)):
                continue
            
            x_temp = x[i]
            y_temp = y[i]
            x_step_counter = -x_space/2
            y_step_counter = -y_space/2
            x_bin_step_counter = 0
            y_bin_step_counter = 0
            for j in range(x_bins):
                x_step_counter = x_step_counter + x_step
                x_bin_step_counter = x_bin_step_counter + x_bin_step
                if x_temp < x_step_counter:
                    x_index = int(np.floor(x_bin_step_counter - 1))
                    break
            for k in range(y_bins):
                    y_step_counter = y_step_counter + y_step
                    y_bin_step_counter = y_bin_step_counter + y_bin_step
                    if y_temp < y_step_counter:
                        y_index = int(np.floor(y_bin_step_counter - 1))
                        break
    
            magnitude_array_u[x_index,y_index] += u_Flux[i]
            magnitude_array_g[x_index,y_index] += g_Flux[i]
            magnitude_array_r[x_index,y_index] += r_Flux[i]
            magnitude_array_i[x_index,y_index] += i_Flux[i]
            magnitude_array_z[x_index,y_index] += z_Flux[i]
            
            SFR_Map[x_index,y_index] += SFRs[i]
                          
            
        return [magnitude_array_u, magnitude_array_g, magnitude_array_r, magnitude_array_i, magnitude_array_z], [Gas_Map, SFR_Map, Weight_Dist]
    
    def Colour_Conversion(fluxes):
        u_flux_array = fluxes[0]
        g_flux_array = fluxes[1]
        r_flux_array = fluxes[2]
        i_flux_array = fluxes[3]
        z_flux_array = fluxes[4]
        
        dims = u_flux_array.shape
        
        u_AB_Array = np.zeros(dims)    
        g_AB_Array = np.zeros(dims)
        r_AB_Array = np.zeros(dims)
        i_AB_Array = np.zeros(dims)
        z_AB_Array = np.zeros(dims)
        
        u_constant = -2.5*np.log10(3631e-23)           #These Lambda zero points were found using: http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
        g_constant = -2.5*np.log10(3631e-23)
        r_constant = -2.5*np.log10(3631e-23)
        i_constant = -2.5*np.log10(3631e-23)
        z_constant = -2.5*np.log10(3631e-23)

        u_AB_Array[u_flux_array > 0] = -2.5*np.log10(u_flux_array[u_flux_array > 0]) - u_constant
        g_AB_Array[g_flux_array > 0] = -2.5*np.log10(g_flux_array[g_flux_array > 0]) - g_constant
        r_AB_Array[r_flux_array > 0] = -2.5*np.log10(r_flux_array[r_flux_array > 0]) - r_constant
        i_AB_Array[i_flux_array > 0] = -2.5*np.log10(i_flux_array[i_flux_array > 0]) - i_constant
        z_AB_Array[z_flux_array > 0] = -2.5*np.log10(z_flux_array[z_flux_array > 0]) - z_constant

        return [u_AB_Array, g_AB_Array, r_AB_Array, i_AB_Array, z_AB_Array]

    def Lupton_Plot(maggies,Name,Time,TU,directory_results):        
        
        test_g = maggies[1]
        test_r = maggies[2]
        test_i = maggies[3]
        
        # Note, the scalings here are those given in SDSS scripts
        scales = np.asarray([45.0,22.0,28.0])#np.asarray([4.9,5.7,7.8]) #np.asarray([5.0,5.0,5.0])#
        scales *= 1.75*1e27
        
        fixed_image_g_filtered = test_g*scales[2]
        fixed_image_r_filtered = test_r*scales[1]
        fixed_image_i_filtered = test_i*scales[0]
        
        r = np.transpose(fixed_image_i_filtered)
        g = np.transpose(fixed_image_r_filtered)
        b = np.transpose(fixed_image_g_filtered)
        
        Noise_Array = 1e-28*np.random.poisson(lam=5,size=[fixed_image_g_filtered.shape[0],fixed_image_g_filtered.shape[1],3])
        
        r += Noise_Array[:,:,0]
        g += Noise_Array[:,:,1]
        b += Noise_Array[:,:,2]
        
        kernel = Gaussian2DKernel(x_stddev=1,y_stddev=1) # NOTE: Array size is 8*x(y)_stddev
        r = convolve(r,kernel)
        g = convolve(g,kernel)
        b = convolve(b,kernel)
        
        rgb_default = make_lupton_rgb(r,g,b,minimum=0, Q=10,stretch=0.05)
        plt.figure()
        plt.imshow(rgb_default,origin='lower')
        plt.title('Image of Arp 256')
        plt.savefig(directory_results+str(Name)+r'\Lupton_Images\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close()        
        
    def Colour_Plotter(maggies,directory_results, Galaxy, Name,Time,TU,N):
        Labels = []
        
        for i in range(N):
            Labels.append('Filter_'+ str(i))
        
        
        for i in range(N):
            plotting = maggies[i]
            plotting[plotting == 0] = np.nan
            fig = plt.figure()
            ax = fig.add_subplot(111)
            color_map = plt.cm.get_cmap('plasma')
            reversed_color_map = color_map.reversed()
            pos = ax.imshow(plotting,cmap=reversed_color_map,interpolation='nearest',origin='lower',vmin=np.nanmin(plotting),vmax=np.nanmax(plotting))
            plt.title(Labels[i],fontsize=14)
            cbar = fig.colorbar(pos)
            cbar.ax.tick_params(labelsize=10)
            cbar.set_label('Magnitude', rotation=270, labelpad=15)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            fig.savefig(directory_results +'\\Colour_Images\\'+Labels[i]+'_Time_'+ '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
            plt.close(fig)