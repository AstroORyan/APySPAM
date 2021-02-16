# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:19:25 2020

@author: oryan
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import make_lupton_rgb
from astropy.convolution import Gaussian2DKernel, convolve
from scipy.interpolate import interp2d


class Plotting:
    def Plotting_Function(Population_Colours_Array, Population_Flux_Array,Weight,position_integer, Particle_Coords, image_scale, Name,
                        n1, n2, Densities, SFRs, Time, TU, DU, r1, r2, Sec_Conversion_Frame,directory_results,Counts):
        
        # LSST_Area = 35e4                # Numbers pulled from the web, my dude
        # SDSS_Area = np.pi*(2.5)**2
        
        # Population_Flux_Array *= SDSS_Area
        
        fluxes, maps = Plotting.Image_Array(Population_Flux_Array,Particle_Coords,Weight,image_scale,DU,Densities,SFRs)
                                  
        maggies = Plotting.Colour_Conversion(fluxes)
        
        Gas_Map = maps[0]
        SFR_Map = maps[1]
        Weight_Map = maps[2]

        Gas_Map = np.transpose(Gas_Map)
        SFR_Map = np.transpose(SFR_Map)
        Weight_Map = np.transpose(Weight_Map)
                
        Gas_Map[Gas_Map == 0] = np.nan
        SFR_Map[SFR_Map == 0] = np.nan
        Weight_Map[Weight_Map == 0] = np.nan
        
        maggies[0] = np.transpose(maggies[0])
        maggies[1] = np.transpose(maggies[1])
        maggies[2] = np.transpose(maggies[2])
        maggies[3] = np.transpose(maggies[3])
        maggies[4] = np.transpose(maggies[4])
        
        maps = [Gas_Map, SFR_Map]
        
        Galaxy = 'Overall'
        Plotting.Colour_Plotter(maggies, directory_results, Galaxy, Name, Time, TU)
        
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
        fig.savefig(directory_results+str(Name)+r'\SFR_Density\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig) 
    
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pos = ax.imshow(Gas_Map,cmap=plt.cm.plasma,interpolation='nearest',origin='lower')#,vmin=0,vmax=25)
        ax.set_xlabel('X from Primary Galactic Center (MPc)')
        ax.set_ylabel('Y from Primary Galactic Center(MPc)')
        plt.title('Initial Gas Density Distribution in Arp 256')
        cbar = fig.colorbar(pos)
        cbar.set_label('Gas Density M$_{\odot}$pc$^{-2}$',rotation=270,labelpad=15)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.savefig(directory_results+str(Name)+r'\Gas_Density\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pos = ax.imshow(Weight_Map,cmap=plt.cm.plasma,interpolation='nearest',origin='lower')#,vmin=0,vmax=25)
        ax.set_xlabel('X from Primary Galactic Center (MPc)')
        ax.set_ylabel('Y from Primary Galactic Center(MPc)')
        plt.title('Total Gas Density in '+str(Name)+' at t = '+ '{:.3f}'.format(round(13.7 + Time*TU,3)))
        cbar = fig.colorbar(pos)
        cbar.set_label('Total Bin Weight',rotation=270,labelpad=15)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.savefig(directory_results+str(Name)+r'\Gas_Density\Weight_Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig)
        
        Plotting.Lupton_Plot(fluxes,Name,Time,TU,directory_results)
        
        Galaxy = 'Primary'
        Plotting.Plotting_Individ_Galaxy(image_scale,Sec_Conversion_Frame,position_integer,r1,DU,maggies,maps,Name,Time,TU,Galaxy,directory_results)
        
        Galaxy = 'Secondary'
        Plotting.Plotting_Individ_Galaxy(image_scale,Sec_Conversion_Frame,position_integer,r2,DU,maggies,maps,Name,Time,TU,Galaxy,directory_results)

    def Image_Array(Population_Flux_Array,Particle_Coords,Weight,image_scale,DU,Densities,SFRs):
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
            
            Gas_Map[x_index,y_index] += Densities[i]
            SFR_Map[x_index,y_index] += SFRs[i]
            Weight_Dist[x_index,y_index] += Weight[i]
                          
        Gas_Map /= (1e3)**2
            
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

    def Plotting_Individ_Galaxy(image_scale,Sec_Conversion_Frame,position_integer,r,DU,maggies,maps,Name,Time,TU,Galaxy,directory_results):  
        u_AB,g_AB,r_AB,i_AB,z_AB = maggies
        Gas_Map, SFR_Map = maps
        
        x_space = image_scale    #We've set this up in galaxy units!
        y_space = image_scale
        
        Resolution = np.ceil(image_scale/10)
        
        x_bins = int(x_space*DU/Resolution)
        y_bins = int(y_space*DU/Resolution)
        
        x_bin_width = x_space/x_bins
        y_bin_width = y_space/y_bins
        
        Simulation_Center_x = np.floor(x_bins/2) - 1
        Simulation_Center_y = np.floor(y_bins/2) - 1
        
        physical_step = x_bin_width*DU

        # Getting a map of just the Primary Galaxy
        if Galaxy == 'Primary':
            Galaxy_coords = [Simulation_Center_x, Simulation_Center_y]
        elif Galaxy == 'Secondary':
            Pert_x_Dist = np.ceil(Sec_Conversion_Frame[position_integer,0]*DU/physical_step)
            Pert_y_Dist = np.ceil(Sec_Conversion_Frame[position_integer,1]*DU/physical_step)
            Galaxy_coords = [Simulation_Center_x + Pert_x_Dist, Simulation_Center_y + Pert_y_Dist]
            
        
        Gal_r_Bins = np.ceil((1*r*DU)/physical_step)
        
        x_min = int(Galaxy_coords[0] - Gal_r_Bins)
        x_max = int(Galaxy_coords[0] + Gal_r_Bins)
        y_min = int(Galaxy_coords[1] - Gal_r_Bins)
        y_max = int(Galaxy_coords[1] + Gal_r_Bins)
            
    
        Reduced_u = u_AB[y_min:y_max, x_min:x_max]
        Reduced_g = g_AB[y_min:y_max, x_min:x_max]
        Reduced_r = r_AB[y_min:y_max, x_min:x_max]
        Reduced_i = i_AB[y_min:y_max, x_min:x_max]
        Reduced_z = z_AB[y_min:y_max, x_min:x_max]
        
        SFR_Reduced_Map = SFR_Map[y_min:y_max,x_min:x_max]
        Gas_Reduced_Map = Gas_Map[y_min:y_max,x_min:x_max]
        
        
        maggies = [Reduced_u, Reduced_g, Reduced_r, Reduced_i, Reduced_z]
        
        Plotting.Colour_Plotter(maggies,directory_results,Galaxy,Name,Time,TU)
        
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pos = ax.imshow(SFR_Reduced_Map,cmap=plt.cm.plasma,interpolation='nearest',origin='lower',vmin = 0, vmax = np.nanmax(SFR_Reduced_Map)/2)
        ax.set_xlabel('X from Primary Galactic Center (MPc)')
        ax.set_ylabel('Y from Primary Galactic Center(MPc)')
        plt.title('SFR Density in '+str(Name)+' Primary at t = '+ '{:.3f}'.format(round(13.7 + Time*TU,3))+' Resolution = '+str(Resolution)+'kpc')
        cbar = fig.colorbar(pos)
        cbar.set_label('M$_{\odot}$yr$^{-1}$',rotation=270,labelpad=15)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.savefig(directory_results+str(Name)+r'\SFR_'+str(Galaxy)+'\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig) 
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pos = ax.imshow(Gas_Reduced_Map,cmap=plt.cm.plasma,interpolation='nearest',origin='lower')#,vmin=0,vmax=25)
        ax.set_xlabel('X from Primary Galactic Center (MPc)')
        ax.set_ylabel('Y from Primary Galactic Center(MPc)')
        plt.title('Gas Density in '+str(Name)+' Primary at t = '+ '{:.3f}'.format(round(13.7 + Time*TU,3))+' Resolution = '+str(Resolution)+'kpc')
        cbar = fig.colorbar(pos)
        cbar.set_label('Gas Density (M$_{\odot}$pc$^{-2}$)',rotation=270,labelpad=15)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.savefig(directory_results+str(Name)+r'\Gas_Density_'+str(Galaxy)+'\Time_' + '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
        plt.close(fig) 
        
        
    def Colour_Plotter(maggies,directory_results, Galaxy, Name,Time,TU):
        Labels = ['u_filter\\','g_filter\\','r_filter\\','i_filter\\','z_filter\\']
        Titles = [Name + ' u Filter (SDSS)', Name + ' u Filter  (SDSS)', Name + ' g Filter  (SDSS)',
                  Name + ' i Filter  (SDSS)', Name + ' z Filter  (SDSS)']
        
        for i in range(5):
            plotting = maggies[i]
            plotting[plotting == 0] = np.nan
            fig = plt.figure()
            ax = fig.add_subplot(111)
            color_map = plt.cm.get_cmap('plasma')
            reversed_color_map = color_map.reversed()
            pos = ax.imshow(plotting,cmap=reversed_color_map,interpolation='nearest',origin='lower',vmin=np.nanmin(plotting),vmax=np.nanmax(plotting))
            plt.title(Titles[i],fontsize=14)
            cbar = fig.colorbar(pos)
            cbar.ax.tick_params(labelsize=10)
            cbar.set_label('Magnitude', rotation=270, labelpad=15)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            fig.savefig(directory_results+str(Name)+'\\Colour_Images\\'+str(Galaxy)+'\\'+Labels[i]+'Time_'+ '{:.3f}'.format(round(13.7 + Time*TU,3)) + '.png')
            plt.close(fig)