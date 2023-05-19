# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 15:22:15 2021

@author: oryan

Quickly creates an image of the system using imshow and a meshgrid.

"""
import matplotlib.pyplot as plt
import numpy as np
import sys

from astropy.cosmology import FlatLambdaCDM
from astropy.nddata import block_reduce
from astropy.convolution import Gaussian2DKernel, convolve, interpolate_replace_nans
import astropy.units as u

# def interpolating(im):
#     mask = np.load(r'C:\Users\oryan\Documents\PySPAM_Original_Python_MCMC_Counts\masks\Arp240-mask-50.npy')
#     pix = np.argwhere(im == 0)
#     for i in pix:
#         x = i[0]
#         y = i[1]
#         if mask[x,y] and np.sum(im[x-1:x+1, y-1:y+1] == 0) < 2:
#             im[x,y] = np.sum(im[x-1:x+1, y-1:y+1]) / (im[x-1:x+1, y-1:y+1].shape[1] * im[x-1:x+1, y-1:y+1].shape[0])

#     return im

# def apply_psf(im):

#     psf = np.load('C:/Users/oryan/Documents/PySPAM_Original_Python/APySPAM/psfs/full-psf.npy')

#     convolved_image = convolve_fft(im, psf, boundary = 'fill')
#     return convolved_image

def assign_pixel(flux, x_pix, y_pix, im, part_coord, x_pix_coord, y_pix_coord, lim):
    
    sep = []
    pix_coord = [[x_pix_coord[x_pix]], [y_pix_coord[y_pix]]]
    sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
    # im[x_pix, y_pix] += sep*flux
    
    pix_coord = [[x_pix_coord[x_pix-1]], [y_pix_coord[y_pix]]]
    sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
    # im[x_pix-1, y_pix] += sep*flux
    
    pix_coord = [[x_pix_coord[x_pix]], [y_pix_coord[y_pix-1]]]
    sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
    # im[x_pix, y_pix-1] += sep*flux
    
    pix_coord = [[x_pix_coord[x_pix-1]], [y_pix_coord[y_pix-1]]]
    sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
    # im[x_pix-1, y_pix-1] += sep*flux
    
    if not x_pix + 1 > lim[0]:
        pix_coord = [[x_pix_coord[x_pix+1]], [y_pix_coord[y_pix]]]
        sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
        # im[x_pix+1, y_pix] += sep*flux
    
    if not y_pix + 1 > lim[0]:
        pix_coord = [[x_pix_coord[x_pix]], [y_pix_coord[y_pix+1]]]
        sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
        # im[x_pix, y_pix+1] += sep*flux
    
    if not x_pix + 1 > lim[0] and not y_pix + 1 > lim[1]:
        pix_coord = [[x_pix_coord[x_pix+1]], [y_pix_coord[y_pix+1]]]
        sep.append(np.sqrt((pix_coord[0] - part_coord[0])**2 + (pix_coord[1] - part_coord[1])**2))
        # im[x_pix+1, y_pix+1] += sep*flux
    
    sep = sep / np.sum(sep)
    
    im[x_pix, y_pix] += sep[0]*flux
    im[x_pix-1, y_pix] += sep[1]*flux
    im[x_pix, y_pix-1] += sep[2]*flux
    im[x_pix-1, y_pix-1] += sep[3]*flux
    if not x_pix + 1 > lim[0]:
        im[x_pix+1, y_pix] += sep[4]*flux
    if not y_pix + 1 > lim[0]:
        im[x_pix, y_pix+1] += sep[5]*flux
    if not x_pix + 1 > lim[0] and not y_pix + 1 > lim[1]:
        im[x_pix+1, y_pix+1] += sep[6]*flux
    
    return im  

def counts_convert(column_si, ct_nmgy):
    jansky_nmgy = (3.631e-6 * u.Jy / u.nanomaggy)
    si_jansky = (1e-23 / u.Jy)

    column_jy = column_si / si_jansky
    column_nmgy = column_jy / jansky_nmgy
    column_ct = column_nmgy / (ct_nmgy * u.nmgy / u.count)

    assert column_ct[0].unit == u.count

    return np.asarray(column_ct / u.count)


class Plotting_Function:
    def binary_plotting(Coordinates, redshift, native_res, reduction, iteration):
        # First, extract the wanted dimensions.  
        x = Coordinates[:,0]
        y = Coordinates[:,1]
                
        # Define Constants
        DU = 15
        cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.275 * u.K, Om0 = 0.308)
        conversion = cosmo.kpc_proper_per_arcmin(redshift)
        #Resolution = float((native_res * u.arcsec)  * conversion.to(u.kpc / u.arcsec) * (block_reduce) / (15 * u.kpc))

        Resolution = float((native_res * u.arcsec)  * conversion.to(u.kpc / u.arcsec) / (15 * u.kpc))
        im_size = int(reduction *  50)
        
        # Define size of observed space
        Image = np.zeros([im_size,im_size])
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

                Image[p,q] = 1

        red_image = block_reduce(Image, reduction, func = np.sum)
            
        plt.figure(figsize=(12,8))
        plt.imshow(Image, origin = 'lower')
        plt.colorbar()
        plt.title('AB Magnitude Map')
        plt.savefig(f'C:/Users/oryan/Documents/PySPAM_Original_Python/APySPAM/im-fold/{iteration}.jpeg')
        plt.close()


    def plotting(Coordinates,part_flux,n_filters,redshift, native_res, reduction,flag,convs):
        '''
        Creates a white light image of the interacting galaxy system. With movie_flag activated, this will be run every 10 time steps to create a movie of interaction. With it false, this will just be called at the end.

        Parameters
        ------------
        Coordinates:
            N_particles x 6 array of particle coordinates and velocities. In simulation units.
        part_flux:
            N_particles x N_filters array of integrated flux for each particle.
        n_filters:
            Integer, number of filters input into simulation.

        Returns
        --------
        None
        '''
        # First, extract the wanted dimensions. 
        x = Coordinates[:,0]
        y = Coordinates[:,1]
        
        total_counts = np.zeros(len(x) - 1)
        
        for i in list(part_flux.keys()):
            total_counts += counts_convert(part_flux[i],convs[i])        
        
        # Define Constants
        DU = 15
        cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.275 * u.K, Om0 = 0.308)
        conversion = cosmo.kpc_proper_per_arcmin(redshift)
        # Resolution = float((native_res * u.arcsec)  * conversion.to(u.kpc / u.arcsec) * (block_reduce) / (15 * u.kpc))

        Resolution = float((native_res * u.arcsec)  * conversion.to(u.kpc / u.arcsec) / (15 * u.kpc))

        im_size = int(reduction *  50)
        
        # Define size of observed space
        Image = np.zeros([im_size,im_size])
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

                Image = assign_pixel(total_counts[i], p, q, Image, [x[i],y[i]], x_pixel_value, y_pixel_value, Image.shape)

        # Image_psf = apply_psf(Image)

        red_image = block_reduce(Image, reduction, func = np.sum)
        # im_exp = interpolating(red_image)

        if flag: 
           # plt.figure(figsize=(12,8))
           # plt.imshow(Image.T, origin='lower')
           # plt.title('White Image')
            
            plt.figure(figsize=(12,8))
            plt.imshow(red_image, origin = 'lower')
            plt.colorbar()
            plt.title('AB Magnitude Map')
            plt.savefig('C:/Users/oryan/Documents/PySPAM_Original_Python/APySPAM/im-fold/tmp1.jpeg')
            plt.close()

        return red_image