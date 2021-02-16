# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:25:20 2020

@author: oryan
"""

import numpy as np

class Colours:
    def Colour_Calculation(Spectral_Density, Wavelength, redshift, n, d_cm, directory):
      ##Setting up the Filter Data
#      LSST_u_filename = directory+r'\Colours\LSST_Filter_Sets_AB\LSST_LSST.u.dat'
#      LSST_g_filename = directory+r'\Colours\LSST_Filter_Sets_AB\LSST_LSST.g.dat'
#      LSST_r_filename = directory+r'\Colours\LSST_Filter_Sets_AB\LSST_LSST.r.dat'
#      LSST_i_filename = directory+r'\Colours\LSST_Filter_Sets_AB\LSST_LSST.i.dat'
#      LSST_z_filename = directory+r'\Colours\LSST_Filter_Sets_AB\LSST_LSST.z.dat'
      
      SDSS_u_filename = directory+r'\Colours\SDSS_Filter_Sets_AB\SLOAN_SDSS.u.dat'
      SDSS_g_filename = directory+r'\Colours\SDSS_Filter_Sets_AB\SLOAN_SDSS.g.dat'
      SDSS_r_filename = directory+r'\Colours\SDSS_Filter_Sets_AB\SLOAN_SDSS.r.dat'
      SDSS_i_filename = directory+r'\Colours\SDSS_Filter_Sets_AB\SLOAN_SDSS.i.dat'
      SDSS_z_filename = directory+r'\Colours\SDSS_Filter_Sets_AB\SLOAN_SDSS.z.dat'
#      
      LSST_u_filter_Array = np.loadtxt(SDSS_u_filename)
      LSST_g_filter_Array = np.loadtxt(SDSS_g_filename)
      LSST_r_filter_Array = np.loadtxt(SDSS_r_filename)
      LSST_i_filter_Array = np.loadtxt(SDSS_i_filename)
      LSST_z_filter_Array = np.loadtxt(SDSS_z_filename)
      
      Wavelength_u = LSST_u_filter_Array[:,0]
      Wavelength_g = LSST_g_filter_Array[:,0]
      Wavelength_r = LSST_r_filter_Array[:,0]
      Wavelength_i = LSST_i_filter_Array[:,0]
      Wavelength_z = LSST_z_filter_Array[:,0]
      
      u_Transmission = LSST_u_filter_Array[:,1]
      g_Transmission = LSST_g_filter_Array[:,1]
      r_Transmission = LSST_r_filter_Array[:,1]
      i_Transmission = LSST_i_filter_Array[:,1]
      z_Transmission = LSST_z_filter_Array[:,1]
      
      Wavelength_u_lr = np.zeros(int(np.ceil(len(Wavelength_u)+100)))
      Wavelength_g_lr = np.zeros(int(np.ceil(len(Wavelength_g)+100)))
      Wavelength_r_lr = np.zeros(int(np.ceil(len(Wavelength_r)+100)))
      Wavelength_i_lr = np.zeros(int(np.ceil(len(Wavelength_i)+100)))
      Wavelength_z_lr = np.zeros(int(np.ceil(len(Wavelength_z)+100)))
      
      u_Transmission_lr = np.zeros(int(np.ceil(len(Wavelength_u)+100)))
      g_Transmission_lr = np.zeros(int(np.ceil(len(Wavelength_g)+100)))
      r_Transmission_lr = np.zeros(int(np.ceil(len(Wavelength_r)+100)))
      i_Transmission_lr = np.zeros(int(np.ceil(len(Wavelength_i)+100)))
      z_Transmission_lr = np.zeros(int(np.ceil(len(Wavelength_z)+100)))
      
      # Will need to add a few lines here which will apply redshift to the measurements. This will be done later.
      Wavelength_Emitted = Wavelength
      Wavelength_Observed = Wavelength_Emitted*(1+redshift)
      
      #Degrades the LSST filter data to be bins of 2 Angstrom.
      counter_lr = 0
      counter_array = 0
      for i in range(len(Wavelength)):
          if Wavelength_Observed[i] >= Wavelength_u[counter_array]:
              Wavelength_u_lr[counter_lr] = Wavelength_Observed[i]
              difference = int(np.ceil(Wavelength_Observed[i] - Wavelength_u[counter_array]))
              if counter_array == 0:
                  u_Padding_Index = i
              counter_array = counter_array + difference
              if counter_array >= Wavelength_u[-1] - Wavelength_u[0]:
                  break
              u_Transmission_lr[counter_lr] = u_Transmission[counter_array]
              counter_lr = counter_lr + 1
      
      u_Transmission_lr = u_Transmission_lr[0:counter_lr]
      
      counter_lr = 0
      counter_array = 0
      for i in range(len(Wavelength_Observed)):
          if Wavelength_Observed[i] >= Wavelength_g[counter_array]:
              Wavelength_g_lr[counter_lr] = Wavelength[i]
              difference = int(np.ceil(Wavelength_Observed[i] - Wavelength_g[counter_array]))
              if counter_array == 0:
                  g_Padding_Index = i
              counter_array = counter_array + difference
              if counter_array >= Wavelength_g[-1] - Wavelength_g[0]:
                  break
              g_Transmission_lr[counter_lr] = g_Transmission[counter_array]
              counter_lr = counter_lr + 1
      
      g_Transmission_lr = g_Transmission_lr[0:counter_lr]
              
              
      counter_lr = 0
      counter_array = 0
      for i in range(len(Wavelength_Observed)):
          if Wavelength_Observed[i] >= Wavelength_r[counter_array]:
              Wavelength_r_lr[counter_lr] = Wavelength_Observed[i]
              difference = int(np.ceil(Wavelength_Observed[i] - Wavelength_r[counter_array]))
              if counter_array == 0:
                  r_Padding_Index = i
              counter_array = counter_array + difference
              if counter_array >= Wavelength_r[-1] - Wavelength_r[0]:
                  break
              r_Transmission_lr[counter_lr] = r_Transmission[counter_array]
              counter_lr = counter_lr + 1
              
      r_Transmission_lr = r_Transmission_lr[0:counter_lr]

      
      counter_lr = 0
      counter_array = 0
      for i in range(len(Wavelength_Observed)):
          if Wavelength_Observed[i] >= Wavelength_i[counter_array]:
              Wavelength_i_lr[counter_lr] = Wavelength_Observed[i]
              difference = int(np.ceil(Wavelength_Observed[i] - Wavelength_i[counter_array]))
              if counter_array == 0:
                  i_Padding_Index = i
              counter_array = counter_array + difference
              if counter_array >= Wavelength_i[-1] - Wavelength_i[0]:
                  break
              i_Transmission_lr[counter_lr] = i_Transmission[counter_array]
              counter_lr = counter_lr + 1
      
      i_Transmission_lr = i_Transmission_lr[0:counter_lr]
              
      counter_lr = 0
      counter_array = 0
      for i in range(len(Wavelength_Observed)):
          if Wavelength_Observed[i] >= Wavelength_z[counter_array]:
              Wavelength_z_lr[counter_lr] = Wavelength_Observed[i]
              difference = int(np.ceil(Wavelength_Observed[i] - Wavelength_z[counter_array]))
              if counter_array == 0:
                  z_Padding_Index = i
              counter_array = counter_array + difference
              if counter_array >= Wavelength_z[-1] - Wavelength_z[0]:
                  break
              z_Transmission_lr[counter_lr] = z_Transmission[counter_array]
              counter_lr = counter_lr + 1
      
      z_Transmission_lr = z_Transmission_lr[0:counter_lr]
      
              
      #Now need to pad data out so that filter data is the same size as Bruzual and Charlot data
      required_length = Spectral_Density.shape[1]
      Transmission_u_Padded = np.zeros(required_length)
      Transmission_g_Padded = np.zeros(required_length)
      Transmission_r_Padded = np.zeros(required_length)
      Transmission_i_Padded = np.zeros(required_length)
      Transmission_z_Padded = np.zeros(required_length)

      #Now will align the transmission curve with the wavelengths.      
      Transmission_u_Padded[u_Padding_Index:u_Padding_Index + len(u_Transmission_lr)] = u_Transmission_lr
      Transmission_g_Padded[g_Padding_Index:g_Padding_Index + len(g_Transmission_lr)] = g_Transmission_lr
      Transmission_r_Padded[r_Padding_Index:r_Padding_Index + len(r_Transmission_lr)] = r_Transmission_lr
      Transmission_i_Padded[i_Padding_Index:i_Padding_Index + len(i_Transmission_lr)] = i_Transmission_lr
      Transmission_z_Padded[z_Padding_Index:z_Padding_Index + len(z_Transmission_lr)] = z_Transmission_lr

      #Quickly Apply an extinction algorithm
      Coefficient = np.zeros(len(Wavelength))
      test_wav = Wavelength_Observed*1e-10/1e-6
      E_B_V = 0.44
      for i in range(len(Wavelength)):              # Equations are from arxiv astro-ph/0109035
          if test_wav[i] < 0.12:
              Coefficient[i] = 0
          elif (test_wav[i] >= 0.12 and test_wav[i] <= 0.63):
              Coefficient[i] = 1.17*(-2.156 + 1.509/test_wav[i] - 0.198/(test_wav[i]**2) + 0.011/(test_wav[i]**3)) + 1.78
          elif (test_wav[i] > 0.63 and test_wav[i] <= 2.2):
              Coefficient[i] = 1.17*(-1.857 + 1.040/test_wav[i]) + 1.78
          else:
              Coefficient[i] = 0
        
      Spectral_Density = Spectral_Density*10**(-0.4*Coefficient*E_B_V)
      
      # Now,we convolve the transmission curves with every population spectral energy distribution
      dims = [n,required_length]
      Population_Colour_Flux_u = np.zeros(dims)
      Population_Colour_Flux_g = np.zeros(dims)
      Population_Colour_Flux_r = np.zeros(dims)
      Population_Colour_Flux_i = np.zeros(dims)
      Population_Colour_Flux_z = np.zeros(dims)
      
      Frequency_Emitted = 2.988e8/(1e-10*Wavelength_Emitted)
      Frequency_Observed = 2.998e8/(1e-10*Wavelength_Observed)
      Frequency_Emitted = np.flip(Frequency_Emitted)
      Frequency_Observed = np.flip(Frequency_Observed)
      
      for i in range(n):
          Population_Colour_Flux_u[i,:] = (((Spectral_Density[i,:]*Transmission_u_Padded)*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
          Population_Colour_Flux_g[i,:] = (((Spectral_Density[i,:]*Transmission_g_Padded)*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
          Population_Colour_Flux_r[i,:] = (((Spectral_Density[i,:]*Transmission_r_Padded)*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
          Population_Colour_Flux_i[i,:] = (((Spectral_Density[i,:]*Transmission_i_Padded)*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
          Population_Colour_Flux_z[i,:] = (((Spectral_Density[i,:]*Transmission_z_Padded)*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
        
      Population_Colour_Flux_Integrated_u = np.zeros(n)
      Population_Colour_Flux_Integrated_g = np.zeros(n)
      Population_Colour_Flux_Integrated_r = np.zeros(n)
      Population_Colour_Flux_Integrated_i = np.zeros(n)
      Population_Colour_Flux_Integrated_z = np.zeros(n)
      
      Convert_Units = 3.826e33  #Flux is in units of Solar Luminosity per Hertz. Have removed the Hertzs from integration, now to do this.
      
      h = 6.63e-34
      
      Counts_g = np.sum(Population_Colour_Flux_g/(h*Frequency_Observed), 1)
      Counts_r = np.sum(Population_Colour_Flux_r/(h*Frequency_Observed), 1)
      Counts_i = np.sum(Population_Colour_Flux_i/(h*Frequency_Observed), 1)
      
      Counts_g *= (Convert_Units/(4*np.pi*d_cm**2))
      Counts_r *= (Convert_Units/(4*np.pi*d_cm**2))
      Counts_i *= (Convert_Units/(4*np.pi*d_cm**2))
      
      Counts = [Counts_g, Counts_r, Counts_i]
         
      for i in range(n):    # This for loop calibrates the found spectra to the distance to the observed galaxy.
          Population_Colour_Flux_Integrated_u[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Colour_Flux_u[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
          Population_Colour_Flux_Integrated_g[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Colour_Flux_g[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
          Population_Colour_Flux_Integrated_r[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Colour_Flux_r[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
          Population_Colour_Flux_Integrated_i[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Colour_Flux_i[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
          Population_Colour_Flux_Integrated_z[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Colour_Flux_z[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
      
      u_Transmission_Integrated = np.trapz(Transmission_u_Padded/Frequency_Emitted, Frequency_Emitted)
      g_Transmission_Integrated = np.trapz(Transmission_g_Padded/Frequency_Emitted, Frequency_Emitted)
      r_Transmission_Integrated = np.trapz(Transmission_r_Padded/Frequency_Emitted, Frequency_Emitted)
      i_Transmission_Integrated = np.trapz(Transmission_i_Padded/Frequency_Emitted, Frequency_Emitted)
      z_Transmission_Integrated = np.trapz(Transmission_z_Padded/Frequency_Emitted, Frequency_Emitted)
      
      Final_Population_Colour_u = Population_Colour_Flux_Integrated_u/u_Transmission_Integrated
      Final_Population_Colour_g = Population_Colour_Flux_Integrated_g/g_Transmission_Integrated
      Final_Population_Colour_r = Population_Colour_Flux_Integrated_r/r_Transmission_Integrated
      Final_Population_Colour_i = Population_Colour_Flux_Integrated_i/i_Transmission_Integrated
      Final_Population_Colour_z = Population_Colour_Flux_Integrated_z/z_Transmission_Integrated
      
      u_AB = -2.5*np.log10(Final_Population_Colour_u) - 48.6
      g_AB = -2.5*np.log10(Final_Population_Colour_g) - 48.6
      r_AB = -2.5*np.log10(Final_Population_Colour_r) - 48.6
      i_AB = -2.5*np.log10(Final_Population_Colour_i) - 48.6
      z_AB = -2.5*np.log10(Final_Population_Colour_z) - 48.6
      
      Population_Colours_Array = np.asarray([u_AB,g_AB,r_AB,i_AB,z_AB])
      Population_Colours_Array = np.transpose(Population_Colours_Array)
      
      Population_Flux_Array = np.asarray([Final_Population_Colour_u,Final_Population_Colour_g,Final_Population_Colour_r,Final_Population_Colour_i,Final_Population_Colour_z])
      Population_Flux_Array = np.transpose(Population_Flux_Array)

      return Population_Colours_Array, Population_Flux_Array, Counts
