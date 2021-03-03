# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:25:20 2020

@author: oryan
"""

import numpy as np

class Colours:
    def Colour_Calculation(filter_data,Spectral_Density, Wavelength_SED, redshift, n, d_cm):
      ##Setting up the Filter Data
      Wavelengths = []
      Transmission = []
      
      N_filters = len(filter_data)
      filter_array_length = Spectral_Density.shape[1]
      
      for i in range(N_filters):
          tmp = filter_data[i]
          Wavelengths.append(tmp[:,0])
          Transmission.append(tmp[:,1])
      del tmp
      
      Wavelength_lr = []
      Transmission_lr = []
      
      for i in range(N_filters):
          tmp = np.zeros(len(Wavelengths[i])+100)
          Wavelength_lr.append(tmp)
          tmp = np.zeros(len(Transmission[i]) + 100)
          Transmission_lr.append(tmp)
      del tmp
      
      Wavelength_Emitted = Wavelength_SED
      Wavelength_Observed = Wavelength_Emitted*(1+redshift)
      
      #Degrades the LSST filter data to be bins of 2 Angstrom.
      for h in range(N_filters):
          counter_lr = 0
          counter_array = 0
          Wavelength_tmp = Wavelengths[h]
          Wavelength_lr_tmp = Wavelength_lr[h]
          Transmission_tmp = Transmission[h]
          Transmission_lr_tmp = Transmission_lr[h]
          Transmission_Padded = np.zeros(filter_array_length)
          for i in range(len(Wavelength_SED)):
              if Wavelength_Observed[i] >= Wavelength_tmp[counter_array]:
                  Wavelength_lr_tmp[counter_lr] = Wavelength_Observed[i]
                  difference = int(np.floor(Wavelength_Observed[i] - Wavelength_tmp[counter_array]))
                  counter_array = counter_array + difference
                  if counter_array == 0:
                      Padding_Index = i
                  if counter_array >= Wavelength_tmp[-1] - Wavelength_tmp[0]:
                      break
                  Transmission_lr_tmp[counter_lr] = Transmission_tmp[counter_array]
                  counter_lr = counter_lr + 1
          Transmission_lr_tmp = Transmission_lr_tmp[0:counter_lr]
          Transmission_Padded[Padding_Index:Padding_Index + len(Transmission_lr_tmp)] = Transmission_lr_tmp
          Transmission_lr[h] = Transmission_Padded

      del Wavelength_tmp, Wavelength_lr_tmp, Transmission_tmp, Transmission_lr_tmp

        #Quickly Apply an extinction algorithm
      Coefficient = np.zeros(len(Wavelength_SED))
      test_wav = Wavelength_Observed*1e-10/1e-6
      E_B_V = 0.44
      for i in range(len(Wavelength_Observed)):              # Equations are from arxiv astro-ph/0109035
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
      Frequency_Emitted = 2.988e8/(1e-10*Wavelength_Emitted)
      Frequency_Observed = 2.998e8/(1e-10*Wavelength_Observed)
      Frequency_Emitted = np.flip(Frequency_Emitted)
      Frequency_Observed = np.flip(Frequency_Observed)
      
      Population_Flux = []
      for h in range(N_filters):
          Population_Flux_tmp = np.zeros(Spectral_Density.shape)
          for i in range(n):
              Population_Flux_tmp[i,:] = (((Spectral_Density[i,:]*Transmission_lr[h])*Wavelength_Observed*(Wavelength_Observed*1e-10)/2.998e8))
          Population_Flux.append(Population_Flux_tmp)
      
      del Population_Flux_tmp
      
      Convert_Units = 3.826e33  #Flux is in units of Solar Luminosity per Hertz. Have removed the Hertzs from integration, now to do this.
       
      h = 6.63e-34
      
      Population_Flux_Integrated = []
      
      for j in range(N_filters):
          Population_Flux_Integrated_tmp = np.zeros(n)
          Population_Flux_tmp = Population_Flux[j]
          for i in range(n):    # This for loop calibrates the found spectra to the distance to the observed galaxy.
              Population_Flux_Integrated_tmp[i] = np.trapz((((1/Frequency_Emitted)*(1+redshift)*Convert_Units*Population_Flux_tmp[i,:])/((4*np.pi*d_cm**2))),Frequency_Emitted)
          Population_Flux_Integrated.append(Population_Flux_Integrated_tmp)
          
      del Population_Flux_tmp, Population_Flux_Integrated_tmp, Population_Flux
      
      Flux_Export_List = []
      
      for i in range(N_filters):
          Population_Flux_Integrated_tmp = Population_Flux_Integrated[i]
          Transmission_tmp = Transmission_lr[i]
          Flux = Population_Flux_Integrated_tmp/np.trapz(Transmission_tmp/Frequency_Emitted, Frequency_Emitted)
          Flux_Export_List.append(Flux)
          
      del Flux, Transmission_tmp, Population_Flux_Integrated_tmp, Population_Flux_Integrated
      
      Flux_Export_Array = np.zeros([n,N_filters])
      
      for i in range(N_filters):
          Flux_Export_Array[:,i] = Flux_Export_List[i]

      return Flux_Export_Array
