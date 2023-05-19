# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 17:07:42 2021

@author: oryan
"""

from numpy import zeros,save, asarray

class IOUtil:

  @staticmethod
  def outputParticles(filename,folder,x0,fluxes,SFRs):
    '''
    Creates a .txt with the specified filename and exports it to a function to write the resultant
    particle data to.

    Parameters
    -------------
    filename:
      Randomnly- or user-generated filename (.txt) to write simulation results to.
    folder:
      User specified folder to save results to.
    x0:
      Nx6 array containing final particle and velocities.
    fluxes:
      Nx(N_filters) array containing the each particles integrated flux in the relevant filter.
    SFRs:
      Nx1 array containing the final star formation rate of each particle.
    
    Returns
    --------
    '''
    fo = open(filename,'w')
    IOUtil.outputParticlesToFile(fo,folder,x0,fluxes,SFRs)


  @staticmethod
  def outputParticlesToFile(fo,folder,x0,fluxes,SFRs):
    '''
    Writes the particle positions, velocities, flux and star formation rates to a .txt file. Also saves
    these results to a .npy file with user specified filename in the relevant folder.

    Parameters
    -----------
    fo:
      Created .txt file from outputParticles.
    folder:
      User specified folder to save results files to.
    x0:
      Nx6 array containing particle positions and velocities.
    fluxes:
      Nx(N_filters) array containing the integrated fluxes of each particle given by the user.
    SFRs:
      An Nx1 array containing the star formation rates of each particle.

    Returns
    --------
    Outputs:
      Saves a .txt and .npy file containing an Nx(6+N_filters+1) array of the results of particle
      positions, velocities, integrated fluxes in each user provided filter and star formation rate.
    '''
    size = len(x0) - 1

    fluxes = asarray(list(fluxes.values()))
        
    output = zeros([size,x0.shape[1] + len(fluxes) + 1])
    output[:,:x0.shape[1]] = x0[:size,:].copy()
    
    for i in range(len(fluxes)):
        output[:,x0.shape[1]+i] = fluxes[i].copy()
        
    output[:,-1] = SFRs.copy()
    
    for i in range(size):
      dtmp = output[i]
      for j in range(x0.shape[1]+len(fluxes)):
        fo.write(IOUtil.formatDouble(dtmp[j]))
      # new line 
      fo.write("\n")
  
    fo.close()
    
    save(folder+'results.npy',output)

  @staticmethod
  def formatDouble(num):
    '''
    Function which converts all results into a 64bit float to be saved in the .txt file.

    Parameters
    -------------
    num:
      Value from the results array.
    
    Returns
    ---------
    float64:
      Float64 version of number from results array.
    '''
    return "%16.8e"%(num)