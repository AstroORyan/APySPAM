# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 17:05:36 2021

@author: oryan
"""
import numpy as np
import sys
import datetime
import matplotlib.pyplot as plt

from Parameters import Parameters
from SetupUtil import SetupUtil
from SPMModel import SPMModel
from MONDModel import MONDModel
from NBIModel import NBIModel
from Integrator import Integrator
from IOUtil import IOUtil
from SFR_Calculation import SFR_Calculations
from Gas_Dist import Gas_Dist
from Plotting_Function import Plotting_Function
from SEDs import SED
from colour import colour
from tqdm import tqdm

class Run:

  def __init__(self):
    self.params = Parameters()
    #self.forceModel = SPMModel()
    #self.forceModel = MONDModel()
    #self.params.potential_type=2
    self.params.potential_type=1
    self.forceModel = NBIModel()
    self.integrator = Integrator(self.forceModel)
    x0 = None
    xout = None

  def initRun(self):
    '''
    Function which initialises the setup and calculation algorithms of the simulation.
    Here is where we link to other scripts like SetupUtil and ForceModel and create the 
    necessary initial conditions.

    '''
    su = SetupUtil()

    su.setHelpers(self.forceModel, self.integrator, self.params)
    su.customCollision()

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,True)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    su.createCollision()

    self.copyParticles(self.integrator.x,self.x0)

    self.forceModel.setParameters(self.params)

  # end initRun


  def initRunP(self,params):
    su = SetupUtil()

    self.params = params
    su.setHelpers(self.forceModel, self.integrator, self.params)

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,True)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    #self.x0 = [[0] * 6 for i in range(n+1)]
    #self.xout = [[0] * 6 for i in range(n+1)]
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    su.createCollision()
    

    self.copyParticles(self.integrator.x,self.x0)

    self.forceModel.setParameters(self.params)
  # end initRunP



  def getMinValues(self,params):
    su = SetupUtil()

    su.setHelpers(self.forceModel, self.integrator, self.params)

    # update the forceModel based upon passed in args
    self.forceModel = su.createForceModel(self.params.potential_type,True)

    self.forceModel.setParameters(self.params)
    self.integrator.initRKVar(self.params.n)

    n = self.params.n
    #self.x0 = [[0] * 6 for i in range(n+1)]
    #self.xout = [[0] * 6 for i in range(n+1)]
    self.x0 = np.zeros((n+1,6))
    self.xout = np.zeros((n+1,6))

    mins = su.perturberPositionVec(self.params.sec_vec, self.params.mass1, self.params.mass2, 
                                   self.params.eps1, self.params.eps2,
                                   self.params.h, self.params.n, self.params.n1, self.params.time, self.integrator.x)

    return mins;
  # end getMinValues

  #Copies the particles from x1 to x2.
  def copyParticles(self,x1,x2):
    n = len(x1)
    np.copyto(x2,x1)
    # assuming 6 members
    #for i in range(n):
      #x2[i][0] = x1[i][0]
      #x2[i][1] = x1[i][1]
      #x2[i][2] = x1[i][2]
      #x2[i][3] = x1[i][3]
      #x2[i][4] = x1[i][4]
      #x2[i][5] = x1[i][5]

  # end copyParticles

  def takeAStep(self):
    '''
    This is the main integration function of the algorithm. This function conducts the 
    steps in time of the particles and calculates the new positions and velocities of
    each of them. 
    
    This function uses the Runge-Kutta integration function found in Integrator.py.

    '''
    self.integrator.rk4(self.x0,self.xout,self.params.h)
    self.copyParticles(self.xout,self.x0)
    self.params.time = self.params.time + self.params.h;



  def getFilename(self,i):
    '''
    Parameters
    ----------
    i :
      Total number of time steps taken in the simulation.
        
    Returns
    -------
    str
        The filename by which to save the results to. Will be of the format a.N_steps.
    '''

    st = str(i)
    while(len(st)<3):
      st = "0" + st

    return "a."+st



  # Run the simulation with current parameters from tstart to tend
  def calculate(self,tstart,tend):
    t0 = 0
    time_interval = 0

    nstep_local = 7500
    self.params.tstart = tstart
    self.params.tend = tend
    t0 = self.params.tstart
    self.params.nstep = ((self.params.tend-t0)/self.params.h)+2
    nstep_local = self.params.nstep
    time_interval = (self.params.tend-t0)*2

    for i in range(nstep_local):
      self.takeAStep()


def main():
  '''
  Main function of the algorithm. Initialises all required functions and runs merger
  simulation. In this function, you will need to specify the path to the folder containing
  all algorithms.
  
  Parameters
  --------
  None.
  
  -------
  Returns
  -------
  None.
  
  Outputs
  Two files:
      a.N_steps - a .txt file containing N_particles x (7 + N_filter) array. The first 6 columns
                      correspond to each particle coordinates and velocities. The following
                      N columns correspond to particle fluxes in each filter with the last
                      column being the star formation rates of each particle.
      results.npy - As above, but in an .npy format. For ease of use with further Python
                    algorithms.
  '''
    
  folder = 'C:\\Users\\oryan\\Documents\\PySPAM_Original_Python\\APySPAM\\'
  filters = colour.get_filters(folder)
      
  run = Run()
  run.initRun()
  params = run.params;
  
  Gas_Masses, Weights = Gas_Dist.MN_Dist(params.rout1,params.rout2,params.n1,params.n,params.Gas_Mass,run.x0,params.Init_Coords)

  t0 = 0
  time_interval = 0

  nstep_local = 7500;

  movie_flag = False

  t0 = params.tstart;
  params.nstep = ((params.tend-t0)/params.h)+2;
  nstep_local = params.nstep;
  time_interval = (params.tend-t0)*2;
  #IOUtil.writeParameterFile(params,"tmp.p")
  
  Spectral_Density_1 = np.loadtxt(folder+'Spectra\\Raw_Spectral_Data_Z_'+str(params.metallicity[0])+'.txt')
  Spectral_Density_2 = np.loadtxt(folder+'Spectra\\Raw_Spectral_Data_Z_'+str(params.metallicity[1])+'.txt')

  for i in tqdm(range(1,int(nstep_local+1))):
    run.takeAStep()
      #run.params.iout = run.params.iout+1
      #print(run.params.iout)
      #IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x)
    if i % 10 == 0 and movie_flag == True:      
      SFRs, SF_Mass = SFR_Calculations.SFR(Gas_Masses,params.mass1,params.mass2,params.rout1,params.rout2,params.Seperation,params.h,time_interval/2,
                                           Weights,params.n1,params.n,params.Ages)
      
      Spectral_Density = SED.getSED(Spectral_Density_1,Spectral_Density_2,params.Ages,params.n1,params.n2,time_interval/2,Weights,[params.mass1,params.mass2],SF_Mass,params.h)
    
      Population_Flux = colour.get_colour(Spectral_Density[0],Spectral_Density[1],filters,params.redshift)
      
      Plotting_Function.plotting(run.x0,Population_Flux,SFRs,len(filters))
      
  print('Sim complete. Computing fluxes. Standby...')
  
  SFRs, SF_Mass = SFR_Calculations.SFR(Gas_Masses,params.mass1,params.mass2,params.rout1,params.rout2,params.Seperation,params.h,time_interval/2,
                              Weights,params.n1,params.n,params.Ages)
  
  Spectral_Density = SED.getSED(Spectral_Density_1,Spectral_Density_2,params.Ages,params.n1,params.n2,time_interval/2,Weights,[params.mass1,params.mass2],SF_Mass,params.h)

  Population_Flux = colour.get_colour(Spectral_Density[0],Spectral_Density[1],filters,params.redshift)
  
  Plotting_Function.plotting(run.x0,Population_Flux,SFRs,len(filters))

  IOUtil.outputParticles(run.getFilename(run.params.iout),folder,run.integrator.x,Population_Flux,SFRs)


if __name__ == "__main__": 
    main()