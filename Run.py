import numpy as np
import sys
import os
import datetime
import matplotlib.pyplot as plt
import pandas as pd

from Parameters import Parameters
from SetupUtil import SetupUtil
from SPMModel import SPMModel
from MONDModel import MONDModel
from NBIModel import NBIModel
from Integrator import Integrator
from IOUtil import IOUtil
from SFR_DT import Delayed_Tau
from SFR_KS import KS_Model
from Gas_Dist import Gas_Dist
from SEDs import SED
from Colour import Colours
from Plotting import Plotting
from Import_Procedure import Imports
 

class Run:

  def __init__(self):
    self.params = Parameters()
    #self.forceModel = SPMModel()
    #self.forceModel = MONDModel()
    self.forceModel = NBIModel()
    self.params.potential_type=1
    self.integrator = Integrator(self.forceModel)
    x0 = None
    xout = None

  def initRun(self,args,Input_Counter):
    su = SetupUtil()

    su.setHelpers(self.forceModel, self.integrator, self.params)
    su.customCollision(args,Input_Counter)
    
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

  def takeAStep(self,i,run,Spectral_Density_Array_1, Spectral_Density_Array_2,Wavelength):
#  def takeAStep(self,i,run):    
    self.integrator.rk4(self.x0,self.xout,self.params.h)
    self.copyParticles(self.xout,self.x0)
    #time_advanced = np.log(time_advanced)
    self.params.time += self.params.h;
    
    self.params.Galactic_Age_1 += self.params.time_in_step
    self.params.Galactic_Age_2 += self.params.time_in_step
    self.params.Galactic_Ages = [self.params.Galactic_Age_1, self.params.Galactic_Age_2]
    
    
    if SFR_Algorithm == 0:
        self.params.SFR,self.Population_Mass[:,i-1],M1,M2 = Delayed_Tau.SFR_DT(self.params.Gas_Mass_Fractions,self.params.mass1,self.params.mass2,self.params.rout1,self.params.rout2,
                                                               self.params.Perturber_Position ,self.params.Distance_per_Unit,self.params.Time_per_Unit,
                                                               self.params.Galactic_Ages,self.params.e_times,self.params.h,i-1,self.Weights,self.params.n1)
    elif SFR_Algorithm == 1:
        self.params.SFR,self.Population_Mass[:,i-1],M1,M2,MT = KS_Model.SFR_KS(self.params.n1,self.params.n2,self.x0,self.params.rout1,self.params.rout2,self.params.Distance_per_Unit,
                                                                               self.Tracer_Mass,self.params.Perturber_Position,self.params.Time_per_Unit,self.params.h)
        self.Tracer_Mass -= MT
        if any(self.Tracer_Mass < 0):
            print('Somethings gone a little wrong.')
            sys.exit()
        del MT
        
    
    self.params.Gas_Mass_Fractions[0] -= M1
    self.params.Gas_Mass_Fractions[1] -= M2
    
    # self.params.SFH[i-1,0] = np.sum(self.params.SFR[:self.params.n1])
    # self.params.SFH[i-1,1] = np.sum(self.params.SFR[self.params.n1:])
    
    self.params.New_Populations_Age[:,:i-1] += self.params.time_in_step
      
  def getFilename(self,i):

    st = str(i)
    while(len(st)<3):
      st = "0" + st

    return directory+"\a"+st+".txt"

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
#def main():
  global directory, directory_results
  global SFR_Algorithm
  SFR_Algorithm = 0         # Note, a value of 0 corresponds to a delayed tau SF method and 1 is a KS one.
  directory = os.getcwd()
  directory_results = directory+'\\Results\\'
  Input_Counter = 0
  
  Inputs = pd.read_csv(r'C:\Users\oryan\Documents\PySPAM_Rewritten\All_Gal_Inputs.csv')
  Input_Data = Inputs.values
  
  filter_data = Imports.Filters()
  
  for p in range(Input_Data.shape[0]):      
      args = sys.argv[1:]
      run = Run()
      run.initRun(args,Input_Counter)
      params = run.params;
      model = 3                 # Note, this line here decides the gas model. 0 = Exponential, 1 = Plummer, 2 = Hernquist and 3 = MN
      
      Spectral_Density_Array_1, Spectral_Density_Array_2, Wavelength = Imports.SSPs([params.metal_1,params.metal_2])
                  
      nstep_local = 7500;
    
      t0 = params.tstart;
      params.nstep = ((params.tend-t0)/params.h)+1;
      nstep_local = params.nstep;
      time_interval = (params.tend-t0)*2;
      
      #IOUtil.writeParameterFile(params,"tmp.p")
      run.Tracer_Mass, run.Weights = Gas_Dist.Gas_Decision(model,params,run.x0)     # Calculates gas masses and weights and stores them in the self object.
      run.Population_Mass, run.Initial_Spectral_Density, params.Avg_Population_Mass_1, params.Avg_Population_Mass_2 = SED.initSED(params.tstart,params.tend,params.h,params.n1,params.n2,run.Weights,params.mass1,
                                                                                                                                  params.mass2,params.Gas_Fraction_1,params.Gas_Fraction_2,params.Mass_per_Unit,len(Wavelength))
      
      for i in range(1,int(nstep_local+1)):
        run.takeAStep(i,run,Spectral_Density_Array_1, Spectral_Density_Array_2,Wavelength)
        if(i % 10 == 5):
          run.params.iout = run.params.iout+1
          print(run.params.iout)
    
          #IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x)
      #sys.exit()
      run.Initial_Spectral_Density = SED.Aging_initSED(params.n1, params.n1+params.n2, params.Galactic_Age_1, params.Galactic_Age_2, 
                                                       params.Avg_Population_Mass_1, params.Avg_Population_Mass_1,Spectral_Density_Array_1,Spectral_Density_Array_2)
      
      Spectral_Density = SED.Final_Mags_Index(run.Initial_Spectral_Density, params.New_Pops_Counter,params.New_Populations_Age,params.n1,
                                              run.Population_Mass,Spectral_Density_Array_1,Spectral_Density_Array_2)
    #  print('Star Formation evaluation completed.')
      
      # run.Diagnostics(Spectral_Density,Wavelength)
      
      Population_Colours_Array,Population_Flux_Array,Counts = Colours.Colour_Calculation(filter_data, Spectral_Density, Wavelength,params.redshift, params.n, params.d_cm)
      
      i = int(nstep_local)
      
      Imports.Export(run.x0,Population_Flux_Array,run.Population_Mass,params.SFR,directory_results)
      
      
      Input_Counter += 1

      del Spectral_Density_Array_1, Spectral_Density_Array_2,Wavelength,run,Spectral_Density,Population_Colours_Array,Population_Flux_Array

  
if __name__ == "__main__": main()