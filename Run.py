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
    
    # if ((i - 1) % 50 == 5):
    #     self.Initial_Spectral_Density = SED.Aging_initSED(self.params.n1, self.params.n1+self.params.n2, self.params.Galactic_Age_1, self.params.Galactic_Age_2, 
    #                                                      self.params.Avg_Population_Mass_1, self.params.Avg_Population_Mass_1,Spectral_Density_Array_1,Spectral_Density_Array_2)
            
    #     Spectral_Density = SED.Final_Mags_Index(self.Initial_Spectral_Density, self.params.New_Pops_Counter,self.params.New_Populations_Age,self.params.n1,
    #                                             self.Population_Mass,Spectral_Density_Array_1,Spectral_Density_Array_2)
            
    #     Population_Colours_Array,Population_Flux_Array,Counts = Colours.Colour_Calculation(Spectral_Density, Wavelength,self.params.redshift, self.params.n, self.params.d_cm, directory)
            
    #     Plotting.Plotting_Function(Population_Colours_Array, Population_Flux_Array,i,self.x0,self.params.display_scale,self.params.Galaxy_Name,
    #                                self.params.n1, self.params.n2,self.Tracer_Mass,self.params.SFR, self.params.time,self.params.Time_per_Unit,
    #                                self.params.Distance_per_Unit, self.params.rout1, self.params.rout2, self.params.Perturber_Position,directory_results,Counts)
    #     dims = self.Initial_Spectral_Density.shape
    #     self.Initial_Spectral_Density = np.zeros(dims)                        
    #     del Spectral_Density
      
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
  
  d1 = datetime.datetime.now()
  filter_data = Imports.Filters()
  for p in range(Input_Data.shape[0]):
#      if Input_Counter == 1:
#          sys.exit()
      
      args = sys.argv[1:]
      run = Run()
      run.initRun(args,Input_Counter)
      params = run.params;
      model = 3                 # Note, this line here decides the gas model. 0 = Exponential, 1 = Plummer, 2 = Hernquist and 3 = MN
      
      if os.path.isdir(directory_results+params.Galaxy_Name):
          Input_Counter += 1
          del run, params
          continue
      else:
          Results_dir = directory_results+params.Galaxy_Name
          os.mkdir(Results_dir)
          folders_make = ['\\Gas_Density', '\\Gas_Density_Primary','\\Gas_Density_Secondary','\\Lupton_Images',
                          '\\SFR_Density','\\SFR_Primary','\\SFR_Secondary','\\Colour_Images','\\Diagnostics']
          for q in range(len(folders_make)):
              os.mkdir(Results_dir+folders_make[q])
          folders_make = ['\\u_filter','\\g_filter','\\r_filter','\\i_filter','\\z_filter']
          Galaxy_Folders = ['\\Overall', '\\Primary', '\\Secondary']
          for p in range(3):
              os.mkdir(Results_dir+'\\Colour_Images'+Galaxy_Folders[p])
          for p in range(3):
              for q in range(len(folders_make)):
                  os.mkdir(Results_dir+'\\Colour_Images'+Galaxy_Folders[p]+folders_make[q])
      
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
        if i == 1:
            run.Initial_Spectral_Density = SED.Aging_initSED(params.n1, params.n1+params.n2, params.Galactic_Age_1, params.Galactic_Age_2, 
                                                        params.Avg_Population_Mass_1, params.Avg_Population_Mass_1,Spectral_Density_Array_1,Spectral_Density_Array_2)
            
            Spectral_Density = SED.Final_Mags_Index(run.Initial_Spectral_Density, params.New_Pops_Counter,params.New_Populations_Age,params.n1,
                                              run.Population_Mass,Spectral_Density_Array_1,Spectral_Density_Array_2)
            
            Population_Colours_Array,Population_Flux_Array,Counts = Colours.Colour_Calculation(filter_data,Spectral_Density, Wavelength,params.redshift, params.n, params.d_cm, directory)
            
            Plotting.Plotting_Function(Population_Colours_Array, Population_Flux_Array,run.Weights,i,run.x0,params.display_scale,params.Galaxy_Name,
                                        params.n1, params.n2,run.Tracer_Mass,params.SFR, params.time,params.Time_per_Unit,
                                        params.Distance_per_Unit, params.rout1, params.rout2, params.Perturber_Position,directory_results,Counts)
            dims = run.Initial_Spectral_Density.shape
            run.Initial_Spectral_Density = np.zeros(dims)                        
            del Spectral_Density
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
      
      Population_Colours_Array,Population_Flux_Array,Counts = Colours.Colour_Calculation(Spectral_Density, Wavelength,params.redshift, params.n, params.d_cm, directory)
      
      i = int(nstep_local)
      Plotting.Plotting_Function(Population_Colours_Array, Population_Flux_Array,run.Weights,i,run.x0,params.display_scale,params.Galaxy_Name,
                                  params.n1, params.n2,run.Tracer_Mass,params.SFR, params.time,params.Time_per_Unit,
                                  params.Distance_per_Unit, params.rout1, params.rout2, params.Perturber_Position,directory_results,Counts)      
      Input_Counter += 1
      
      # Note, the following lines are to save SFR stuff.
      # Diagnostics.SFR_Record(run.Population_Mass,params.SFH,params.e_times,[params.Galactic_Age_1, params.Galactic_Age_2],params.Time_per_Unit,params.h,params.tmin,params.n1,params.mass1,params.mass2,t0,directory_results,params.Galaxy_Name)

      del Spectral_Density_Array_1, Spectral_Density_Array_2,Wavelength,run,Spectral_Density,Population_Colours_Array,Population_Flux_Array
      d2 = datetime.datetime.now()
      print(d1)
      print(d2)
      
      sys.exit()

  
if __name__ == "__main__": main()