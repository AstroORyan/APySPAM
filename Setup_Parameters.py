# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:06:45 2020

@author: oryan
This function is the input script for the APySPAM algorithm. Here, all of the parameters are set into the Parameter object in self. This Parameter
object is defined in the Run.py script as a section of self. self.params will only contain parameters intrinsic to the interaction.

Note, there will be a brief description of each parameter next to it but for a full breakdown of what they do see the MkDocs page for APySPAM
(*Insert link here*)

Inputs: self - The main object which will contain the different branches of variables. Here, the underlying parameters are defined as self.params.variable.

Outputs: None, but self is considered a global object of the class Run.py. Therefore, to call these parameters in main(), use params.variable.

"""

import numpy as np

class Setup_Parameters:
  def Setup_Parameters(self):      
      self.params.n = 2000  # Total number of particles in the simulation.
      
      self.params.phi1 = 97.78523#301            # Both phi and theta are the initial orienation angles of the galaxies. Theta is about the y axis and phi is about the z axis.
      self.params.theta1 = 144.6012#310.8075
      self.params.rscale1 = [0.0,0.0,1.0]       # A scaling relation used in the initial distribution of particles.
      self.params.rout1 = 0.94075#2.97913       # rout1 and rout2 are the initial radii of the primary and secondary galaxy respectively. Note, this is in Galaxy Units of 15kpc.
      self.params.mass1 = 1.3042#34.41671      # mass1 and mass2 are the total initial mass of the primary and secondary galaxy respectively. These are in Galaxy Units of 10^11M_{Solar}
      self.params.epsilon1 = 0.3
      self.params.eps1 = self.params.epsilon1*self.params.epsilon1   # The softening length used in the SPMModel
      self.params.heat1 = 0.0           # The initial heat of the galaxy. Applied as a factor for a random z velocity and position. 
      self.params.opt1 = 1              # Initial Particle distribution method. 1 = 1/r, 2 = exp(-r/rscale), 3 = exp(r*r*rscale[0] - rscale[1]*r - rscale[2])
      self.params.Galactic_Age_1 = 10   # Galaxy Age in Gyrs.
      self.params.Gas_Fraction_1 = 0.17 # The Gas Fraction of the galaxy. Assumed to be a direct percentage of the total stellar mass.
      self.params.metal_1 = 0.005       # The metallicity of the initial population as a factor of Z_{Solar}. Note, default range is 0.0001 - 0.05. Cannot be outwith this range.
      
      self.params.phi2 = 60.52239#35.5
      self.params.theta2 = 216.5012#321.9876
      self.params.rscale2 = [0.0,0.0,0.0]
      self.params.rout2 = 0.54373#4.2501
      self.params.mass2 = 1.25831#29.49607
      self.params.epsilon2 = 0.3
      self.params.eps2 = self.params.epsilon2*self.params.epsilon2
      self.params.heat2 = 0.0
      self.params.opt2 = 1
      self.params.Galactic_Age_2 = 8
      self.params.Gas_Fraction_2 = 0.14
      self.params.metal_2 = 0.05
      
      self.params.n1 = int(self.params.n*(self.params.mass1/(self.params.mass1 + self.params.mass2)))    # Number of particles in the primary and secondary galaxy.
      self.params.n2 = int(self.params.n*(self.params.mass2/(self.params.mass1 + self.params.mass2)))    # Note, these have been distributed based on galaxy mass.
      
      while (self.params.n1 + self.params.n2) < self.params.n:                  # A safety loop to make sure all particles are assigned a galaxy.
          random = np.random.random()
          if random <= 0.5:
              self.params.n1 += 1
          elif random > 0.5:
              self.params.n2 += 1
      
      r_x = -0.40739#-9.93853                # The FINAL position and velocity vectors of the SECONDARY galaxy.
      r_y = -1.92518#-4.5805                 # Note, as algorithm used backwards integration, these vectors will form the basis of how the interaction
      r_z = 1.95772#15.43348                # is calculated.

      v_x = -1.14918#-1.54025
      v_y = -1.02439#-2.97899
      v_z = 0.77963#4.69234
      
      self.params.sec_vec = [r_x, r_y, r_z, v_x, v_y, v_z]
      self.params.use_sec_vec = True    # Leave as True. This means APySPAM will use the specified velocity vector rather than trying to calculate it from
                                        # Keplarian parameters.
      self.params.h = 0.01              # Time step. 1 time step = approx. 256Myrs.
      # self.params.tstart = -5         # These lines allow for a user specified timeframe. tstart must be negative and tend must be positive or 0.
      # self.params.tend = 0#
      # self.params.time = self.params.tstart - self.params.tend
      self.params.tIsSet = False  #TRUE = User has set time, and code will use that. FALSE = algorithm will calculate time for interaction on the fly.
      
      self.params.redshift = 0.0273   #Galaxy redshift.
      
      self.params.Galaxy_Name = 'Arp_256'#'Arp_240'   # Name that results will be saved under. 
      
      self.params.display_scale = 6  # In Galaxy units (x15kpc). This will change the scale of the overall interaction plot at the end.
                                        # WARNING!!! Making this too small may crash the code. Start big, and work down to better scales.