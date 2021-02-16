# -*- coding: utf-8 -*-
"""
This function creates a gas distribution based upon the user inputted preference. The four options are:
Exponential, Plummer, Hernquist and Miyamoto-Nagai gas distributions. The gas distribution will also 
feed into the stellar material distribution, and reflect how that is distributed.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

class Gas_Dist:
    
  def Gas_Decision(model,params,x0):
      r1 = params.rout1
      r2 = params.rout2
      n1 = params.n1
      n2 = params.n2
      Gas_Mass1, Gas_Mass2 = params.Gas_Mass_Fractions[0], params.Gas_Mass_Fractions[1]
      DU = params.Distance_per_Unit
      Sec_CM_Coords = params.Perturber_Position
      
      if model == 0:
          Tracer_Mass, Tracer_Weight = Gas_Dist.Init_Gas_Fraction_EXP(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU)
      elif model == 1:
          Tracer_Mass, Tracer_Weight = Gas_Dist.Init_Gas_Fraction_Plummer(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU)
      elif model == 2:
          Tracer_Mass, Tracer_Weight = Gas_Dist.Init_Gas_Fraction_Hernquist(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU)
      elif model == 3:
          Tracer_Mass, Tracer_Weight = Gas_Dist.Init_Gas_Fraction_Miyamoto_Nagai(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU,params.theta1,
                                                                                 params.theta2,params.phi1,params.phi2)
      else:
          print('ERROR: Please input a valid decision for the gas model (0,1,2,3).')
          sys.exit()
          
      return Tracer_Mass, Tracer_Weight
      
  def Init_Gas_Fraction_EXP(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU):
    ############################ EXTRACT FROM SELF ##############################
    n = n1+n2
    Prim_Coords = x0[:n1,:]
    Sec_Coords = x0[n1:n,:]
    Sec_Frame_Conversion = Sec_CM_Coords[0,:]
    
    ############ ALGORITHM #####################################################
    Particle_Weights = np.zeros(n)
    Tracer_Mass = np.zeros(n)
    Radii = np.zeros(n)
    
    Radii[:n1] = np.sqrt(Prim_Coords[:,0]*Prim_Coords[:,0] + Prim_Coords[:,1]*Prim_Coords[:,1] + Prim_Coords[:,2]*Prim_Coords[:,2])*DU
    
    x_sec = Sec_Coords[:,0] - Sec_Frame_Conversion[0]
    y_sec = Sec_Coords[:,1] - Sec_Frame_Conversion[1]
    z_sec = Sec_Coords[:,2] - Sec_Frame_Conversion[2]
    
    Radii[n1:] = np.sqrt(x_sec*x_sec + y_sec*y_sec + z_sec*z_sec)*DU 
    
    a1 = 0.5*r1*DU/1.69
    a2 = 0.5*r2*DU/1.69
    
    Particle_Weights[:n1] = np.exp(-Radii[:n1]/a1)
    Particle_Weights[n1:] = np.exp(-Radii[n1:]/a2) 
    
    Particle_Weights[:n1] /= np.sum(Particle_Weights[:n1])
    Particle_Weights[n1:] /= np.sum(Particle_Weights[n1:])

    Tracer_Mass[:n1] = Gas_Mass1*Particle_Weights[:n1]
    Tracer_Mass[n1:] = Gas_Mass2*Particle_Weights[n1:]
            
    return Tracer_Mass, Particle_Weights
    
  def Init_Gas_Fraction_Plummer(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU):
    ##############EXTRACT FROM OBJECT##########################################
    n = n1+n2
    Prim_Coords = x0[:n1,:]
    Sec_Coords = x0[n1:n,:]
    Sec_Frame_Conversion = Sec_CM_Coords[0,:]
    
    ############# ALGORITHM #################################################
    Particle_Weights = np.zeros(n)
    Tracer_Mass = np.zeros(n)
    Radii = np.zeros(n)

    Radii[:n1] = np.sqrt(Prim_Coords[:,0]*Prim_Coords[:,0] + Prim_Coords[:,1]*Prim_Coords[:,1] + Prim_Coords[:,2]*Prim_Coords[:,2])*DU
    
    x_sec = Sec_Coords[:,0] - Sec_Frame_Conversion[0]
    y_sec = Sec_Coords[:,1] - Sec_Frame_Conversion[1]
    z_sec = Sec_Coords[:,2] - Sec_Frame_Conversion[2]
    
    Radii[n1:] = np.sqrt(x_sec*x_sec + y_sec*y_sec + z_sec*z_sec)*DU 
    
    a1 = 0.5*(r1*DU)/1.69
    a2 = 0.5*(r2*DU)/1.69
    
    Particle_Weights[:n1] = ((3)/(4*np.pi*a1**3))*(1 + (Radii[:n1]**2/a1**2))**(-5/2)
    Particle_Weights[n1:] = ((3)/(4*np.pi*a2**3))*(1 + (Radii[n1:]**2/a2**2))**(-5/2)
    
    Particle_Weights[:n1] /= np.sum(Particle_Weights[:n1])
    Particle_Weights[n1:] /= np.sum(Particle_Weights[n1:])
    
    Tracer_Mass[:n1] = Gas_Mass1*Particle_Weights[:n1]
    Tracer_Mass[n1:] = Gas_Mass2*Particle_Weights[n1:]
        
    return Tracer_Mass, Particle_Weights

  def Init_Gas_Fraction_Hernquist(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU):
    n = n1+n2
    Prim_Coords = x0[:n1,:]
    Sec_Coords = x0[n1:n,:]
    Sec_Frame_Conversion = Sec_CM_Coords[0,:]
    
    ############## ALGORITHM ##################################################
    
    Particle_Weights = np.zeros(n)
    Tracer_Mass = np.zeros(n)
    Radii = np.zeros(n)

    Radii[:n1] = np.sqrt(Prim_Coords[:,0]*Prim_Coords[:,0] + Prim_Coords[:,1]*Prim_Coords[:,1] + Prim_Coords[:,2]*Prim_Coords[:,2])*DU 
    
    x_sec = Sec_Coords[:,0] - Sec_Frame_Conversion[0]
    y_sec = Sec_Coords[:,1] - Sec_Frame_Conversion[1]
    z_sec = Sec_Coords[:,2] - Sec_Frame_Conversion[2]
    
    Radii[n1:] = np.sqrt(x_sec*x_sec + y_sec*y_sec + z_sec*z_sec)*DU 
    
    a1 = (0.5*(r1*DU)/1.69)
    a2 = (0.5*(r2*DU)/1.69)
    
    Particle_Weights[:n1] = (1/(2*np.pi*a1**3))*((Radii[:n1]/a1)**-1)*((1 + (Radii[:n1]/a1))**-3)
    Particle_Weights[n1:] = (1/(2*np.pi*a2**3))*((Radii[n1:]/a2)**-1)*((1 + (Radii[n1:]/a2))**-3)
    
    Particle_Weights[:n1] /= np.sum(Particle_Weights[:n1])
    Particle_Weights[n1:] /= np.sum(Particle_Weights[n1:])
    
    Tracer_Mass[:n1] = Gas_Mass1*Particle_Weights[:n1]
    Tracer_Mass[n1:] = Gas_Mass2*Particle_Weights[n1:]
    
    return Tracer_Mass, Particle_Weights
    
    
  def Init_Gas_Fraction_Miyamoto_Nagai(r1,r2,n1,n2,Gas_Mass1,Gas_Mass2,x0,Sec_CM_Coords,DU,theta1,theta2,phi1,phi2):
    ################## EXTRACT REQUIRED PARAMETERS FROM OBJECT ##################
    n = n1+n2
    Prim_Coords= x0[:n1,:]
    Sec_Coords = x0[n1:n,:]
    Sec_Frame_Conversion = Sec_CM_Coords[0,:]
    
    ################### ALGORITHM ############################################    
    Particle_Weights = np.zeros(n)
    Tracer_Mass = np.zeros(n)
    Radii = np.zeros(n)
    
    # Will need the global densities of the disk galaxies to use vs the weights.
    Disk_Height = 0.1   # kpc
    Primary_Global_Density = Gas_Mass1/(np.pi*Disk_Height*(r1*DU)**2)     # Gives the global gas densities in solar masses/kpc3
    Secondary_Global_Density = Gas_Mass2/(np.pi*Disk_Height*(r2*DU)**2)
    
    # First, we must transform all coordinates into the reference plane of the disk
    ctheta = np.cos(np.deg2rad(-theta1))
    cphi = np.cos(np.deg2rad(-phi1))
    stheta = np.sin(np.deg2rad(-theta1))
    sphi = np.sin(np.deg2rad(-phi1))
    
    #Undo rotation about the y-axis:
    x_temp = Prim_Coords[:n1,0]*cphi - Prim_Coords[:n1,1]*sphi
    y_temp = Prim_Coords[:n1,0]*sphi + Prim_Coords[:n1,1]*cphi
    z_temp = Prim_Coords[:n1,2]
    
    #Undo rotation about the z-axis
    x = x_temp*ctheta + z_temp*stheta
    y = y_temp
    z = -x_temp*stheta + z_temp*ctheta
    
    # Previous step usually has some weird rounding errors in it. Removing those for safety
    for i in range(len(z)):
        if np.abs(z[i]) <= 1e-10:
            z[i] = 0
        
    # Convert into physical units
    x *= DU
    y *= DU
    z *= DU
    
    # Find the radii in reference frame of the galaxy
    Radii[:n1] = np.sqrt(x*x + y*y)
    z = z
    
    # Conduct Distribution Calculation for Primary
    a1 = 0.5*r1*DU/1.69
    b1 = 0.2*a1
    Particle_Weights[:n1] = ((b1**2)*(a1*Radii[:n1]**2 + (a1 + 3*(z**2 + b1**2)**0.5)*(a1 + (z**2+b1**2)**0.5)**2))/(4*np.pi*((Radii[:n1]**2 + (a1 + (z**2 + b1**2)**0.5)**2)**(5/2))*(z**2+b1**2)**(3/2))
    Particle_Weights[:n1] /= np.sum(Particle_Weights[:n1])
    
    # Do as above for the secondary galaxy
    # First, convert into the secondary's system of coordinates:
    x_sec = Sec_Coords[:,0] - Sec_Frame_Conversion[0]
    y_sec = Sec_Coords[:,1] - Sec_Frame_Conversion[1]
    z_sec = Sec_Coords[:,2] - Sec_Frame_Conversion[2]
    
    
    # Find the transformation angles:
    ctheta = np.cos(np.deg2rad(-theta2))
    cphi = np.cos(np.deg2rad(-phi2))
    stheta = np.sin(np.deg2rad(-theta2))
    sphi = np.sin(np.deg2rad(-phi2))
    
    #Undo the rotation about the y-axis
    x_temp = x_sec*cphi - y_sec*sphi
    y_temp = x_sec*sphi + y_sec*cphi
    z_temp = z_sec
    
    #Undo the rotation about the z-axis:
    x = x_temp*ctheta + z_temp*stheta
    y = y_temp
    z = -x_temp*stheta + z_temp*ctheta
    
    # Previous step usually has some weird rounding errors in it. Removing those for safety
    for i in range(len(z)):
        if np.abs(z[i]) <= 1e-10:
            z[i] = 0
    
    # Convert to physical Units
    x *= DU
    y *= DU
    z *= DU
    
    # Find Radii:
    Radii[n1:n] = np.sqrt(x*x + y*y)
    z = z
    
    # Conduct density distribution for secondary
    a2 = 0.5*r2*DU/1.69
    b2 = 0.2*a2
    Particle_Weights[n1:] = ((b2**2)*(a2*Radii[n1:]**2 + (a2 + 3*(z**2 + b2**2)**0.5)*(a2 + (z**2+b2**2)**0.5)**2))/(4*np.pi*((Radii[n1:]**2 + (a2 + (z**2 + b2**2)**0.5)**2)**(5/2))*(z**2+b2**2)**(3/2))
    Particle_Weights[n1:] /= np.sum(Particle_Weights[n1:])
    
    Tracer_Mass[:n1] = Gas_Mass1*Particle_Weights[:n1]
    Tracer_Mass[n1:] = Gas_Mass2*Particle_Weights[n1:]
    
    
    return Tracer_Mass, Particle_Weights