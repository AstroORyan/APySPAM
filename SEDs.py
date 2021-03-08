# -*- coding: utf-8 -*-
"""
This file holds all of the functions related to initializing, updating, or displaying the SEDs.
"""

import numpy as np

class SED:
    def initSED(tstart,tend,time_step,n1,n2,Weights,mass1,mass2,Gas_Fraction1,Gas_Fraction2,MU,length):
        # Updated
        n = n1+n2
        
        Baryonic_Fraction = (1+0.3333)/7.1333
      
        temp_1 = int(np.ceil(abs((tstart - tend)/time_step) + 2))
        Population_Mass = np.zeros([n, temp_1])  #These lines are required to store the newly created populations initial age.
        Population_Ages = np.zeros([n,temp_1])
        
    
        # First, will construct the SED Distribution of the Primary Galaxy for all particles.
        Initial_Spectral_Density = np.zeros([n,length])
      
        Avg_Population_Mass1 = (1 - Gas_Fraction1)*(Baryonic_Fraction*mass1*MU)*Weights[:n1] 
        Avg_Population_Mass2 = (1 - Gas_Fraction2)*(Baryonic_Fraction*mass2*MU)*Weights[n1:]
        
        return Population_Mass, Initial_Spectral_Density, Avg_Population_Mass1, Avg_Population_Mass2, Population_Ages
     
        
    def Aging_initSED(n1, n, Age_1, Age_2, Pop1_Mass, Pop2_Mass,Spectral_Density_Array_1,Spectral_Density_Array_2):
        # Updated
        Spectral_Age_1 = Spectral_Density_Array_1[0,1:]/1e9
        Spectral_Age_2 = Spectral_Density_Array_2[0,1:]/1e9
        
        Initial_Spectral_Density_temp = np.zeros([n,(Spectral_Density_Array_1.shape[0] - 1)])
      
        current_age_index_1 = next(x[0] for x in enumerate(Spectral_Age_1) if (x[1] >= Age_1))
        current_age_index_2 = next(x[0] for x in enumerate(Spectral_Age_2) if (x[1] >= Age_2))
        
        for q in range(n1):
            Initial_Spectral_Density_temp[q,:] = Spectral_Density_Array_1[1:,current_age_index_1]*Pop1_Mass[q]
            Initial_Spectral_Density_temp[n1+q:,:] = Spectral_Density_Array_2[1:,current_age_index_2]*Pop2_Mass[q]  
      
        return Initial_Spectral_Density_temp
  
    def Final_Mags_Index(Initial_Spectral_Density,Population_Counter,Population_Ages,n1,Population_Mass,Spectral_Density_Array_1, Spectral_Density_Array_2):
        #Updated  
        ages_SSP = Spectral_Density_Array_1[0,:]/1e9
        Total_Addition = np.zeros(Initial_Spectral_Density.shape)
        cutoff = int(n1)
       
        addition_array = np.zeros(Initial_Spectral_Density.shape)
        Spectral_Density = np.zeros(Initial_Spectral_Density.shape)
      
        Prior_Age_Index = 100
        for j in range(Population_Counter):
            temp = Population_Ages[0,j]
            New_Age_Index = next(x[0] for x in enumerate(ages_SSP) if x[1] > temp)
            if New_Age_Index == Prior_Age_Index and j > 0:
                Total_Addition += Addition
            else:
                Prior_Age_Index = New_Age_Index
                temp = Population_Ages[0,j]
                Spectral_Density_New_1 = SED.Age_Interpolation(New_Age_Index, Spectral_Density_Array_1, ages_SSP, temp)
                Spectral_Density_New_2 = SED.Age_Interpolation(New_Age_Index, Spectral_Density_Array_2, ages_SSP, temp)
                Addition = SED.Star_Formation_Flux_Calculation(addition_array, Population_Mass[:,j], cutoff, Spectral_Density_New_1, Spectral_Density_New_2)
                Total_Addition += Addition
              
    #      Total_Addition = np.transpose(Total_Addition)
        Spectral_Density = Initial_Spectral_Density + Total_Addition
      
        return Spectral_Density
  
#   @jit(forceobj=True)
    def Star_Formation_Flux_Calculation(addition_array, temp_Masses, cutoff, Spectral_Density_New_1, Spectral_Density_New_2):    
        #Updated        
    
        addition_array[0:int(cutoff),:] = Spectral_Density_New_1
        addition_array[(int(cutoff +1)):,:] = Spectral_Density_New_2
        
        Addition = temp_Masses[:,np.newaxis]*(addition_array)
        
        return Addition

    #This function will aim to do a MANUAL interpolation between the potential fluxes of two different ages.    
    def Age_Interpolation(age_index, Spectral_Density_Array,Age_Arrays,temp):
        Max_Age_Point_Array = Age_Arrays[age_index]
        Min_Age_Point_Array = Age_Arrays[age_index-1]
        
        Age_Point = temp
        
        Age_Difference_Array = Max_Age_Point_Array - Min_Age_Point_Array
        Age_Difference_Actual = Age_Point - Min_Age_Point_Array
        
        Percentage_Increase = Age_Difference_Actual/Age_Difference_Array
        
        Spectral_Density_Max = Spectral_Density_Array[1:,age_index]
        Spectral_Density_Min = Spectral_Density_Array[1:,age_index-1]
        
        Spectral_Density_Difference = Spectral_Density_Max - Spectral_Density_Min
        
        Spectral_Density_Addition = Percentage_Increase*Spectral_Density_Difference
        
        Spectral_Density_New = Spectral_Density_Min + Spectral_Density_Addition
        
        return Spectral_Density_New
      
    def Metal_Interpolation(metal,directory):
        # Updated
        Possible_Metallicity_files = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
        
        if (metal > 0.5 or metal < 0.0001):
            print('Specified metal out of bounds! (Must be between 0.0001 and .05)')
            exit()
            
        if metal in Possible_Metallicity_files:
            filename_Spectral_Density_Array = directory+'Spectral_Data_Z_'+str(metal)+'.txt'
            Spectral_Density_Array = np.loadtxt(filename_Spectral_Density_Array)
            return Spectral_Density_Array
            
        Max_Metal_Index = next(x[0] for x in enumerate(Possible_Metallicity_files) if x[1] > metal)
        Min_Metal_Index = Max_Metal_Index - 1
        
        Max_Metal = Possible_Metallicity_files[Max_Metal_Index]
        Min_Metal = Possible_Metallicity_files[Min_Metal_Index]
        
        filename_Spectral_Density_Max = directory+r'Spectral_Data_Z_'+str(Max_Metal)+'.txt'
        filename_Spectral_Density_Min = directory+r'Spectral_Data_Z_'+str(Min_Metal)+'.txt'
    
        Spectral_Density_Array_Max = np.loadtxt(filename_Spectral_Density_Max,comments='#')
        Spectral_Density_Array_Min = np.loadtxt(filename_Spectral_Density_Min,comments='#')
        
        File_Metal_Difference = Max_Metal - Min_Metal
        Spec_Metal_Difference = metal - Min_Metal
        
        Percentage_Increase = Spec_Metal_Difference/File_Metal_Difference
        
        Spectral_Density_Difference = Spectral_Density_Array_Max[1:,1:] - Spectral_Density_Array_Min[1:,1:]
        Increment = Percentage_Increase*Spectral_Density_Difference
        
        Spectral_Density_Array = Spectral_Density_Array_Min
        
        Spectral_Density_Array[1:,1:] = Spectral_Density_Array_Min[1:,1:] + Increment
        
        return Spectral_Density_Array
    