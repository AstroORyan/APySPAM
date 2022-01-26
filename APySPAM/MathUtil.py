# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 17:06:40 2021

@author: oryan
"""
import math
import numpy as np

class MathUtil:
  '''
  Builds a rotation matrix for rotating about the specified axis by the specified angle in radians.
  Utilises the angles defined in SetupUtil to rotate both galaxies with respect to the plane of the sky.
  Zero degrees is in the plane of the sky.
  '''
  def buildRot(angle, axis):
    '''
    This function builds the rotation matrix which is used to initially rotate the primary and secondary
    galaxies with respect to the sky.

    Parameters
    ------------
    angle:
      A Tuple of [Phi, Theta] used to build the rotation matrix.
    axis:
      An array of [1, 2, 3]. Used to define the x,y,z axis to rotate the galaxy about and create the 
      rotation matrix.

    Returns
    ----------
    mat:
      A 3x3 matrix used to rotate the galaxies, and update the positional information of each particle.
    '''
    X_AXIS = 1
    Y_AXIS = 2
    Z_AXIS = 3
       
    mat = np.zeros([3,3]) #[[0] * 3 for i in range(3)]
    cosa = math.cos(angle)
    sina = math.sin(angle)
        
    if axis == X_AXIS:
      mat[0][0]=1.0
      mat[1][1]=cosa
      mat[1][2]=sina
      mat[2][1]=-sina
      mat[2][2]=cosa
    elif axis == Y_AXIS:
      mat[0][0]=cosa
      mat[0][2]=-sina
      mat[1][1]=1.0
      mat[2][0]=sina
      mat[2][2]=cosa
    elif axis == Z_AXIS:
      mat[0][0]=cosa
      mat[0][1]=sina
      mat[1][0]=-sina
      mat[1][1]=cosa
      mat[2][2]=1.0
        
    return mat

    
  # Calculates cross product of two vectors.
  def cross(v1,v2):
    """
    Conducts the cross product between two vectors.

    Parameters
    --------------
    v1:
      Primary velocity vector.
    v2:
      Secondary velocity vector.

    Returns
    ----------
    v:
      Cross product of the two input vectors.
    """
    v=[0,0,0]
        
    v[0] = v1[1]*v2[2]-v1[2]*v2[1]
    v[1] = v1[2]*v2[0]-v1[0]*v2[2]
    v[2] = v1[0]*v2[1]-v1[1]*v2[0]
        
    return v
   
 
  # Calculates the dot product of two vectors.
  def dot(v1,v2):
    '''
    Computes the dot product of two vecotrs. Should be an artifact of the old version of JSPAM, and not be used anymore. Have updated to use numpy function instead.

    Parameters
    -----------
    v1:
      Primary vector.
    v2:
      Secondary vector.
    
    Returns
    ---------
    cp: The dot product of two vectors.
    '''
    size = len(v1)
    cp = 0.0
        
    for i in range(size):
      cp += v1[i]*v2[i]
        
    return cp

  # Calculates the magnitude of the vector.
  def mag(v):
    '''
    Calculates the magnitude between two vectors.
    '''
    return math.sqrt(np.dot(v,v))
   
 
  # Add two vectors.
  def add(v1, v2):
    '''
    Adds two vectors together. Again, an artifact of old JSPAM. Have removed as Python can now handle this without the need for a for loop.

    Parameters
    -----------
    v1:
      Primary vector.
    v2:
      Secondary vector.
    
    Returns
    --------
    v:
      The sum of the two vectors.
    '''
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = v1[i]+v2[i]
        
    return v


  # Sub two vectors.
  def sub(v1, v2):
    '''
    Subtracts two vectors. Should be an artifact of old JSPAM, and should be removed from code.

    Parameters
    -----------
    v1:
      Primary vector.
    v2:
      Secondary vector.
    
    Returns
    --------
    v:
      The subtracted vector.
    '''
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = v1[i]-v2[i]
        
    return v

    
  # Scale the vector by the specified multiplier.
  def scale(sc, v1):
    '''
    Scales a vector by a given scalar. Again, should be an artifact of old JSPAM and not used anymore.

    Parameters
    -----------
    sc:
      Scale to alter vector by.
    v1:
      Vector to scale/
    
    Returns
    --------
    v:
      Newly scaled vector.
    '''
    
    size = len(v1)
    v = [0]*size
   
    for i in range(size):
      v[i] = sc*v1[i]
        
    return v


  # Determine the angle of rotation between two vectors.
  def angleBetween(v1,v2):
    '''
    Uses the equation of the dot product to find the angle between two given vectors' direction. Again, should be deprecated and not used in code.
    Replaced by math or numpy modules.
    
    Parameter
    ----------
    v1:
      Primary vector.
    v2:
      Secondary vector.
    
    Returns
    ---------
    ang:
      Angle between vectors v1 and v2.
    '''
    m1 = MathUtil.mag(v1)
    m2 = MathUtil.mag(v2)
    ang = np.dot(v1,v2)/(m1*m2)
    ang = math.acos(ang)

    return ang
   
 
  # Multiplies a 3x3 matrix by a 3 vector.
  # It is mildly optimized by unrolling loops.
  def mult3(m,v):
    '''
    Multiplies a matrix and a vector. Done by unrolling loops. I believe numpy or math functions are faster than using a pre-written function like this.
    Should be depracated and not used in code.

    Parameters
    -----------
    m:
      Matrix to multiply vector by.
    v:
      Vector to multiply by.

    Returns
    -------
    b:
      New 3x3 matrix found by multiplying matrix m by vector v.
    '''
    b = [0,0,0]

    if len(v) < 3 :
      return v

    b[0] = m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2]
    b[1] = m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2]
    b[2] = m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]
        
    return b
   
 
  # Multiplies an N vector by an MxN matrix.
  def mult(m,v):
    '''
    Multiplies a N vector by an MxN matrix. Should be deprecated, and not used in main code anymore as numpy or math module functions will be more 
    optimised.

    Parameters
    ----------
    m:
      MxN matrix to be multiplied.
    v:
      N vector to multiply by matrix.
    
    Returns
    --------
    b:
      Solution to above multiplication and an MxN matrix.
    '''
    col = len( m[0])
    row = len(m)
    
    b = [0]*row
    
    for i in range(row):
      for j in range(col):
        b[i]+=m[i][j]*v[j]
        
    return b