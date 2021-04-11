import numpy as np
import constants
import sys

def omega2w(omega,temp,pressure,**kwargs):

   '''
    This function converts omega (pressure velocity) to vertical velocity (cartesian) assuming hydrostatic balance.

    INPUTS
    ------
    omega : np.ndarray
      array of pressure velocities (arbitrary dimension) in Pa/s
    temp : np.ndarray
      array of temperatures (arbitrary dimension) in K
    pressure : np.ndarray
      array of pressure values (arbitrary dimension) in Pa

    OUTPUTS
    -------
    w : np.ndarray
      array (same size as the input arrays) of vertical velocities (m/s)

    KWARGS
    ------
    q : np.ndarray
      array of specific humidities (kg/kg), if given, will use to calculate virtual temp

    NOTES
    -----
    - All arrays for positional arguments must be of the same shape, otherwise the function will quit
    - If q keyword arg is passed, virtual temp will be used instead of temp passed

   '''

   #Process kwargs
   known_kwargs=["q"]
   kwargs = {k:kwargs.pop(k) for k,v in list(kwargs.items()) if k in known_kwargs}

   #Make sure shapes match
   shapes = [ np.shape(omega), np.shape(temp), np.shape(pressure) ]
   allsame = all(elem == shapes[0] for elem in shapes)
   if not allsame:
      print("omega2w: Not all passed arguments have same shape, exiting...")
      sys.exit()

   #Establish constants
   g = constants.g
   R = constants.R

   #See if q is kwargs, if so calculate virutal temp
   if "q" in kwargs.keys():
      q = kwargs['q']
      temp = temp*(1+0.608*q)
   
   w = -R*temp*omega/(pressure*g)

   return w
