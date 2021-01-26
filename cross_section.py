import numpy as np
from scipy.interpolate import RegularGridInterpolator,griddata
import sys

def cross_section(lons,lats,verts,variable,lonpoints,latpoints):

   '''
   This function takes the coordinates of a grid, the variable on that grid, and some lat/lon points and creates a 2d slice cross section.
 
   INPUTS
   ------
   lons : numpy.ndarray
     1 or 2D  array of longitudes for the original grid
   lats : numpy.ndarray
     1 or 2D  array of latitudes for the original grid
   verts : numpy.ndarray
     1d array of values of the vertical coordinate
   variable :: numpy.ndarray
     array shape (nvert,nlat,nlon) with the variable to be sliced
   lonpoints :: numpy.ndarray
     1d array of the longitude points to take slice at (interpolate to)
   latpoints :: numpy.ndarray
     1d array of the latitude points to take slice at (interpolate to)

   OUTPUTS
   -------
   slice : numpy.ndarray
    array shape (nlonpoints/nlatpoints,nvert) (2d slice) through the inputed latpoints and lonpoints 

   NOTES
   -----
   If the passed coordinate variables for lons/lats are 1D, it will use the scipy Regular Grid Interpolator's linear interpolation. If the passed coordinate variables are 2D, it will use scipy's griddata interpolator which does linear interpolation on with unstructured grids. This would be applicable for a lambert confromal grid, while the 1D coordinate variables would be for a lat-lon grid.

   '''

   #Make sure the varible dimensions match the coordinate dimensions
   if lats.ndim == 2 and lons.ndim == 2:
      numvert = len(verts)
      numlat = np.shape(lats)[0]
      numlon = np.shape(lons)[1]
      if not np.shape(variable) == (numvert,numlat,numlon):
         sys.exit("cross_section: passed variable does not match dimensions of passed lats, lons, and vertical coordinate... exiting")

   elif lats.ndim == 1 and lons.ndim == 1:
      numvert = len(verts)
      numlat = len(lats)
      numlon = len(lons)
      if not np.shape(variable) == (numvert,numlat,numlon):
         sys.exit("cross_section: passed variable does not match dimensions of passed lats, lons, and vertical coordinate... exiting")
   
   #Get interpolation ready
   if lats.ndim == 2 and lons.ndim == 2:
      points = np.array( (lats.flatten(), lons.flatten()) ).T
      twoD=True
   else:
      #Scipy interpolating function
      interp_func = RegularGridInterpolator((verts,lats,lons),variable,'linear')
      twoD=False
 
   #Slice points 
   slice = np.zeros((len(lonpoints),len(verts)),dtype="float")

   #Loop through vertical levels and interpolate
   print("cross_section.cross_section: starting cross section calculation")
   for zi in range(len(verts)):
      if twoD:
         values = variable[zi].flatten()
         slice[:,zi] = griddata(points,variable[zi].flatten(),(latpoints,lonpoints))
      else:
         slicepoints = [ [verts[zi],latpoints[i],lonpoints[i]] for i in range(len(latpoints))]
         slice[:,zi] = interp_func(slicepoints)

   return slice

