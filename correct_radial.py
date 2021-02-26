import numpy as np
from haversine import haversine
from create_cross_points import create_cross_latlon_points,correct_distance

def correct_radial_winds(lonpoints,latpoints,radial_slice,tc_lons,tc_lats):

   '''
    This function corrects the radial winds so that they can be correctly quivered in a cross section. If they are to the "left" of the tc center at each vertical level the sign of the winds is reveresed.

    INPUTS
    ------
    lonpoints : numpy.ndarray
      1d array of longitudes of the cross section path
    latpoints : numpy.ndarray
      1d array of latitudes of the cross section path
    radial_slice : numpy.ndarray
      2d array of the radial winds cross section (nlatlon,nvert)
    tc_lons : numpy.ndarray
      1d array of tc center longtidues (nvert)
    tc_lats : numpy.ndarray
      1d array of tc center latitudes (nvert)

    OUTPUTS
    -------
    radial_slice : numpy.ndarray
      2d array the same size as the original slice passed in with corrected radial winds
 
   '''

   #Create boolean array to see which radial winds will be reversed
   change_to_neg = np.zeros_like(radial_slice,dtype=bool)

   #Loop through vertical levels and see which locations need to be reversed
   for i,(tc_lon,tc_lat) in enumerate(zip(tc_lons,tc_lats)):
      distance = haversine(lonpoints,latpoints,tc_lons[i],tc_lats[i]) 
      this_dist = correct_distance(lonpoints,latpoints,distance,tc_lons[i],tc_lats[i])
      where_neg = np.where(this_dist<0)
      change_to_neg[where_neg,i] = True

   #Reverse 
   radial_slice[np.where(change_to_neg==True)] = -1.*radial_slice[np.where(change_to_neg==True)]

   return radial_slice


