#Module containing functions used in the finding of vortex center and calculation of tangential winds

#Imports 
from math import radians, cos, sin, asin, sqrt
import numpy as np

# Function to find index of element closest to specified value:
def find_nearest(array, value):
    array = np.asarray(array)
    X = np.abs(array - value)
    idx = np.where( X == np.nanmin(X) )
    return array[idx]

def find_nearest_ind(array, value):
    array = np.asarray(array)
    X = np.abs(array - value)
    idx = np.where( X == np.nanmin(X) )
    return idx

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def calc_distance_from_point(lons,lats,lon0,lat0):
   [I,J] = np.shape(lons)
   dist = np.zeros_like(lons)
   for i in range(I):
      for j in range(J):
         dist[i,j] = haversine(lons[i,j],lats[i,j],lon0,lat0)

   return dist

def cut_grid(lons,lats,lon0,lat0,npoints):
   #cuts a grid (lons,lats, and var) npoints in each direction from some lat and lon point
   lons2d,lats2d = np.meshgrid(lons,lats)
   r = calc_distance_from_point(lons2d,lats2d,lon0,lat0)
   [I,J] = np.where(r == np.min(r))
   I=I[0]
   J=J[0] 

   lons_bool = np.zeros_like(lons,dtype=bool)
   lons_bool[J-npoints:J+npoints] = 1
   lats_bool = np.zeros_like(lats,dtype=bool)
   lats_bool[I-npoints:I+npoints]=1

   return lons_bool,lats_bool




 
