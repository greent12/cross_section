from recenter_utils import haversine
import numpy as np
import matplotlib.pyplot as plt
import sys
import constants
from geopy.distance import great_circle

def get_valid_cross_points(lons1d,lats1d,lonpoints,latpoints):

   minlon = np.min(lons1d)
   maxlon = np.max(lons1d)
   minlat = np.min(lats1d)
   maxlat = np.max(lats1d)

   validpoints = np.zeros(len(lonpoints),dtype=bool)

   for i,(lon,lat) in enumerate(zip(lonpoints,latpoints)):
      if lon<=maxlon and lon>=minlon and lat<=maxlat and lat>=minlat:
        validpoints[i] = 1

   return validpoints

def create_start_end_cross(center,inbound_outbound,angle):
   """
    This function takes the TC center, inbound/outbound distance in km, and angle relative to a meridian, and calculates the cross section start and end points

   INPUTS:
   1: center, tuple if tc center, (lat,lon)
   2: inbound_outbound, how far in km you want the cross section to go out on each side of the TC center
   3: angle, angle in degrees relative to a meridian that you want the cross section to go, can be between -90 and 90

   """

   #Input checks
   if type(center) is not tuple:
      print("'center' input needs to be a tuple, exitting...")
      sys.exit()
   else:
      if len(center) != 2:
         print("'center' tuple must be of length 2, exitting...")
         sys.exit()

   if np.abs(angle) > 90.:
      print("Angle has to be between -90 and 90, exitting...")

   if angle == -90:
      angle = 90
   
   dang=0.01 #Degrees
   angd=0.0
   #right side of cross section
   while True:
      dlat=angd*np.sin(np.deg2rad(angle))
      dlon=angd*np.cos(np.deg2rad(angle))
      lat = center[0]+dlat
      lon = center[1]+dlon

      dist = haversine(center[1],center[0],lon,lat)

      if dist > inbound_outbound:
         break 
      angd+=dang

   end=(lat,lon)

   #left side of cross section
   angle =  angle-180.
   if angle < -180.:
      angle = angle + 360
   while True:
      dlat=angd*np.sin(np.deg2rad(angle))
      dlon=angd*np.cos(np.deg2rad(angle))
      lat = center[0]+dlat
      lon = center[1]+dlon

      dist = haversine(center[1],center[0],lon,lat)

      if dist > inbound_outbound:
         break
      angd+=dang   

   start=(lat,lon)

   return start,end


def create_cross_latlon_points(lons2d,lats2d,lon_center,lat_center,inbound_outbound,angle):

   ''' 

     This function calculates the lat/lon points for a cross section as well as the radial distance of these points from the center point of the cross section. 

   INPUTS
   ------
   lons2d : numpy.ndarray
     2d mesh of longitudes
   lats2d : numpy.ndarray
     2d mesh of latitudes
   lon_center : float
     longitude of the center that the cross section will go through
   lat_center : float
     latitude of the center that the cross section will go through
   inbound_outbound : float
     distance from the center on both sides of the center in km that the cross section will go out (roughly)
   angle : float
     angle relative to a meridian that the ross section will go through, should be between -90 and 90

   OUTPUTS
   -------
   lonpoints :  list
     list of longitude points for the cross section path
   latpoints : list
     list of latitude points for the cross section path
   distance_from_center : list
     list of distances (in km) of the latpoint/lonpoint pairs from the center of the cross section
 
   '''

   #Get 1d arrays of lats and lons
   lats1d = lats2d[:,0]
   lons1d = lons2d[0,:]

   #Find start and end points for cross section. This will be through the center and 
   # a certain distance (inbound_outbound) from the center on both sides
   start,end = create_start_end_cross((lat_center,lon_center),inbound_outbound,angle)

   #Determine number of points in the cross section depending on grid spacing
   dy = np.max(lats2d[1::]-lats2d[0:-1])*2*np.pi*constants.RE/360.0
   npoints = int(np.floor(great_circle(start,end).meters/dy))

   #Set up points for cross section
   distance = np.sqrt((end[0]-start[0])**2+(end[1]-start[1])**2)
   dr = distance/(npoints-1)

   theta = np.arctan2(end[0]-start[0],end[1]-start[1])
   dlat = dr*np.sin(theta)
   dlon = dr*np.cos(theta)

   #Lat and lon points
   latpoints = np.array([ start[0]+i*dlat for i in range(npoints)])
   lonpoints = np.array([ start[1]+i*dlon for i in range(npoints)])

   #Check that all lats and lons are within domain
   valid_points = get_valid_cross_points(lons1d,lats1d,lonpoints,latpoints)
   if not len(lonpoints) == np.sum(valid_points):
      print("inbound_outbound distance wanted for cross section went off the available grid, cutting it to be shorter so it will fit on grid")
  
   #Only take lat and lon points that are within the domain, otherwise the interpolation will fail
   lonpoints = lonpoints[np.ix_(valid_points)]
   latpoints = latpoints[np.ix_(valid_points)]
    
   #Calculate the distance of each lat/lon point form the cross section center
   distance_from_center = np.array([ haversine(lonpoints[i],latpoints[i],lon_center,lat_center) for i in range(len(latpoints)) ])
   distance_from_center[0:int(np.floor(len(distance_from_center)/2))] = -distance_from_center[0:int(np.floor(len(distance_from_center)/2))]

   return lonpoints,latpoints,distance_from_center
