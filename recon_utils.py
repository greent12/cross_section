import pandas as pd
import numpy as np
from haversine import haversine
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator

def read_file(file):
   ''' 
   This function reads recon obs in from a GSI diag file
   NOTE: There should only be recon obs in this file.
   '''

   #Read file using pandas
   names = ['vartype','ID','type','time','lat','lon','pressure','usag','val1','inc1','val2','inc2']
   data = pd.read_csv(file,usecols=[0,2,4,5,6,7,8,9,10,11,12,13],header=None,names=names,delim_whitespace=True)

   data = data[(data.usag==1)]

   #Longitude correction if they are from 0-360
   if any(data.lon.values>180.):
      lons=data.lon.values-180
   else:
      lons=data.lon.values
  
   #Variables that will be returned
   vartype  = data.vartype.values
   time = data.time.values
   lats = data.lat.values
   lons = lons
   pres = data.pressure.values
   usag = data.usag.values
   val1 = data.val1.values
   inc1 = data.inc1.values
   val2 = data.val2.values
   inc2 = data.inc2.values

   return vartype,time,lats,lons,pres,usag,val1,inc1,val2,inc2 

def save_recon_obs(radius,heights,variable):
   ''' The point of this routine is to save the tranformed recon observations' in height/radius coords'''
   f=open("recon_obs_save.txt","w")
   for r,height,var in zip(radius,heights,variable):
      f.write("{0:7.2f} {1:7.2f} {2:7.2f}\n".format(r,height,var))
   f.close()
   return
 
def find_close_recon_obs(cross_lons,cross_lats,ob_lons,ob_lats,dthresh):

   '''
   This function will find all obs that are within a certain distance from a cross section path
   INPUTS
   1 : cross_lons : 1d numpy array
       array of the longitude values on the cross section path
   2 : cross_lats : 1d numpy array
       array of the latitude values on the cross section path
   3 : ob_lons : list or 1d numpy array 
       longitude values of the recon observation points
   4 : ob_lats : list or 1d numpy array
       latitude values of the recon observation points
   5 : dthresh : float
       represents the max distance an observation point can be away from the cross section path

   OUTPUTS
   1 : obs_lon_keep : list
       the longitudes of the observations that were kept
   2 : obs_lat_keep : list
       the latitudes of the observations that were kept
   3 : indx_keep : list
       indicies corresponding to the original lists of observations that indicate which were kept

   '''

   #Make sure ob_lons and ob_lats are both numpy arrays
   ob_lons = np.array(ob_lons)
   ob_lats = np.array(ob_lats)

   #going to loop through every observation lat/lon pair and see if they are within the threshold, if they are, were going to keep them
   obs_lon_keep=[]
   obs_lat_keep=[]
   indx_keep=[]

   for i,(ob_lon,ob_lat) in enumerate(zip(ob_lons,ob_lats)):
      iob_dist_from_cross = haversine(cross_lons,cross_lats,ob_lon,ob_lat)

      if np.min(iob_dist_from_cross) <= dthresh:
         obs_lon_keep.append(ob_lon)
         obs_lat_keep.append(ob_lat)
         indx_keep.append(i)

   obs_lon_keep = np.array(obs_lon_keep)
   obs_lat_keep = np.array(obs_lat_keep)
  
   return obs_lon_keep,obs_lat_keep,indx_keep
         
def recon_convert_pres_to_height(lons2d,lats2d,plevs,heightOnPresSurf,ob_lons,ob_lats,ob_press):

   ''' 
   This function will take recon observations (ob_lon, ob_lat, and ob_press) and convert the pressure to a height so that it can be plotted in height coords
   '''

   ob_heights=[]
   #Loop through each observation
   for ob_lon,ob_lat,ob_pres in zip(ob_lons,ob_lats,ob_press):

      #Find what point in the gird this ob is closest to 
      dist = haversine(lons2d,lats2d,ob_lon,ob_lat)
      I,J = np.where(dist == np.min(dist))
      I=I[0]
      J=J[0]
    
      #Get heights of pressure surfaces in column above this point 
      heights_pcol = heightOnPresSurf[:,I,J] 
       
      f = interpolate.interp1d(plevs,heights_pcol)

      ob_heights.append(f(ob_pres))

   ob_heights = np.array(ob_heights) 

   return ob_heights

def project_recon_onto_cross_section(ob_lons,ob_lats,lon0,lat0,angle):

   recon_r = []
   for ob_lon,ob_lat in zip(ob_lons,ob_lats):
      ob_lon_rel = ob_lon-lon0
      ob_lat_rel = ob_lat-lat0

      angle_ob = np.rad2deg(np.arctan2(ob_lat_rel,ob_lon_rel))
      dist_ob = haversine(ob_lon,ob_lat,lon0,lat0)
      dangle = np.deg2rad(angle_ob-angle)

      recon_r.append(dist_ob*np.cos(dangle))

   recon_r = np.array(recon_r)
   
   if angle == -90 or angle == 90:
      recon_r = np.flip(recon_r)

   return recon_r 

def recon_trim_convert(cross_lons,cross_lats,lons2d,lats2d,plevs,heightOnPresSurf,ob_lons,ob_lats,ob_press,dthresh,lon0,lat0,angle):
   '''
   This function will take the lats and lons of a cross section path, along with the coordinates of the orginal grid, and do 2 things:
   1 : retain on the recon obs that lie within dthresh km of the cross section path
   2 : convert the vertical coordinate from pressure tho height
   3 : convert the lat/lon coordinates into radius coordinates
   '''

   #Trim
   obs_lon_trim,obs_lat_trim,indx_keep=find_close_recon_obs(cross_lons,cross_lats,ob_lons,ob_lats,dthresh)

   #Trim recon pressures
   pres_array = np.array(ob_press)
   obs_pres_trim = pres_array[indx_keep]

   #Get pressure values
   obs_heights=recon_convert_pres_to_height(lons2d,lats2d,plevs,heightOnPresSurf,obs_lon_trim,obs_lat_trim,obs_pres_trim)

   #Convert into radius coords
   obs_r = project_recon_onto_cross_section(obs_lon_trim,obs_lat_trim,lon0,lat0,angle)

   return obs_lon_trim,obs_lat_trim,obs_heights,obs_r,indx_keep

def interp_model_cross_section_to_recon_obs(verts,distance_from_center,var,recon_radii,recon_heights):
   interp_func = RegularGridInterpolator((distance_from_center,verts),var,'linear')

   interpoints = [ [ recon_radii[i],recon_heights[i]] for i in range(len(recon_radii)) ]   

   interpolated_var = interp_func(interpoints)

   return interpolated_var
