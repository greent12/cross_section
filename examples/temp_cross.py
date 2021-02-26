#Temperature Cross section

#INPUTS to the script are as follows:
#trwind_file : path to the tr wind file holding center lons/lats and tangential and radial winds
#grib_file : the original grib file corresponding to the trwind file
# inbound_outbound_dist : the distance in km to make cross section from the center of the TC. So if it 100km the cross section will be 200km wide
# angle : the angle relative to a parallel to cut throught the storm. the code automatically finds the center of the vortex and the bottom most level and makes sure to cut through it
# hor_contour_lev : the vertical level to make an inset plot to show where the cross section is through
# date : string of the date uses for saving the figure
# vert_coord : name of vertical coordinate in original grib file
# lat_coord : name of latitue coordinate
# lon_coord : name of longitude coordinate
# plot_recon : whether or not to plot recon obs on the plot inset/cross section
# recon_file : path to the file with recon obs, if plot_recon is False this won't matter
# recon var : recon variable to plot in the cross section
# recon_ob_diff : if True, it will plot the difference, (model interpolated value at the obs points) - (the obs)

####################################################
# IMPORTS
####################################################
import sys
sys.path.insert(0, '../colors')
sys.path.insert(0, '..')
sys.path.insert(0, '../constants')

import numpy as np
import matplotlib.pyplot as plt
import constants
from haversine import haversine
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import color_bars
from recon_utils import *
from read_trwind import read_trwind_file
from create_cross_points import create_cross_latlon_points
from cross_section import cross_section
from create_output_dir import create_output_directory
from read_grib import read_grib
from pygrib_util import *

#colorbars
cmap,norm,ticks = color_bars.Tyler_colors('wind_speed')
cmap_t,norm_t,ticks_t = color_bars.Tyler_colors('temp')

####################################################
# INPUTS
####################################################
trwind_file="../test/matthew_201610061900_anl_trwind.nc"
grib_file="../test/wrfout_201610061900_anal.grib"
inbound_outbound_dist=100.
angle=45
hor_contour_lev = 1.0 #km
date="201610061900"
vert_coord  = "heightAboveSea"
latcoord = "latitude"
loncoord = "longitude"
plot_recon=True
recon_file="../test/recon_obs.txt"
recon_var="t"
recon_ob_diff=True

####################################################
# READ IN DATA FROM FILES
####################################################

#Read trwind file
lons2d,lats2d,heights,twind,rwind,tc_lons,tc_lats = read_trwind_file(trwind_file,latcoord,loncoord,vert_coord)

#Read in grib file
data_grib = read_grib(grib_file,vert_coord)

####################################################
#CREATE CROSS SECTION
####################################################

#Lat-lon grids-=========
#########################dd###########################
#Get the 1d version of lats and lons for interpolating function
lons1d = lons2d[0,:]
lats1d = lats2d[:,0]
#=---------------------

#Get points for the cross section and their distance from the center of cross section
lonpoints,latpoints,distance_from_center = create_cross_latlon_points(lons2d,lats2d,tc_lons[0],tc_lats[0],inbound_outbound_dist,angle)

#Calculate cross section

#lat-lon grid
temp = data_grib['t'].values
myslice = cross_section(lons1d,lats1d,heights,temp,lonpoints,latpoints)

#LCC grid
#temp = data_grib['t'].values
#myslice = cross_section(lons2d,lats2d,heights,temp,lonpoints,latpoints)

####################################################
#RECON OBSERVATIONS
####################################################
if plot_recon:
   print("Doing some calculations for recon obs.. one second")

   #Get heights on pressure surfaces, cant use xarray to get this variable
   lats_grib,lons_grib,plevs,heightsOnP = pygrib_get_3d_var(grib_file,"Geopotential Height","isobaricInhPa")

   #Read data out of recon file
   vartype,rtimes,rlats,rlons,rp,rusag,val1,inc1,val2,inc2=read_file(recon_file)

   #Only keep observations within 10.0 km of cross section and convert them into polar coords
   rlons,rlats,rheights,rradius,rindx=recon_trim_convert(lonpoints,latpoints,lons2d,lats2d,plevs,heightsOnP,rlons,rlats,rp,10.0,tc_lons[0],tc_lats[0],angle)

   #lons,lats were trimmed when they came out of above function, trim other variables
   vartype=vartype[rindx]
   rtimes=rtimes[rindx]
   rusag=rusag[rindx]
   val1=val1[rindx]
   inc1=inc1[rindx]
   val2=val2[rindx]
   inc2=inc2[rindx]

   #Only keep recon obs that match the wanted recon variable
   where_vartype = vartype == recon_var
   
   #Trim again
   vartype=vartype[where_vartype]
   rtimes=rtimes[where_vartype]
   rusag=rusag[where_vartype]
   val1=val1[where_vartype]
   inc1=inc1[where_vartype]
   val2=val2[where_vartype]
   inc2=inc2[where_vartype]
   rradius=rradius[where_vartype]
   rheights=rheights[where_vartype]

   #If difference between model and recon obs is wanted
   if recon_ob_diff:
      where_inrange = np.abs(rradius)<inbound_outbound_dist
      vartype=vartype[where_inrange]
      rtimes=rtimes[where_inrange]
      rusag=rusag[where_inrange]
      val1=val1[where_inrange]
      inc1=inc1[where_inrange]
      val2=val2[where_inrange]
      inc2=inc2[where_inrange]
      rradius=rradius[where_inrange]
      rheights=rheights[where_inrange]

      #Interpolate model to ob points
      slice_interp = interp_model_cross_section_to_recon_obs(heights,distance_from_center,myslice,rradius,rheights)

      #model interpolated value - ob value
      val_diff = slice_interp - val1

####################################################
#PLOT
####################################################
distancemesh,heightmesh = np.meshgrid(distance_from_center,heights,indexing='ij')

fig,ax = plt.subplots(figsize=(12,8))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.5)
cf=ax.contourf(distancemesh,heightmesh,myslice,ticks_t,cmap=cmap_t,norm=norm_t,alpha=0.9)
ax.set_ylabel("Height (km)")
ax.set_xlabel("Distance from TC Center (km)")
ax.set_xlim(left=-inbound_outbound_dist-5,right=inbound_outbound_dist+5)
ax.set_xticks(np.arange(-inbound_outbound_dist,inbound_outbound_dist+1,25))
cb = plt.colorbar(cf,cax=cax)
cb.set_label("Temperature (K)")

if plot_recon:
   if recon_ob_diff:
      scat=ax.scatter(rradius,rheights,c=val_diff,cmap="seismic",vmin=-5,vmax=5)
      cax2 = divider.append_axes("right",size="2%",pad=1)
      cb2 = plt.colorbar(scat,cax=cax2)
      cb2.set_label("Temperature Difference (K)")
   else:
      scat=ax.scatter(rradius,rheights,c = val1,cmap=cmap_t,norm=norm_t,alpha=0.9)

      

#inset plot of horizonta tangetial wind
axins = ax.inset_axes([0.75, 0.7, 0.28, 0.3])
m = Basemap(projection = 'cyl',llcrnrlon=tc_lons[0]-1.5,llcrnrlat=tc_lats[0]-1.5,urcrnrlon=tc_lons[0]+1.5,urcrnrlat=tc_lats[0]+1.5,resolution='i',ax=axins)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'black', linewidth = 0.5)
m.drawstates(color='black',linewidth = 0.3)
hindex = np.where(heights==hor_contour_lev*1000.)[0][0]
m.contourf(lons2d,lats2d,twind[hindex,:,:]*constants.msToKnots,ticks,cmap=cmap,norm=norm,alpha=0.9,extend='both')

m.plot(lonpoints,latpoints,'b-')
m.plot(lonpoints[0],latpoints[0],'go',markersize=4)
m.plot(lonpoints[-1],latpoints[-1],'ro',markersize=4)

#Read recon obs
if plot_recon:
   axins.scatter(rlons,rlats,c="black",s=3)

axins.text(0.2, 0.025, '{} km Temp'.format(heights[hindex]/1000.), horizontalalignment='center',verticalalignment='center', transform=axins.transAxes)
greendot_dist = -1.0*haversine(lonpoints[0],latpoints[0],tc_lons[0],tc_lats[0])
reddot_dist = haversine(lonpoints[-1],latpoints[-1],tc_lons[0],tc_lats[0])

ax.plot(greendot_dist,0.2,'go',markersize=10)
ax.plot(reddot_dist,0.2,'ro',markersize=10)

r = haversine(lons2d,lats2d,tc_lons[0],tc_lats[0])
m.contour(lons2d,lats2d,r,np.arange(50,200,50),colors="black",linestyles="dashed")

ax.set_title("Temperature Cross Section")

#Save
twind_filepath = create_output_directory(trwind_file)
plt.savefig(twind_filepath+"/temp_cross_{}_{}deg.png".format(date,str(angle)),dpi=500)
