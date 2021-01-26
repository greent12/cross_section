#This script takes a nc file with tangential/radial winds as well as the tc center locations for all vertical coordinates

#INPUTS to the script are as follows:
# twind_filep = the path to tangential/radial wind .nc file
# inbound_outbound_dist = the distance in km to make cross section from the center of the TC. So if it 100km the cross section will be 200km wide
# angle = the angle relative to a meridian to cut throught the storm. the code automatically finds the center of the vortex and the bottom most level and makes sure to cut through it
# hor_contour_lev = the vertical level to make an inset plot to show where the cross section is through
# date = string of the date uses for saving the figure
# plot_recon = whether or not to plot recon obs on the plot inset
# recon_file = path to the file with recon obs, if plot_recon is False this won't matter
# vert_coord = name of vertical coordinate in original grib file
import sys
sys.path.insert(0, '../colors')
sys.path.insert(0, '..')
sys.path.insert(0, '../constants')

import numpy as np
import matplotlib.pyplot as plt
import constants
from geopy.distance import great_circle
from recenter_utils import haversine,calc_distance_from_point
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
#import color_bars
#cmap,norm,ticks = color_bars.Tyler_colors('wind_speed')
import read_recon_obs
from read_trwind import read_trwind_file
from create_cross_points import create_cross_latlon_points
from cross_section import cross_section
from create_output_dir import create_output_directory
from read_grib import read_grib

####################################################
# INPUTS
####################################################
trwind_file="/Users/tylergreen/ben/test.nc"
grib_file="/Users/tylergreen/ben/wrfprs19anal_d03.04.06"
inbound_outbound_dist=100.
angle=45.0
hor_contour_lev = 1.0 #km
date="20170904"
plot_recon=False
vert_coord  = "heightAboveSea"
latcoord = "latitude"
loncoord = "longitude"
####################################################

#Read trwind file
lons2d,lats2d,heights,twind,rwind,tc_lons,tc_lats = read_trwind_file(trwind_file,latcoord,loncoord,vert_coord)

#Read in grib file
data_grib = read_grib(grib_file,vert_coord)

#Lat-lon grids-=========
#Get the 1d version of lats and lons for interpolating function
#lons1d = lons2d[0,:]
#lats1d = lats2d[:,0]
#=---------------------

#Get points for the cross section and their distance from the center of cross section
lonpoints,latpoints,distance_from_center = create_cross_latlon_points(lons2d,lats2d,tc_lons[0],tc_lats[0],inbound_outbound_dist,angle)

#Calculate cross section

#lat-lon grid
#myslice = cross_section(lons1d,lats1d,heights,twind*constants.msToKnots,lonpoints,latpoints)

#LCC grid
temp = data_grib['t'].values
myslice = cross_section(lons2d,lats2d,heights,temp,lonpoints,latpoints)
myslice_radial = cross_section(lons2d,lats2d,heights,rwind,lonpoints,latpoints)

#Plot
distancemesh,heightmesh = np.meshgrid(distance_from_center,heights,indexing='ij')

fig,ax = plt.subplots(figsize=(12,8))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.5)
cf=ax.contourf(distancemesh,heightmesh,myslice)
ax.set_ylabel("Height (km)")
ax.set_xlabel("Distance from TC Center (km)")
ax.set_xlim(left=-inbound_outbound_dist-5,right=inbound_outbound_dist+5)
ax.set_xticks(np.arange(-inbound_outbound_dist,inbound_outbound_dist+1,25))
#cb = plt.colorbar(cf,cax=cax,ticks = [7,25,40,52,64,96,125,155])
#cb.set_label("Tangential Wind Speed (kts)")

#inset plot of horizonta tangetial wind
axins = ax.inset_axes([0.75, 0.7, 0.28, 0.3])
m = Basemap(projection = 'cyl',llcrnrlon=tc_lons[0]-1.5,llcrnrlat=tc_lats[0]-1.5,urcrnrlon=tc_lons[0]+1.5,urcrnrlat=tc_lats[0]+1.5,resolution='i',ax=axins)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'black', linewidth = 0.5)
m.drawstates(color='black',linewidth = 0.3)
hindex = np.where(heights==hor_contour_lev*1000.)[0][0]
m.contourf(lons2d,lats2d,temp[hindex,:,:])
m.plot(lonpoints,latpoints,'b-')
m.plot(lonpoints[0],latpoints[0],'go',markersize=4)
m.plot(lonpoints[-1],latpoints[-1],'ro',markersize=4)

#Read recon obs
if plot_recon:
   rtimes,rlats,rlons,rp = read_recon_obs.read_file(recon_file)
   axins.scatter(rlons,rlats,c="black",s=3)

axins.text(0.2, 0.025, '{} km Temp'.format(heights[hindex]/1000.), horizontalalignment='center',verticalalignment='center', transform=axins.transAxes)
greendot_dist = -1.0*haversine(lonpoints[0],latpoints[0],tc_lons[0],tc_lats[0])
reddot_dist = haversine(lonpoints[-1],latpoints[-1],tc_lons[0],tc_lats[0])

ax.plot(greendot_dist,0.2,'go',markersize=10)
ax.plot(reddot_dist,0.2,'ro',markersize=10)

r = calc_distance_from_point(lons2d,lats2d,tc_lons[0],tc_lats[0])
m.contour(lons2d,lats2d,r,np.arange(50,200,50),colors="black",linestyles="dashed")

ax.set_title("Temperature Cross Section")


#Save
twind_filepath = create_output_directory(trwind_file)
plt.savefig(twind_filepath+"/temp_cross_{}_{}deg.png".format(date,str(angle)),dpi=500)
