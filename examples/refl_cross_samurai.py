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
import color_bars
cmap_z,norm_z,ticks_z = color_bars.Tyler_colors('z')
cmap_wind,norm_wind,ticks_wind = color_bars.Tyler_colors('wind_speed')
import read_recon_obs
from read_trwind import read_trwind_file
from create_cross_points import create_cross_latlon_points
from cross_section import cross_section
from create_output_dir import create_output_directory
from read_nc_variable import read_nc_variables
from correct_radial import correct_radial_winds

####################################################
# INPUTS
####################################################
twind_file="/Users/tylergreen/Cha_triple_dop_anl/pass_4_trwind/Cha_pass_4.nc"
samurai_file="/Users/tylergreen/Cha_triple_dop_anl/pass4_triple.nc"
inbound_outbound_dist=100.
angle=-45.
hor_contour_lev = 1.0 #km
date="201610062325"
vert_coord  = "altitude"
latcoord = "latitude"
loncoord = "longitude"
####################################################

#Read trwind file
lons2d,lats2d,heights,twind,rwind,tc_lons,tc_lats = read_trwind_file(twind_file,latcoord,loncoord,vert_coord)

#For lat-lon grid-------
#Get the 1d version of lats and lons for interpolating function
lons1d = lons2d[0,:]
lats1d = lats2d[:,0]
#------------------------------

#Get points for the cross section and their distance from the center of cross section
lonpoints,latpoints,distance_from_center = create_cross_latlon_points(lons2d,lats2d,tc_lons[0],tc_lats[0],inbound_outbound_dist,angle)

#Calculate cross section
lats_samuri,lons_samuri,vardict = read_nc_variables(samurai_file,latcoord,loncoord,vert_coord,"DBZ","W")

reflectivity = np.squeeze(vardict["DBZ"])
w=np.squeeze(vardict["W"])

#Lat-lon grid
myslice_z = cross_section(lons1d,lats1d,heights,reflectivity,lonpoints,latpoints)
myslice_twind = cross_section(lons1d,lats1d,heights,twind*constants.msToKnots,lonpoints,latpoints)
myslice_rwind = cross_section(lons1d,lats1d,heights,rwind*constants.msToKnots,lonpoints,latpoints)
myslice_w = cross_section(lons1d,lats1d,heights,w*constants.msToKnots,lonpoints,latpoints)

#Correct radial winds so that the sign of it will be correct for a quiver plot
rwind_slice_corrected = correct_radial_winds(lonpoints,latpoints,myslice_rwind,tc_lons,tc_lats)

distancemesh,heightmesh = np.meshgrid(distance_from_center,heights,indexing='ij')

fig,ax = plt.subplots(figsize=(15,9))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.4)
cf=ax.contourf(distancemesh,heightmesh,myslice_z,ticks_z,cmap=cmap_z,norm=norm_z,alpha=0.9,extend='both')
ax.set_ylabel("Height (km)")
ax.set_xlabel("Distance from TC Center (km)")
ax.set_xlim(left=-inbound_outbound_dist,right=inbound_outbound_dist)
ax.set_xticks(np.arange(-inbound_outbound_dist,inbound_outbound_dist+1,25))
ax.set_yticks(np.arange(0,20,2))
ax.set_ylim(bottom=0,top=15)
ax.grid(b=True,linestyle="--",alpha=0.9)
cb = plt.colorbar(cf,cax=cax,ticks = [-10,0,10,20,30,40,50,60])
cb.set_label("Reflectivity (dBz)")
ax.contour(distancemesh,heightmesh,myslice_twind,np.arange(60,130,10),colors="black")

skipx=2
skipy=1

#Scale w wind
w_slice_corrected = myslice_w*3
ax.quiver(distancemesh[::skipx,::skipy],heightmesh[::skipx,::skipy],rwind_slice_corrected[::skipx,::skipy],w_slice_corrected[::skipx,::skipy],pivot="mid")

#inset plot of horizonta tangetial wind
axins = ax.inset_axes([0.75, 0.7, 0.28, 0.3])
m = Basemap(projection = 'cyl',llcrnrlon=tc_lons[0]-1.5,llcrnrlat=tc_lats[0]-1.5,urcrnrlon=tc_lons[0]+1.5,urcrnrlat=tc_lats[0]+1.5,resolution='i',ax=axins)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'black', linewidth = 0.5)
m.drawstates(color='black',linewidth = 0.3)
hindex = np.where(heights==hor_contour_lev)[0][0]
m.contourf(lons2d,lats2d,twind[hindex,:,:]*constants.msToKnots,ticks_wind,cmap=cmap_wind,norm=norm_wind,alpha=0.9,extend='both')
m.plot(lonpoints,latpoints,'b-')
m.plot(lonpoints[0],latpoints[0],'go',markersize=4)
m.plot(lonpoints[-1],latpoints[-1],'ro',markersize=4)

axins.text(0.2, 0.025, '{} km Vt'.format(heights[hindex]/1000.), horizontalalignment='center',verticalalignment='center', transform=axins.transAxes)
greendot_dist = -1.0*haversine(lonpoints[0],latpoints[0],tc_lons[0],tc_lats[0])
reddot_dist = haversine(lonpoints[-1],latpoints[-1],tc_lons[0],tc_lats[0])

ax.plot(greendot_dist,0.2,'go',markersize=10)
ax.plot(reddot_dist,0.2,'ro',markersize=10)

r = calc_distance_from_point(lons2d,lats2d,tc_lons[0],tc_lats[0])
m.contour(lons2d,lats2d,r,np.arange(50,200,50),colors="black",linestyles="dashed")

ax.set_title("Reflectivity Cross Section")


#Save
twind_filepath = create_output_directory(twind_file)
plt.savefig(twind_filepath+"/refl_cross_{}_{}deg.png".format(date,str(angle)),dpi=500)
#plt.show()

