import pygrib 
import numpy as np

def open_grib_file(file):
   grbs=pygrib.open(file)
   return grbs

def close_grib_file(grbs):
   grbs.close()
   return

def pygrib_get_vertical_levs(grbs,name,levtype):
    message=grbs.select(name=name,typeOfLevel=levtype)
    lev_list=[]
    for submsg in message:
        lev_list.append(submsg['level'])
    return np.array(lev_list)

def pygrib_get_lat_lon(grbs,name,levtype):
    message=grbs.select(name=name,typeOfLevel=levtype)[0]
    return message.latlons()

def pygrib_get_3d_var(file,name,levtype):
    print("pygrib_util: reading {} from {}".format(name,file))

    grbs = open_grib_file(file)
    lats,lons = pygrib_get_lat_lon(grbs,name,levtype)
    levs = pygrib_get_vertical_levs(grbs,name,levtype)
    nz=len(levs)
    nlat,nlon = np.shape(lats)
    array = np.zeros((nz,nlat,nlon),dtype="float") 
    mask = np.full((nz,nlat,nlon),False,dtype=bool)

    for i,lev in enumerate(levs):
        grb=grbs.select(name=name,typeOfLevel=levtype,level=lev)[0]
        array[i] = grb.values
        mask[i] = grb.values.mask

    close_grib_file(grbs)

    output_array = np.ma.MaskedArray(array,mask)

    print("pygrib_util: all done")

    return lats,lons,levs,output_array


