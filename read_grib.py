import xarray as xr
import sys

def read_grib(file,vertcoord):
   try:
      data_array=xr.open_dataset(file,engine='cfgrib',
                                 backend_kwargs={'filter_by_keys':
                                       {'typeOfLevel':vertcoord}})
   except:
      print("Error when reading {} with xarray...exiting.".format(file))
      sys.exit()

   return data_array
