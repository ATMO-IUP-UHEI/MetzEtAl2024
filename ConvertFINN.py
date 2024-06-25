#!/usr/bin/env python
# Convert FINN data units
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz


import numpy as np
import xarray as xr
from netCDF4 import Dataset
import argparse

#tutorial on https://unidata.github.io/netcdf4-python/

datapath = "."

# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to put FINN in one file and rename the variables")
    args = parser.parse_args()
    return args


#settings:
args = parse_arguments()
#date = args.date
#date = "201901"

inputpath = datapath + "/FINN/netcdf_files/"

#settings for the first file
#newday = True
filename = datapath + "/FINN/FINN1.5_CO2_2009_2020_daily_1x1" #output file

    
rootgrp = Dataset(filename+".nc", "w", format="NETCDF4_CLASSIC")

# set dimensions
time = rootgrp.createDimension("time", None) #unlimited size
latitude = rootgrp.createDimension("latitude", 180)
longitude = rootgrp.createDimension("longitude", 360)

# set varaibles
times = rootgrp.createVariable(varname = "time", datatype = "i4", dimensions=("time"))
times.units = "days since 1990-01-01 00:00:00.0"
times.calendas = "gregorian"

latitude = rootgrp.createVariable(varname = "latitude", datatype = "f4", dimensions=("latitude"))
latitude.units = "degrees_north"
latitude.comments = "Center latitude of the 1°x1° grid cell"
latitude.long_name = "latitude"

longitude = rootgrp.createVariable(varname = "longitude", datatype = "f4", dimensions=("longitude"))
longitude.units = "degrees_east"
longitude.comments = "Center longitude of the 1°x1° grid cell"
longitude.long_name = "longitude"

fire = rootgrp.createVariable(varname = "fire", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
fire.units = "kg m**-2 s**-1"
#fire.scale_factor = 1.0276405974730713E-9
#fire.add_offset = 3.367167181680266E-5
fire.missing_value = -32767
fire.long_name = "Wildfire flux of Carbon Dioxide"

date = rootgrp.createVariable(varname = "date", datatype = "i4", dimensions=("time"))
date.units = "YYYYMMDD"
date.long_name = "Date"


# ---------------------------------------------------------------------------------------

# empty the parameter lists
adate = np.array([])
atime = np.array([])

longitude[:] = np.array(range(5,3600,10))/10
latitude[:] = np.array(range(5,1805,10))/10-90
# for each timestep
for i in range(2009,2021):
    #read GFAS file
    if i == 2020:
        inputfile = inputpath + 'emissions-finnv1.5_daily_CO2_bb_surface_'+str(i)+'0101-'+str(i)+'0630_1x1.nc'
    else:
        inputfile = inputpath + 'emissions-finnv1.5_daily_CO2_bb_surface_'+str(i)+'0101-'+str(i)+'1231_1x1.nc'
    DS = xr.open_mfdataset(inputfile, combine='by_coords',concat_dim='None',decode_times=False)
    print(str(i))

    if len(adate) == 0:
        aco2fire = DS.fire.values*7.3079494*10**(-22)
    else:
        aco2fire = np.append(aco2fire,DS.fire.values*7.3079494*10**(-22),axis=0) # converted from molecules/cm^2*s to kg/(m^2*s) -> *44g/(6.022*10^23)*10000cm^2/m^2*1/1000 kg/g
    #atime = np.append(atime, int(str((DS.time.values/24 - 32872).astype('int'))))
    #print((DS.time.values/24 - 32872).astype('int'))
    atime = np.append(atime, DS.time.values)
    adate = np.append(adate, DS.date.values)


#write in variables
times[:] = atime
date[:] = adate
fire[:] = aco2fire


