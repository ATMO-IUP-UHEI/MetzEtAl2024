#!/usr/bin/env python
# regridd the 0.1x0.1° GFAS data to 1°x1° resolution
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz


import numpy as np
import pandas as pd
import datetime
import xarray as xr
from netCDF4 import Dataset
import argparse
from skimage.measure import block_reduce

#tutorial on https://unidata.github.io/netcdf4-python/netCDF4/index.html#section6

# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to put gefas in one file and rename the variables")
    args = parser.parse_args()
    return args


#settings:
args = parse_arguments()
datapath = "."

inputpath = datapath + "/GFAS/Daily_1x1_regridded/"

#settings for the first file
#newday = True
filename = datapath + "/GFAS/cams_gfas_2009_2020_daily_1x1_noscal" #output file

    
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

#number of days of month
lensy = [0,31,29,31,30,31,30,31,31,30,31,30,31]
lenny = [0,31,28,31,30,31,30,31,31,30,31,30,31]

# empty the parameter lists
adate = np.array([])
atime = np.array([])

longitude[:] = np.array(range(5,3600,10))/10
latitude[:] = np.array(range(5,1805,10))/10-90
# for each timestep
for i in range(2009,2021):
    for j in range(1,13):
        #read GFAS file
        inputfile = inputpath + 'cams_gfas_daily_regrid1x1_'+str(i)+str(j).zfill(2)+'.nc'
        DS = xr.open_mfdataset(inputfile, combine='by_coords',concat_dim='None',decode_times=False)

        if i == 2012 or i == 2016 or i == 2020:
            length = lensy[j]
        else:
            length = lenny[j]
        print(str(i)+str(j))

        #flip latitude axis (axis1)
        if len(adate) == 0:
            aco2fire = np.flip(DS.co2fire.values[0:length,:,:],axis=1)
        else:
            aco2fire = np.append(aco2fire,np.flip(DS.co2fire.values[0:length,:,:],axis=1),axis=0) 
        
        for k in range(length):
            dt = int(str((DS.time.values[k]/24 - 32872).astype('int')))
            adate = np.append(adate, 
                    int(str((datetime.datetime(1990,1,1) + datetime.timedelta(days=dt)).year).zfill(2)
                   +str((datetime.datetime(1990,1,1) + datetime.timedelta(days=dt)).month).zfill(2)
                   +str((datetime.datetime(1990,1,1) + datetime.timedelta(days=dt)).day).zfill(2)))
            atime = np.append(atime, dt)
        if i == 2020 and j == 6:
            break

#write in variables
times[:] = atime
date[:] = adate
fire[:] = aco2fire


