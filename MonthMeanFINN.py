#!/usr/bin/env python
# Create FINN files with monthly fire emissions
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz


import numpy as np
import pandas as pd
import sys
import datetime
import geopandas
import xarray as xr
import h5py 
from netCDF4 import Dataset
import read_remotec_out
import time
import calendar
import argparse
from skimage.measure import block_reduce

#tutorial on https://unidata.github.io/netcdf4-python/netCDF4/index.html#section6

datapath = "."

# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to put FINN in one file and calculate the monthly sum")
    args = parser.parse_args()
    return args


#settings:
args = parse_arguments()

inputpath = datapath + "/FINN/FINN1.5_CO2_2009_2020_daily_1x1.nc"

#settings for the first file
#newday = True
filename = datapath + "/FINN/FINN1.5_CO2_2009_2020_daily_1x1_monthMeans" #output file
   
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
DS = xr.open_mfdataset(inputfile, combine='by_coords',concat_dim='None',decode_times=False)


for year in range(2009,2021):
    print(year)
    if year == 2020:
        maxm = 7
    else:
        maxm = 13
    for month in range(1,maxm):
        if year == 2009 and month ==1:
            aco2fire = DS.where((DS.date < int(str(year)+str(month+1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).fire.sum(axis=0)
            aco2fire = aco2fire.expand_dims({'time':1})
            timen = DS.where((DS.date < int(str(year)+str(month+1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).time.values[14]
        else:
            if month <12:
                aco2fire2 = DS.where((DS.date < int(str(year)+str(month+1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).fire.sum(axis=0)
                aco2fire2 = aco2fire2.expand_dims({'time':1})
                aco2fire = np.append(aco2fire,aco2fire2,axis=0)
                timen = DS.where((DS.date < int(str(year)+str(month+1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).time.values[14]
            else:
                aco2fire2 = DS.where((DS.date < int(str(year+1)+str(1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).fire.sum(axis=0)
                aco2fire2 = aco2fire2.expand_dims({'time':1})
                aco2fire = np.append(aco2fire,aco2fire2,axis=0)
                timen = DS.where((DS.date < int(str(year+1)+str(1).zfill(2)+'01'))&(DS.date >= int(str(year)+str(month).zfill(2)+'01')),drop = True).time.values[14]

        atime = np.append(atime, timen)
        adate = np.append(adate, int(str(year)+str(month).zfill(2)+'15'))


#write in variables
times[:] = atime
date[:] = adate
fire[:] = aco2fire


