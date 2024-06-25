#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Read ERA5 precipitation data
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import datetime
import numpy as np
import pandas as pd
import glob, os
import geopandas
import xarray as xr
from skimage.measure import block_reduce
from RegionParam import getRegion

#IMPORTANT see: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Monthlymeans


def CreateDF(DS, nam):
    #for each timestep
    DS = DS.sel(latitude=slice(0,-35.9),longitude=slice(7,56.9))
        
    for d in range(len(DS[nam])):
        dfyCO2 = pd.DataFrame(data=block_reduce(
                                    DS.variables[nam][d].values,
                                    block_size=(10,10),
                                    func=np.mean),
                                index = list(np.array(range(-5,-365,-10))/10),#-95,-505,-10))/10), 
                                columns = list(np.array(range(75,575,10))/10),#1125,1805,10))/10),
                                dtype='float' )
        dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                    value_name =nam,var_name = 'Long')
        dfyCO2["Lat"] = dfyCO2['index']
        dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
        dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
        dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
        if d == 0:
            df = dfyCO2.copy()
        else:
            df = df.append(dfyCO2)                         
    return df

def CreateDataFrameERA5P(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Data_path, nam):
       
    if not os.path.isfile(Data_path+"/DataFrames/DF1_"+nam+"_"+RegionName+str(Num)+".pkl"):
        # check if dataframe without time variable already has been created 
        if not os.path.isfile(Data_path+"/DataFrames/DFD0_"+nam+"_"+RegionName+str(Num)+".pkl"):
            print("Start reading data:")
            filepath = Data_path+"/SoilMoisture_Temp.nc"               
            DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None')#decode_times=False)
   
            #create Dataframe
            df3= CreateDF(DS, nam)
            # preselect region as rectangle
            df = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)]
        
            print("finished reading data") 
            df.to_pickle(Data_path+"/DataFrames/DFD0_"+nam+"_"+RegionName+str(Num)+".pkl")
        else:
            df = pd.read_pickle(Data_path+"/DataFrames/DFD0_"+nam+"_"+RegionName+str(Num)+".pkl")
        #create date variable
        print("create timestamp")
        date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
        df.insert(loc=1,column='Date',value=date)
        df.to_pickle(Data_path+"/DataFrames/DF1_"+nam+"_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    else:
        df = pd.read_pickle(Data_path+"/DataFrames/DF1_"+nam+"_"+RegionName+str(Num)+".pkl") 


    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
 
       
    if Num >= 700: 
        #geodataframe (gdf) with borders of the Transcom regions
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        elif Num >= 700:
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
        
        #not needed for the interpolated gdf as gdfGrid already cut by Transcom border
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
        

    gdf.to_pickle(Data_path+"DataFrames/GDF1_"+nam+"_"+RegionName+str(Num)+".pkl")

if __name__=='__main__':
    Numm = 756
    nam = 'swvl1'
    Data_path = '/ERA5/Africa/'
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    CreateDataFrameERA5P(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, Data_path, nam)
