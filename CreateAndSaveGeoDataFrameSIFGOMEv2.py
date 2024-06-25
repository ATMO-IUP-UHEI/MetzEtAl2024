#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create Pandas dataframe from GOME-2 SIF data
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import xarray as xr
import pandas as pd
import numpy as np
import datetime
import geopandas
from RegionParam import getRegion
import glob



def CreateDF(DS):
    d = {'SIF_daily_av': DS.Daily_Averaged_SIF.values, 
         'SIF_740': DS.SIF_740.values, 
         'Lat':DS.Latitude.values,
         'Long':DS.Longitude.values,
         'time':DS.Delta_Time.values,
         }
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameSIFGOMEv2(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    for num, filepath in enumerate(glob.glob('/SIF/GOME2_v2/*_all.nc')):
        DS0 = xr.open_dataset(filepath)
        #print('read data done')
        #Quality filter
        DS = DS0.where(DS0.Quality_Flag == 2).dropna(dim='obs')
        #print('QF done')
        
        #create Dataframe
        df3= CreateDF(DS)
        df2 = df3[(df3.Long >= Long_min)
                & (df3.Long <= Long_max)
                & (df3.Lat >= Lat_min)
                & (df3.Lat <= Lat_max)]
        if num == 0:
            df = df2.copy(deep = True)
        else:      
            df = df.append(df2, ignore_index=True)
            
    print("finished reading data") 
    print(df.keys())
    #create date variable
    print("create timestamp")
    year = df.apply(lambda x: int(str(x.time)[0:4]),axis = 1)
    month = df.apply(lambda x: int(str(x.time)[5:7]),axis = 1)
    day = df.apply(lambda x: int(str(x.time)[8:10]),axis = 1)
    df.insert(loc=1,column='Year',value=year)
    df.insert(loc=1,column='Month',value=month)
    df.insert(loc=1,column='Day',value=day)
    
    date = df.apply(lambda x: datetime.date(int(x.Year), int(x.Month),int(x.Day)),axis=1)
    df.insert(loc=1,column='Date',value=date)
            
    df.to_pickle("/SIF/GOME2_v2/DataFrames/DF2_"+RegionName+str(Num)+".pkl")

    df.insert(1, column = 'Lat_round', value = df.apply(lambda x: int(x.Lat),axis = 1))
    df.insert(1, column = 'Long_round', value = df.apply(lambda x: int(x.Long),axis = 1))
    df1 = df.groupby(['Lat_round','Long_round','Date','Year','Month','Day'])['SIF_740','SIF_daily_av'].mean().reset_index()
    df1.rename(columns={'Lat_round':'Lat','Long_round':'Long'},inplace=True)

    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df1, geometry=geopandas.points_from_xy(df1.Long, df1.Lat))
    gdf.crs = {'init' :'epsg:4326'}

    if Num >= 700:
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        else:
            Transcom = pd.read_pickle("/CTTranscomMask1x1Borders.pkl")
        igdf = gdf.intersects(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        #igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]

    gdf.to_pickle("/SIF/GOME2_v2/DataFrames/GDF2_"+RegionName+str(Num)+".pkl")

if __name__=='__main__':
    Numm  = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    CreateDataFrameSIFGOMEv2(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
