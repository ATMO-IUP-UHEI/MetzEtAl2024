#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript to crate dataframe based on transcom region definition
out of the monthly CT 1x1° datasets
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import datetime
import numpy as np
import pandas as pd
from functions import getNumDayOfMonth, getAreaOfGrid
import geopandas
import xarray as xr
from RegionParam import getRegion

def CreateDF(DS):
    for d in range(len(DS['bio_flux_opt'])):
        for nam in ['bio_flux_opt','fire_flux_imp','fossil_flux_imp','ocn_flux_opt']:        
            dfyCO2 = pd.DataFrame(data=DS.variables[nam][d].values,
                                  index = DS.variables['latitude'].values, 
                                  columns = DS.variables['longitude'].values, 
                                  dtype='float' )
            dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                      value_name =nam,var_name = 'Long')
            dfyCO2["Lat"] = dfyCO2['index']
            dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
            dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
            dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
            
            if nam in 'bio_flux_opt':
                dfy = dfyCO2.copy()
            else:
                dfy = pd.merge(dfyCO2, dfy, on=['Lat','Long'])        
                dfy = dfy.drop(columns=['Month_y','Year_y','Day_y'])
                dfy = dfy.rename(columns={"Day_x":"Day","Month_x": "Month","Year_x":"Year"})
            del dfyCO2
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
                     
    return df

def CreateDataFrameCTFluxMonthly1x1(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num, Version ='2019B'):
    print("Start reading data")
    if Version == '2019B':
        filepath = "/CT2019/fluxMonthly/CT2019B.flux1x1-monthly.nc"    
        savepath = "/CT2019/dataframes"
    elif Version == '2022':
        filepath = "/CT2022/Flux/CT2022.flux1x1-monthly.nc"
        savepath = "/CT2022/dataframes"
        
    DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_components')
   
    #create Dataframe
    df3= CreateDF(DS)
    df = df3[(df3.Long >= Long_min)
     & (df3.Long <= Long_max)
     & (df3.Lat >= Lat_min)
     & (df3.Lat <= Lat_max)]
    
    print("create timestamp")
    date = df.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(x.Day)),axis=1)
    df.insert(loc=1,column='Date',value=date)
    df.insert(loc=1,column='MonthDate',value=date)

    #gdf with latitude dependent area of a 1x1° cell
    Area = getAreaOfGrid()
    df = pd.merge(df, Area, on=['Lat'])     
    #       sec/min * min/h * h/day * day/month * g/mol * m²/gridcell
    biot = df.apply(lambda x: (x.bio_flux_opt * 60 * 60 * 24* 
                               getNumDayOfMonth(int(x.Year),int(x.Month))*
                               44 * x.Area),axis=1)
    print('works')
    fft = df.apply(lambda x: (x.fossil_flux_imp * 60 * 60 * 24* 
                               getNumDayOfMonth(int(x.Year),int(x.Month))*
                               44 * x.Area),axis=1)
    ocnt = df.apply(lambda x: (x.ocn_flux_opt * 60 * 60 * 24* 
                               getNumDayOfMonth(int(x.Year),int(x.Month))*
                               44 * x.Area),axis=1)
    firet = df.apply(lambda x: (x.fire_flux_imp * 60 * 60 * 24* 
                               getNumDayOfMonth(int(x.Year),int(x.Month))*
                               44 * x.Area),axis=1)
    #       g/mol
    biots = df.apply(lambda x: (x.bio_flux_opt * 44),axis=1)
    ffts = df.apply(lambda x: (x.fossil_flux_imp * 44),axis=1)
    ocnts = df.apply(lambda x: (x.ocn_flux_opt * 44),axis=1)
    firets = df.apply(lambda x: (x.fire_flux_imp * 44),axis=1)
    
    df.insert(loc=1, column = 'Biotot', value = biot)
    df.insert(loc=1, column = 'Ocntot', value = ocnt)
    df.insert(loc=1, column = 'fftot', value = fft)
    df.insert(loc=1, column = 'firetot', value = firet)
    df.insert(loc=1, column = 'Bio', value = biots)
    df.insert(loc=1, column = 'Ocn', value = ocnts)
    df.insert(loc=1, column = 'ff', value = ffts)
    df.insert(loc=1, column = 'fire', value = firets)
    df.drop(columns= ['bio_flux_opt','fossil_flux_imp','ocn_flux_opt','fire_flux_imp'])

    df.to_pickle(savepath+"/DF2FLUXmonthly1x1_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
    
    if Num >= 700: 
        #geodataframe (gdf) with borders of the Transcom regions
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        elif Num >= 700:
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
                
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(savepath+"/GDF2FLUXmonthly1x1_"+RegionName+str(Num)+".pkl")

    del gdf,df
    
#main
if __name__=='__main__':
    Numm  = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    CreateDataFrameCTFluxMonthly1x1(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'2022')    
