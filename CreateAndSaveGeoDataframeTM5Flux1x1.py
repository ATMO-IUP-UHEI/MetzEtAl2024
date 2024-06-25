#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript to crate geopandas dataframe based on transcom region definition
out of the monthly 1x1° datasets of TM5-4DVar
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import datetime
#from datetime import timedelta
import numpy as np
import pandas as pd
from functions import getNumDayOfMonth
import geopandas
import xarray as xr
from RegionParam import getRegion
from shapely.geometry import Polygon

def CreateDF(DS,Lat_min,Lat_max,Long_min,Long_max):
    DS = DS.assign_coords(longitude=np.round(DS.longitude - 179.5,1))
    #DS = DS.assign_coords(latitude=np.round(-DS.latitude + 90,1))
    DS = DS.assign_coords(latitude=np.round(DS.latitude - 89.5,1))
    DS = DS.sortby(['longitude','latitude'])
    DS = DS.sel(latitude=slice(Lat_min,Lat_max),longitude=slice(Long_min,Long_max))
    DS = DS.drop(['CO2_flux_fire_global_total',
                  'CO2_flux_oce_global_total',
                  'CO2_flux_fos_global_total',
                  'CO2_flux_nee_global_total'])
    DS2 = DS.drop(['month_tuple',
                  'grid_cell_area'])
    DSmonth = DS.drop(['CO2_flux_fire',
                  'CO2_flux_oce',
                  'CO2_flux_fos',
                  'CO2_flux_nee',
                  'grid_cell_area'])
    DSgrid = DS.drop(['CO2_flux_fire',
                  'CO2_flux_oce',
                  'CO2_flux_fos',
                  'CO2_flux_nee',
                  'month_tuple'])
    df = DS2.to_dataframe().reset_index()
    dfmonth = pd.DataFrame(data={'Year' : [x[0] for x in DSmonth.month_tuple.values],
                                 'Month' : [x[1] for x in DSmonth.month_tuple.values],
                                 'months' : DSmonth.months.values})
    dfgrid = DSgrid.to_dataframe().reset_index()
    df = pd.merge(df,dfmonth,on='months')
    df = pd.merge(df,dfgrid,on=['latitude','longitude'])
    df = df.rename(columns = {'latitude':'Lat','longitude':'Long'})
    df.insert(loc = 1, 
              column = 'CO2_flux', 
              value = df.CO2_flux_nee + df.CO2_flux_oce + df.CO2_flux_fire + df.CO2_flux_fos)
    # from gC/m^2/day to absolute fluxes (gCO2)
    df.insert (loc = 1, column='DaysInMonth', value = df.apply(lambda x: getNumDayOfMonth(int(x.Year),int(x.Month)),axis = 1))
    total_flux = df.apply(lambda x: (x.CO2_flux * x.DaysInMonth * 
                                     44/12 * x.grid_cell_area),axis=1)
    total_flux_nee = df.apply(lambda x: (x.CO2_flux_nee * x.DaysInMonth * 
                                     44/12 * x.grid_cell_area),axis=1)
    total_flux_oce = df.apply(lambda x: (x.CO2_flux_oce * x.DaysInMonth * 
                                     44/12 * x.grid_cell_area),axis=1)
    total_flux_fire = df.apply(lambda x: (x.CO2_flux_fire * x.DaysInMonth * 
                                     44/12 * x.grid_cell_area),axis=1)
    total_flux_fos = df.apply(lambda x: (x.CO2_flux_fos * x.DaysInMonth * 
                                     44/12 * x.grid_cell_area),axis=1)                                                                                                        
    df.insert(loc = 1, column = 'CO2_flux_total',value = total_flux)
    df.insert(loc = 1, column = 'CO2_flux_nee_total',value = total_flux_nee)
    df.insert(loc = 1, column = 'CO2_flux_oce_total',value = total_flux_oce)
    df.insert(loc = 1, column = 'CO2_flux_fire_total',value = total_flux_fire)
    df.insert(loc = 1, column = 'CO2_flux_fos_total',value = total_flux_fos)
                   
    return df

def CreateDataFrameTM5FluxMonthly1x1(dataset,Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    if True:#if False use dataframe already created to create GeoDataFrame
        print("Start reading data")
        filepath = "/"+dataset+".nc"
        print(filepath)
        DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_components')
    
        #create Dataframe
        df= CreateDF(DS,Lat_min,Lat_max,Long_min,Long_max)
        date = df.apply(lambda x: datetime.date(int(x.Year),
                                                        int(x.Month),
                                                        15),axis=1)
        df.insert(loc=1,column='Date',value=date)
        df.insert(loc=1,column='MonthDate',value=date)     
        print('calculate gdf')
        df.to_pickle("/dataframes/DF2FLUXmonthly1x1_"+dataset+RegionName+str(Num)+".pkl")
    else:
        df = pd.read_pickle("/dataframes/DF2FLUXmonthly1x1_"+dataset+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
    geom = gdf.apply(lambda x: Polygon(zip([x.Long-0.5,x.Long-0.5,x.Long+0.5,x.Long+0.5],[x.Lat-0.5,x.Lat+0.5,x.Lat+0.5,x.Lat-0.5])),axis = 1)
    gdf.insert(loc = 1, column = 'geomPoly',value = geom)
    gdf = gdf.set_geometry(gdf.geomPoly)

    if Num >= 700:
        if Numm >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        else:
            Transcom = pd.read_pickle("/CTTranscomMask1x1Borders.pkl")
        igdf = gdf.intersects(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        #igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
    
    gdf.to_pickle("/dataframes/GDF2FLUXmonthly1x1_"+dataset+RegionName+str(Num)+"V2.pkl")
    
    del gdf,df
    
if __name__=='__main__':
    Numm  = 784
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    for ds in ["flux_1x1_RemoTeC_2.4.0+IS","flux_1x1_ACOS+IS","flux_1x1_IS"]:
        dataset = ds#"flux_1x1_RemoTeC+IS"
        CreateDataFrameTM5FluxMonthly1x1(dataset,Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)    
    
