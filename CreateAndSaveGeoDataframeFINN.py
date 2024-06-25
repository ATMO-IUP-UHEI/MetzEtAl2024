#!/usr/bin/env python
# Create a GeoDataFrame out of FINN files
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

import datetime
import numpy as np
import pandas as pd
import glob
import geopandas
import xarray as xr
import math

def CreateDF(DS,DS1):
    year = []
    month = []
    day = []
    date = []

    for d in range(len(DS['fire'])):
        try:
            dfyCO2 = pd.DataFrame(data=DS.variables['fire'][d].values,index = DS.variables['lat'].values, columns = DS.variables['lon'].values, dtype='float' )
        except:
            dfyCO2 = pd.DataFrame(data=DS.variables['fire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
        dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                  value_name ='CO2fire',var_name = 'Long')
        dfyCO2["Lat"] = dfyCO2['index']
        dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
        dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
        dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
        try:
            dfyCO = pd.DataFrame(data=DS1.variables['fire'][d].values,index = DS1.variables['lat'].values, columns = DS1.variables['lon'].values, dtype='float' )
            dfyCO = pd.melt(dfyCO.reset_index(), id_vars='index',
                                  value_name ='COfire',var_name = 'Long')
            dfyCO["Lat"] = dfyCO['index']

            dfy = pd.merge(dfyCO2, dfyCO, on=['Lat','Long'])
        
        except:
            dfy = dfyCO2

        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
        
    return df

def CreateDataFrameFINN(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype = ''):
    print("Start reading data:")
    if 'Month' in Datatype:
        dp = "/FINN/FINN1.5_CO2_2009_2020_daily_1x1_monthMeans.nc"
    else:
        dp = "/FINN/netcdf_files/emissions-finnv1.5_daily_CO2_*.nc"
    for num, filepath in enumerate(glob.glob(dp)):
        #for num, filepath in enumerate(glob.glob("/FINN/*.nc")):
        DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
        filepath = filepath.replace("CO2","CO")       
        try:
            DS1 = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
        except:
            DS1 = ''
        #create Dataframe
        df3= CreateDF(DS, DS1)
        if Long_min >= 0:
            df2 = df3[(df3.Long >= Long_min)
                & (df3.Long <= Long_max)
                & (df3.Lat >= Lat_min)
                & (df3.Lat <= Lat_max)]
        elif Long_max < 0:
            df2 = df3[(df3.Long >=(360 + Long_min))
                & (df3.Long <= (360 + Long_max))
                & (df3.Lat >= Lat_min)
                & (df3.Lat <= Lat_max)]
        else:
            df2 = df3[(((df3.Long >=(360 + Long_min)) & (df3.Long <= 360))|
                  ((df3.Long >=(0)) & (df3.Long <= (Long_max))))
                & (df3.Lat >= Lat_min)
                & (df3.Lat <= Lat_max)]

        if num == 0:
            df = df2.copy()
        else:           
            df = df.append(df2, ignore_index=True)
        
       
    #only mentain data with good quality flag
    
    print("finished reading data") 

    #create date variable
    print("create timestamp")
 
    date = []
    for i,row in df.iterrows():
        date.append(datetime.date(int(row.Year),int(row.Month),int(row.Day)))    
    df.insert(loc=1,column='Date',value=date)

    df.to_pickle("/FINN/dataframes/DF1_"+RegionName+str(Num)+".pkl")

    long2 = []
    COf = []
    CO2f = []
    if 'Month' in Datatype:
        fac = 1  # already in kg/m^2*s
    else:
        fac = 7.3079494*10**(-22)  #*44g/(6.022*10^23)*10000cm^2/m^2*1/1000 kg/g
    for i,row in df.iterrows():
        #        (kg/mol*mol/molecules*cm**2/m**-2   * kg*m**-2*s**-1*s/month*m**2/cell*g/kg)
        try:
            COf.append(row.COfire*fac*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
        except:
            pass
        CO2f.append(row.CO2fire*fac*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
        if row.Long > 180:
            long2.append(row.Long - 360)
        else:
            long2.append(row.Long)

    df = df.drop(columns = 'Long')
    df.insert(loc=1, column='Long',value = long2)
    try:
        df.insert(loc=1, column = 'COfireE', value = COf)
    except:
        pass
    df.insert(loc=1, column = 'CO2fireE', value = CO2f)


    df.to_pickle("/FINN/dataframes/DF1_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 700:
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        else:
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle("/FINN/dataframes/GDF1_"+RegionName+str(Num)+".pkl")


