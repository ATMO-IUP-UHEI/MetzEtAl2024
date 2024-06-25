#!/usr/bin/env python
# Create a GeoDataFrame out of GFAS files
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

import datetime
import numpy as np
import pandas as pd
from functions import getAreaOfGrid 
from RegionParam import getRegion
import glob
import geopandas
import xarray as xr

def CreateDF(DS,Datatype):
    nam = 'co2fire'
    if 'Month' in Datatype:
        nam = 'fire'

    for d in range(len(DS[nam])):
        try:
            dfyCO2 = pd.DataFrame(data=DS.variables['co2fire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
        except:
            dfyCO2 = pd.DataFrame(data=DS.variables['fire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
        dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                  value_name ='CO2fire',var_name = 'Long')
        dfyCO2["Lat"] = dfyCO2['index']
        dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
        dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
        dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
        if 'Month' not in Datatype:
            dfyCO = pd.DataFrame(data=DS.variables['cofire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyCO = pd.melt(dfyCO.reset_index(), id_vars='index',
                                  value_name ='COfire',var_name = 'Long')
            dfyCO["Lat"] = dfyCO['index']

            dfyC = pd.DataFrame(data=DS.variables['cfire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyC = pd.melt(dfyC.reset_index(), id_vars='index',
                                  value_name ='Cfire',var_name = 'Long')
            dfyC["Lat"] = dfyC['index']

            dfy = pd.merge(dfyCO2, dfyCO, on=['Lat','Long'])
            dfy = pd.merge(dfy, dfyC, on=['Lat','Long'])
        else:
            dfy = dfyCO2

        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
        
    return df

def CreateDataFrameGFAS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    if True: #not os.path.isfile("/GFAS/dataframes/DF1_"+RegionName+str(Num)+".pkl"):
        print("Start reading data:")
        if 'Month' in Datatype:
            dp = "/GFAS/cams_gfas_CO2_2009_2020_daily_1x1_monthMeans.nc"
        else:
            dp = "/GFAS/Daily_1x1_regridded/*.nc"
        for num, filepath in enumerate(glob.glob(dp)):
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
       
            #create Dataframe
            df3= CreateDF(DS,Datatype)
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
        
       
        print("finished reading data") 

        #create date variable
        print("create timestamp")
 
        date = []
        for i,row in df.iterrows():
            date.append(datetime.date(int(row.Year),int(row.Month),int(row.Day)))    
        df.insert(loc=1,column='Date',value=date)

        Area = getAreaOfGrid()
        df = pd.merge(df, Area, on=['Lat'])     
        #        ((kg*m**-2*s**-1 days/month) -> *s/month*m**2/cell*g/kg -> gCO2/(month*gridell))
        
        CO2f = df.apply(lambda x: (x.CO2fire * 60 * 60 * 24* 
                               x.Area * 1000),axis=1)
        df.insert(loc=1, column = 'CO2fireE', value = CO2f)

        if 'Month' not in Datatype:
            Cf = df.apply(lambda x: (x.Cfire * 60 * 60 * 24* 
                                x.Area * 1000),axis=1)
            COf = df.apply(lambda x: (x.COfire * 60 * 60 * 24* 
                                x.Area * 1000),axis=1)
            df.insert(loc=1, column = 'CfireE', value = Cf)
            df.insert(loc=1, column = 'COfireE', value = COf)

        if df.Long.max() > 180:
            df = df.drop(columns = 'Long')
            df.insert(loc=1, column='Long',value = df.Long - 360)

        df.to_pickle("/GFAS/dataframes/DF2_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame

    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 700: 
        #geodataframe (gdf) with borders of the Transcom regions
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        elif Num >= 700:
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")

        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle("/GFAS/dataframes/GDF2_"+RegionName+str(Num)+".pkl")


#main
if __name__=='__main__':
    Numm  = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    CreateDataFrameGFAS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Month')