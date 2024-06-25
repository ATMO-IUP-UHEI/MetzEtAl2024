#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of CAMS files
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

import datetime
import numpy as np
import pandas as pd
from functions import getGdfOfGrid, getAreaOfGrid 
import glob
import geopandas
import xarray as xr
from RegionParam import getRegion
from shapely.geometry import Polygon

def CreateDF(DS2,DS,Datatype,year,month):
     
    dArea = pd.DataFrame(data=DS2.variables['area'].values,index = DS.variables['latitude'].values,
                                  columns =  DS.variables['longitude'].values, dtype='float' )
    dArea = pd.melt(dArea.reset_index(), id_vars='index',
                                 value_name ='area',var_name = 'Long')
    dArea["Lat"] = dArea['index']
    
    dArea2 = pd.DataFrame(data=DS2.variables['lsf'].values,index = DS.variables['latitude'].values,
                                  columns =  DS.variables['longitude'].values, dtype='float' )
    dArea2 = pd.melt(dArea2.reset_index(), id_vars='index',
                                 value_name ='lsf',var_name = 'Long')
    dArea2["Lat"] = dArea2['index']
    dArea = pd.merge(dArea,dArea2,on=['Lat','Long'],how='left')
    #print(dfy2)
        
    for d in [0]:#range(len(DS['flux_apos_bio'])):
        for nam in ['flux_apos_bio','flux_apri_bio','flux_apos_ocean','flux_apri_ocean','flux_foss']:        
            
            #dfyCO2 = pd.DataFrame(data=block_reduce(DS.variables[nam][d].values,block_size=(4,5),func=np.sum),index = DS.variables['latitude'].value, columns = DS.variables['longitude'].values , dtype='float' )
            dfyCO2 = pd.DataFrame(data=DS.variables[nam].values*1000*(44/12),index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            #dfyCO2 = pd.DataFrame(data=DS.variables[nam][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                      value_name =nam+'0',var_name = 'Long')
            dfyCO2["Lat"] = dfyCO2['index']
            dfyCO2 = pd.merge(dfyCO2, dArea, on=['Lat','Long'], how = 'left')
            #dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
            #dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
            dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*year)        
            dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*month)
            #dfyCO2.insert(loc=1,column=nam+'_totlsf',value=dfyCO2[nam+'0']*dfyCO2['area']*dfyCO2.lsf)
            #dfyCO2.insert(loc=1,column=nam+'_tot',value=dfyCO2[nam+'0']*dfyCO2['area'])

            if nam in 'flux_apos_bio':
                dfy = dfyCO2.copy()
            else:
                dfy = pd.merge(dfyCO2, dfy, on=['Lat','Long'])        
                dfy = dfy.drop(columns=['Month_y','Year_y'])
                dfy = dfy.rename(columns={"Month_x": "Month","Year_x":"Year"})
            del dfyCO2
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
                     
    return df

def CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    if Datatype in 'Sat':
        starty = 2014
        paths = "/CAMS/Flux/cams73_latest_co2_flux_satellite_mm_"
        DS2 = xr.open_mfdataset("/CAMS/Flux/cams73_latest_co2_flux_satellite_mm_201905.nc",combine='by_coords',concat_dim='None',decode_times=False)
    
    else:
        starty = 2009
        paths = "/CAMS/Flux/cams73_latest_co2_flux_surface_mm_"
        DS2 = xr.open_mfdataset("/CAMS/Flux/cams73_latest_co2_flux_surface_mm_201905.nc",combine='by_coords',concat_dim='None',decode_times=False)
    
    #for y in range(starty,2020): #not os.path.isfile("/GFAS/dataframes/DF21x1_"+RegionName+str(Num)+".pkl"):
    for y in range(starty,2020): #not os.path.isfile("/GFAS/dataframes/DF21x1_"+RegionName+str(Num)+".pkl"):
        print("Start reading data for : " +str(y))
        dp = paths+str(y)+"*.nc"
        for num, filepath in enumerate(glob.glob(dp)):
            #print(filepath)
            DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_components')#decode_times=False)
            mon = int(filepath[-5:-3])
            #create Dataframe
            df3= CreateDF(DS2,DS,Datatype,y,mon)
            del DS
            df2 = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)]

            if num == 0:
                df = df2.copy()
            else:           
                df = df.append(df2, ignore_index=True)
        
        print("finished reading data") 

        #create date variable
        print("create timestamp")
        date = df.apply(lambda x: datetime.date(int(x.Year),
                                                  int(x.Month),
                                                  15),axis=1)
        df.insert(loc = 1, value = date, column = 'Date') 
        df.to_pickle("/CAMS/dataframes/DF21x1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")
        df.drop(columns = ['index_x_x','index_y_x','index_y_y','index_x_y','lsf_x','lsf_y'],inplace=True)
        #create GeoDataFrame
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)

        print('Interpolate on 1°x1° grid')
        #get reference 1x1° grid geodataframe of the region
        try:
            gdfGrid = pd.read_pickle("/Grid1_1"+RegionName+".pkl")
        except:
            gdfGrid = getGdfOfGrid(Num)
        #get grid cell size
        dLong = (DS2.variables['longitude'].values[2]
                -DS2.variables['longitude'].values[1])
        dLat = (DS2.variables['latitude'].values[2]
                -DS2.variables['latitude'].values[1])
        #create polygon geomatry for geodataframe
        geom = gdf.apply(lambda x: Polygon(zip(
                [x.Long-dLong/2,x.Long-dLong/2,x.Long+dLong/2,x.Long+dLong/2],
                [x.Lat-dLat/2,x.Lat+dLat/2,x.Lat+dLat/2,x.Lat-dLat/2])),axis = 1)
        gdf.insert(loc = 1, value = geom, column = 'geomPoly')
        gdf = gdf.set_geometry(gdf.geomPoly)
        #cut the gdf grid with the 1x1° grid and keep all subgridcells of the intersection
        inter = geopandas.overlay(gdfGrid,gdf,how='intersection')
        #fig2,ax2=plt.subplots(figsize=(12,9)); inter[(inter.Year == 2012)&(inter.Month == 1)].boundary.plot(color='blue',ax = ax2); inter[(np.isnan(inter.gpp))&(inter.Year == 2012)&(inter.Month == 1)].boundary.plot(color='red',ax = ax2)
        #fig2,ax2=plt.subplots(figsize=(12,9)); gdf[(gdf.Year == 2012)&(gdf.Month == 1)].boundary.plot(color='red',ax = ax2);inter[(inter.Year == 2012)&(inter.Month == 1)].boundary.plot(color='blue',ax = ax2)
        # to m² insteadt of degree to get area of the subgridcells
        inter = inter.to_crs({'proj':'cea'}) 
        inter.insert(loc = 1, value = inter.area,column = 'subarea')
        inter = inter.to_crs(crs=4326) # back to degrees
        firstVar = True
        for DataVar in ['flux_apos_bio0','flux_apri_bio0','flux_apos_ocean0','flux_apri_ocean0','flux_foss0']:
            nonanindex = ~np.isnan(inter[DataVar])
            #calculate weighted mean on the 1x1° grid (using GridID)
            newdf2 = inter.loc[inter.index[nonanindex]].groupby(
                        ['GridID','Month','Year'])[DataVar].agg(
                        lambda x: np.average(
                        x, weights=inter.loc[x.index, "subarea"])).reset_index()
            if firstVar:
                newdf = newdf2.copy()
            else:
                newdf = newdf.merge(newdf2,on=['GridID','Month','Year'],how='outer')
            firstVar = False
        gdfGrid = gdfGrid.drop(columns='geometry')
        gdfGrid = gdfGrid.rename(columns={'geomPoint':'geometry'})
        gdfGrid = gdfGrid.set_geometry(gdfGrid.geometry)
        del gdf
        #get  Lat, Lon and geometry in the new gdf 'newdf'
        #gdf = pd.merge(gdfGrid[['GridID','geometry','Lat','Long','AreaGrid']],newdf,on=['GridID'])
        gdf = pd.merge(gdfGrid[['GridID','geometry','Lat','Long']],newdf,on=['GridID'])
        gdf = gdf.drop(columns = 'GridID')
        Area = getAreaOfGrid()
        gdf = pd.merge(gdf, Area, on=['Lat'])     
    
        for DataVar in ['flux_apos_bio','flux_apri_bio','flux_apos_ocean','flux_apri_ocean','flux_foss']:
            #gdf.insert(loc = 1, value = gdf[DataVar+'0']*gdf['AreaGrid'],column = DataVar+'_tot')
            gdf.insert(loc = 1, value = gdf[DataVar+'0']*gdf['Area'],column = DataVar+'_tot')
        
        date = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(15)),axis=1)
        gdf.insert(loc=1,column='Date',value=date)
        Monthdate = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                        int(x.Month),
                                                        15),axis=1)
        gdf.insert(loc=1,column='MonthDate',value=Monthdate)
        #this should change nothing as the 1x1° gdf is already on Australian region
        
        if Num >= 700: 
            #geodataframe (gdf) with borders of the Transcom regions
            if y == starty:
                if Num >= 900:
                    Transcom = pd.read_pickle("/Transcom_Regions.pkl")
                elif Num >= 700:
                    Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
                    
            igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
            gdf = gdf.loc[igdf]
       
        gdf.to_pickle("/CAMS/dataframes/GDF21x1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")

        del gdf,df
        
def CombineGDFCAMS(RegionName, Num,Datatype):
    if 'Sat' in Datatype:
        starty = 2014
    else:
        starty = 2009
    for y in range(starty,2020):
        df2 = pd.read_pickle("/CAMS/dataframes/GDF21x1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")
        if y == starty:
            df = df2.copy()
        else:           
            df = df.append(df2, ignore_index=True)
            
    df.to_pickle("/CAMS/dataframes/GDF21x1FLUX"+Datatype+"_"+RegionName+str(Num)+".pkl")

#main
if __name__=='__main__':
    Numm  = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Sat')   
    CombineGDFCAMS(RegionName, Numm,'Sat') 
