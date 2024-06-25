#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript to crate geopandas dataframe based on transcom region definition
out of the monthly 1x1° datasets of MIP
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import datetime
import numpy as np
import pandas as pd
from functions import getAreaOfGrid
import geopandas
import xarray as xr
from RegionParam import getRegion

def CreateDF(DS):
    for d in range(len(DS['land'])): #over time
           
        dfy = pd.DataFrame(data=DS.variables['land'][d].values,
                                index = DS.variables['latitude'].values, 
                                columns = DS.variables['longitude'].values, 
                                dtype='float' )
        dfy = pd.melt(dfy.reset_index(), id_vars='index',
                                    value_name ='land',var_name = 'Long')
        dfy.loc[:,"Lat"] = dfy['index']
        if 'start_date' in list(DS.keys()): #this is the case for all but LoFI
            dfy.insert(loc=1,column='Year',value= np.ones(len(dfy["Lat"]))*int(DS.start_date[d][0]))        
            dfy.insert(loc=1,column='Month',value= np.ones(len(dfy["Lat"]))*int(DS.start_date[d][1]))
        else:
            DS2 = xr.open_dataset("/OCO-2_v10_MIP/LNLGIS/CAMS_gridded_fluxes_LNLGIS.nc4")
            dfy.insert(loc=1,column='Year',value= np.ones(len(dfy["Lat"]))*int(DS2.start_date[d][0]))        
            dfy.insert(loc=1,column='Month',value= np.ones(len(dfy["Lat"]))*int(DS2.start_date[d][1]))
        dfy.insert(loc=1,column='Day',value= np.ones(len(dfy["Lat"]))*15)
        
            
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
        del dfy         
    return df

def CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,version, ModelNames):
    ModelNamesl = ['LoFI','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames','UT','WOMBAT','OU','CSU','NIES','CAMS']#,'LoFI',
    ModelNamesIS = [x+'_IS' for x in ModelNamesl]
    ModelNamesprior = [x+'_prior' for x in ModelNamesl]
    
    
    print("Start reading data")
    # from 20230727 on changed regular to regularOld and regularnew to regular
    if version == 'regular':
        filepath = "/OCO-2_v10_MIP/LNLGIS/EnsMean_gridded_fluxes_LNLGIS.nc4"
        filepathStd = "/OCO-2_v10_MIP/LNLGIS/EnsStd_gridded_fluxes_LNLGIS.nc4"
    elif version == 'comp':
        filepath = "/OCO-2_v10_MIP/MIP_v10_Sourishnew/gridded_fluxes_original_20221106/flux_mip/gridded_fluxes/EnsMean_gridded_fluxes_LNLGIS.nc4"
        filepathStd = "/OCO-2_v10_MIP/MIP_v10_Sourishnew/gridded_fluxes_original_20221106/flux_mip/gridded_fluxes/EnsStd_gridded_fluxes_LNLGIS.nc4"
    elif version == 'regularIS':
        filepath = "/OCO-2_v10_MIP/IS/EnsMean_gridded_fluxes_IS.nc4"
        filepathStd = "/OCO-2_v10_MIP/IS/EnsStd_gridded_fluxes_IS.nc4"
    elif version == 'TM5-4DVar':
        filepath = "/OCO-2_v10_MIP/LNLGIS/TM5-4DVAR_gridded_fluxes_LNLGIS.nc4"
    elif version in ModelNames:
        filepath = "/OCO-2_v10_MIP/LNLGIS/"+version+"_gridded_fluxes_LNLGIS.nc4"
    elif version in ModelNamesIS:
        filepath = "/OCO-2_v10_MIP/IS/"+version.split('_')[0]+"_gridded_fluxes_IS.nc4"
    elif version in ModelNamesprior:
        filepath = "/OCO-2_v10_MIP/Prior/"+version.split('_')[0]+"_gridded_fluxes_Prior.nc4"
    elif version == 'regularOld':
        filepath = "/OCO-2_v10_MIP/LNLGIS_new/EnsMean_gridded_fluxes_LNLGIS.nc4"
        filepathStd = "/OCO-2_v10_MIP/LNLGIS_new/EnsStd_gridded_fluxes_LNLGIS.nc4"
    else:
        print('version does not exist')

    DSMean = xr.open_dataset(filepath)
    df3Mean= CreateDF(DSMean)
    df3Mean.drop(columns='index',inplace = True)
    dfMean = df3Mean[(df3Mean.Long >= Long_min) & (df3Mean.Long <= Long_max) & (df3Mean.Lat >= Lat_min) & (df3Mean.Lat <= Lat_max)]
    
    if version not in ModelNames+ModelNamesIS+ModelNamesprior:
        DSStd = xr.open_dataset(filepathStd)
        df3Std= CreateDF(DSStd)
        df3Std.drop(columns='index',inplace = True)
        dfStd = df3Std[(df3Std.Long >= Long_min) & (df3Std.Long <= Long_max) & (df3Std.Lat >= Lat_min) & (df3Std.Lat <= Lat_max)]
        dfStd.rename(columns = {'land':'landStd'},inplace=True)
        df = pd.merge(dfMean,dfStd,on =['Year','Month','Day','Lat','Long'])
    else:
        df = dfMean.copy(deep=True)

    print("create timestamp")
    date = df.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(x.Day)),axis=1)
    df.insert(loc=1,column='Date',value=date)
    df.insert(loc=1,column='MonthDate',value=date)

    #gdf with latitude dependent area of a 1x1° cell
    Area = getAreaOfGrid()
    df = pd.merge(df, Area, on=['Lat'])     
    # former unit gC/(m^2*year)  -> TgC/(gridcell*month) by multiplying with * m²/gridcell * 1/12 *1//10**12
    landt = df.apply(lambda x: (x.land * 1/12 * x.Area/10**12),axis=1)
    df.insert(loc=1, column = 'Landtot', value = landt)
    print('works')
    if version not in ModelNames+ModelNamesIS+ModelNamesprior:
        
        landStdt = df.apply(lambda x: (x.landStd * 1/12 * x.Area/10**12),axis=1) 
        df.insert(loc=1, column = 'LandStdtot', value = landStdt)
    
    if version == 'regular':
        df.to_pickle("/OCO-2_v10_MIP/dataframes/DF2MIP_"+RegionName+str(Num)+".pkl")
    elif version == 'regularIS':
        df.to_pickle("/OCO-2_v10_MIP/dataframes/DF2MIP_IS_"+RegionName+str(Num)+".pkl")
    elif version == 'regularOld':
        df.to_pickle("/OCO-2_v10_MIP/dataframes/DF2MIP_old_"+RegionName+str(Num)+".pkl")
    else:
        df.to_pickle("/OCO-2_v10_MIP/MIP_v10_Sourishnew/dataframes/DF2MIP_"+RegionName+str(Num)+".pkl")
    
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
    if Num >= 700: 
        if Num >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        elif Num >= 700:
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
        
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]

    if version == 'regular':
        gdf.to_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+RegionName+str(Num)+".pkl")
    elif version == 'regularIS':
        gdf.to_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_IS_"+RegionName+str(Num)+".pkl")
    elif version == 'regularOld':
        gdf.to_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_old_"+RegionName+str(Num)+".pkl")
    elif version == 'TM5-4DVar':
        gdf.to_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_TM5-4DVar_"+RegionName+str(Num)+".pkl")
    elif version in ModelNames+ModelNamesIS+ModelNamesprior:
        gdf.to_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+version+"_"+RegionName+str(Num)+".pkl")
        
    else:
        df.to_pickle("/OCO-2_v10_MIP/MIP_v10_Sourishnew/dataframes/GDF2MIP_"+RegionName+str(Num)+".pkl")
    
    
    del gdf,df
    
#main
if __name__=='__main__':
    Numm  = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    ModelNames = ['LoFI','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames','UT','WOMBAT','OU','CSU','NIES','CAMS']#,'LoFI',
    ModelNamesIS = [x+'_IS' for x in ModelNames]
    ModelNamesprior = [x+'_prior' for x in ModelNames[1:]]#not available for LoFI
    for Mname in ModelNamesIS[8:]:
        print(Mname)
        CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Mname,ModelNames)    
