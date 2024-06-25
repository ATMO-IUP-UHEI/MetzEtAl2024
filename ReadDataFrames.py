#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Script to read dataframes and produce MonthMeanFrames, to be combined with Plot_TimeSeries
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

import glob, os, sys
import numpy as np
import pandas as pd
import geopandas
import logging
import datetime
import xarray as xr

from functions2 import getTimeMeans, DetrendMonthdate, getReferenceDate, getTimeSums, addMonthDate, getDataInfo
from CreateSingleDataFrame import CreateDF
from RegionParam import getRegion


def getData(Dataset ,ValueName, timeperiod,RegionName, Numm,FrameType, offset, minNum,rate='Noaa', **selection):
    '''
    # function to read gdfs and return dataframes with monthly mean values
    # arguments:
    #           Dataset: 'DataName_resolution' as given in /DataFrames/Data_Overview.csv
    #           ValueName: name of variable in gdf (needed for month aggregation), look in Data_Overview.csv for available vars
    #           timeperiod: list containing the start and end date of the chosen period as datetime.dates
    #           RegionName: name of the investigated region
    #           Numm: ID of the investiagted region
    #           offsets: offset which shall be used for detrending
    #           minNum: threshold for the monthly number of values
    #           FrameType: wether to return the Monthly Means ('MM') or the whole GeoDataFrame/xarray ('all')
    #           rate: which rate to use for detrending
    #           optional keywords: [newDFversion,DFversion,meas_geom,quality,commentNewDf ]
    # returns:
    #           MM: pandas DataFrame containing detrended monthly mean values
    '''
    dfPath = "/DataFrames/"
    DataName, DataRes, datafiles, DataDir, DFversion, Datainfo = getDataInfo(Dataset)
    commentNewDf = ''
    if 'DFversion' in selection.keys():
        DFversion = selection['DFversion'] #if older version requested set, DFversion to this old version
        del selection['DFversion']
  
    if 'newDFversion' in selection.keys():
        # check if newversion > than current versions
        if selection['newDFversion'] <= DFversion:
            sys.exit('A newDFversion for '+Dataset+' is requested in getData but is < the current DFversion, please increase the newDFversion or delete attribute.')
        if 'commentNewDf' in selection.keys():
            commentNewDf = selection['commentNewDf']
         
        DFversion = selection['newDFversion']
        del selection['newDFversion'], selection['commentNewDf'] 
        


    DateRef = getReferenceDate(timeperiod[0].year,timeperiod[1].year,timeperiod[0].month,timeperiod[1].month)
    
    #create dataframe of dataset if it doesn't extist
    if not os.path.isfile(dfPath +"DFv"+str(DFversion)+"_"+Dataset+"_"+RegionName+str(Numm)+".pkl") or 'newDFversion' in selection.keys():
        print('a new dataframe will be created for '+Dataset+"_"+RegionName+str(Numm))

        availkeys = CreateDataFrame(Numm,DataName,DataRes,DataDir,datafiles,DFversion,commentNewDf)
        #log new informations
        Datainfo.loc[(Datainfo.DataName == DataName)&(Datainfo.DataRes ==DataRes),['DFversion']] = DFversion
        Datainfo.loc[(Datainfo.DataName == DataName)&(Datainfo.DataRes ==DataRes),['AvailVars']] = str(availkeys)
        Datainfo.to_csv("/DataFrames/Data_Overview.csv")
    
    if FrameType == 'MS': #return monthly sum of fluxes
        #check if dataframe already exists
        availMM = os.path.isfile(dfPath + "MonthFrames/MonthSums_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+RegionName+str(Numm)+".pkl")
        if availMM:
            print('read MonthSums '+ Dataset )
            MS = pd.read_pickle(dfPath + "MonthFrames/MonthSums_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+RegionName+str(Numm)+".pkl")
            
        else:
            gdf = pd.read_pickle(dfPath+"DFv"+str(DFversion)+"_"+Dataset+'_'+RegionName+str(Numm)+".pkl")
            MS = DateRef.copy(deep=True)
            for valname in ValueName:
                MSadd = getTimeSums(gdf,valname,'Month')
                MSadd = addMonthDate(MSadd)
                MS = pd.merge(MS, MSadd, on=['MonthDate'], how = 'left')
            MS.to_pickle(dfPath + "MonthFrames/MonthSums_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+RegionName+str(Numm)+".pkl")
        MS = pd.merge(DateRef, MS, on=['MonthDate'], how = 'left')
        result = MS.copy(deep = True)
        
        
        del MS
        
    elif FrameType == 'MM': #return monthly mean values

        flags = '' # used for filename of monthly means
        for key in selection.keys():
            flags += key + '-'+str(selection[key])+'_'

        #check if dataframe already exists
        availMM = os.path.isfile(dfPath + "MonthFrames/MonthMeans_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+flags+RegionName+str(Numm)+".pkl")
        print(DFversion)
        availkeys = []
        if availMM:
            print('read MonthMeans '+ Dataset )
            MM = pd.read_pickle(dfPath + "MonthFrames/MonthMeans_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+flags+RegionName+str(Numm)+".pkl")
            availkeys = MM.keys()
        
        #check if needed variables are available
        if not ValueName+'_Detrend'+str(offset)+'_'+str(minNum)+'_rate'+rate in availkeys or not availMM:
            print('read gdf')
            gdf = pd.read_pickle(dfPath+"DFv"+str(DFversion)+"_"+Dataset+'_'+RegionName+str(Numm)+".pkl")

            #filter according to flags
            print('select according to flags')
            for key in selection.keys():
                gdf = gdf[(gdf[key] == selection[key])] #attribute selection according to keyword arguments (e.g. quality=0)          
            
            if ValueName+'_error' in gdf.keys():
                if ValueName+'_uncorr' in gdf.keys():
                    MMadd = getTimeMeans(gdf,ValueName,'Month',minNum,ValueName+'_error',ValueName+'_uncorr')
                else:
                    MMadd = getTimeMeans(gdf,ValueName,'Month',minNum,ValueName+'_error')
            else:
                MMadd = getTimeMeans(gdf,ValueName,'Month',minNum)
            MMaddD = DetrendMonthdate(MMadd, offset, minNum,ValueName,rate)

            if availMM: # monthly means already existed, but variable was missing
                #get all keys in new month means and add it to former month means
                detrends = [name for name in MMaddD.keys().values if 'Detrend' in name]
                detrends = detrends + [name for name in MMaddD.keys().values if (ValueName in name and name not in detrends and name not in MM.keys().values)]
                if 'level_0' in MM.keys().values:
                    MM.drop(columns=['level_0'], inplace=True)
                MM = pd.merge(MM, MMaddD[detrends+['Year','Month']], on = ['Year','Month']).reset_index()
            else:
                MM = MMaddD.copy()
            MM.to_pickle(dfPath + "MonthFrames/MonthMeans_DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+ValueName+"_"+flags+RegionName+str(Numm)+".pkl")
            del gdf
        
        MM = pd.merge(DateRef, MM, on=['MonthDate'], how = 'left')
        result = MM.copy(deep = True)
        del MM

    else: # get whole dataset, no month means
        print('read gdf')
        if DataRes == 'grid':
            gdf = xr.open_dataset(dfPath+"DF0v"+str(DFversion)+"_"+Dataset+'_'+RegionName+str(Numm)+".nc")
        else:
            gdf = pd.read_pickle(dfPath+"DFv"+str(DFversion)+"_"+Dataset+'_'+RegionName+str(Numm)+".pkl")
        result = gdf.copy(deep = True)
        del gdf

    return result





def CreateDataFrame(Numm,DataName,DataRes,DataDir,datafiles,DFversion,commentNewDf):
    '''
    # function to create and save a pandas Dataframe out of spatial data, successful create is logged
    # arguments:
    #           Numm: region ID, needed to get Region Parameter
    #           DataName: Name of dataset, used for df file name
    #           DataRes: Temporal resolution of dataset, choose from ['sounding','grid'], used for df file name  
    #           DataDir: directory to folder containing the data
    #           datafiles: file name, can contain wildcard to load multiple data (e.g. '*.nc' = all nc files in folder)
    #           DFversion: version of GDF as string or int
    '''   
    savepath = "Savepath"
    # Create and configure logger
    logging.basicConfig(filename = savepath + "DataFramesList.log",
                        format='%(asctime)s [%(levelname)s] %(message)s')
    logger = logging.getLogger('Read_DataFrames') # Creating an object with name Read_DataFrames
    logger.setLevel(logging.INFO) # Setting the threshold of logger to INFO

    
    tstart = datetime.datetime.now()
    print("Start reading data:")
    
    if DataRes == 'sounding':
        treadend, availVars = CreateDataFrameSounding(Numm,DataName,DataRes,DataDir,datafiles,DFversion, savepath)
    elif DataRes == 'grid':
        treadend, availVars = CreateDataFrameGrid(Numm,DataName,DataRes,DataDir,datafiles,DFversion, savepath)
    elif 'cs' in DataRes:
        treadend, availVars = CreateDataFrameSounding(Numm,DataName,DataRes,DataDir,datafiles,DFversion, savepath)
    else:
        sys.exit('Not yet implemented')
        
    dtime = datetime.datetime.now() - tstart
    dtimeread = treadend - tstart
    logger.info(', '+str(Numm)
                +', '+DataName
                +', '+DataRes
                +', '+DataDir+datafiles
                +', '+str(DFversion)
                +', '+str(dtime.total_seconds())
                +', '+str(dtimeread.total_seconds())
                +', Comment: '+commentNewDf
                )
    return availVars


def CreateDataFrameSounding(Numm,DataName,DataRes,DataDir,datafiles,DFversion, savepath):  
    '''
    # function to create and save a pandas Dataframe and GeoDataframe out of spatial vector data, successful create is logged
    # arguments:
    #           Numm: region ID, needed to get Region Parameter
    #           DataName: Name of dataset, used for df file name
    #           DataRes: Temporal resolution of dataset, choose from ['sounding','grid'], used for df file name  
    #           DataDir: directory to folder containing the data
    #           datafiles: file name, can contain wildcard to load multiple data (e.g. '*.nc' = all nc files in folder)
    #           DFversion: version of GDF as string or int
    #           savepath: where to save the dataframes
    # saves:
    #           DataFrame of all vector data within the coarse rectengular region
    #           Geodataframe of all vector data within the finer (for Numm > 900) region
    '''   
    #get region parameters
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    print(glob.glob(DataDir+datafiles))
    for num, FilePath in enumerate(glob.glob(DataDir+datafiles)):
        print(FilePath)
        #create Dataframe
        df3 = CreateDF(FilePath,Numm,DataName,DataRes)
        # cut out the coarse region
        df2 = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)].copy()  
        if num == 0:        
            df = df2.copy()            
        else:            
            df = pd.concat([df,df2], ignore_index=True)
     
    print("finished reading data") 
    treadend = datetime.datetime.now()
    
    #create date variable
    print("create timestamp")
    if 'Date' not in df.keys():
        date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
        df.insert(loc=1,column='Date',value=date)

    #create calender week number variable: deleted if needed check CreateAndSaveGeoDataframe.py in Skript folder 

    df.to_pickle(savepath+"DF0v"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+RegionName+str(Numm)+".pkl")

    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 'EPSG:4326')
    
    # all numbers >=700 are not rectangle but have a shapefile as border
    if Numm >= 700:
        if Numm >= 900:
            Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        else:
            print('read CT Transcom')
            Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")

        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf2 = gdf.loc[igdf]
        del gdf
        gdf = gdf2.copy(deep=True)
        print(len(gdf))

    gdf.to_pickle(savepath+"DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+RegionName+str(Numm)+".pkl")
    return treadend, list(df.keys())


def CreateDataFrameGrid(Numm,DataName,DataRes,DataDir,datafiles,DFversion, savepath):  
    '''
    # function to create and save a regional xarray on 1°x1° and a regional mean pandas Dataframe out of grid data, successful create is logged
    # arguments:
    #           Numm: region ID, needed to get Region Parameter
    #           DataName: Name of dataset, used for df file name
    #           DataRes: Temporal resolution of dataset, choose from ['sounding','grid'], used for df file name  
    #           DataDir: directory to folder containing the data
    #           datafiles: file name, can contain wildcard to load multiple data (e.g. '*.nc' = all nc files in folder)
    #           DFversion: version of GDF as string or int
    #           savepath: where to save the dataframes
    # saves:
    #           is possible xarray dataset on 1x1 cut to the region
    #           Dataframe with regional means but daily time resolution
    '''   
    #get region parameters
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
         
    for num, FilePath in enumerate(glob.glob(DataDir+datafiles)):
        #print("read: "+FilePath)
        #create xarray dataset for region, it needs:
        #  to be on 1x1, 
        #  CO2 in ppm,
        #  time as string YYYYMMDD or datetime64 
        #  daily means or coarser
        DS2 = CreateDF(FilePath,Numm,DataName,DataRes)
        if num == 0:        
            DS = DS2.copy(deep = True)            
        else:            
            DS = xr.concat([DS,DS2], dim = 'time')
 
    DS.to_netcdf(savepath+"DS0v"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+RegionName+str(Numm)+".nc")
    
    print("finished reading data") 
    treadend = datetime.datetime.now()
    
    #create df with mean region concentration or sum of regional fluxes
    if 'flux' in DataName:
        DS_mean = DS.sum(dim=['latitude','longitude'])
        df = pd.DataFrame(data = {'time':DS_mean.time.values,'CO2':DS_mean.CO2_total.values})
    else:
        DS_mean = DS.mean(dim=['latitude','longitude'])
        df = pd.DataFrame(data = {'time':DS_mean.time.values,'CO2':DS_mean.CO2.values})
    
    #create date variable, set all dates to middle of month (day = 15)!
    print("create timestamp")
    df.insert(loc = 1, column = 'Year', value = df.apply(lambda x: int(str(x.time).replace('-','')[0:4]),axis = 1))#replace needed as time can be numpy datetime64 or string
    df.insert(loc = 1, column = 'Month', value = df.apply(lambda x: int(str(x.time).replace('-','')[4:6]),axis = 1))
    #df.insert(loc = 1, column = 'Day', value = np.ones(len(df.Year)).astype(int)*15)
    df.insert(loc = 1, column = 'Day', value = df.apply(lambda x: int(str(x.time).replace('-','')[6:8]),axis = 1))
    
    date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
    df.insert(loc=1,column='Date',value=date)
    df.drop(columns= ['time'],inplace = True)

    
    df.to_pickle(savepath+"DFv"+str(DFversion)+"_"+DataName+"_"+DataRes+"_"+RegionName+str(Numm)+".pkl")
    return treadend, list(DS.keys()) + list(DS.coords)
        