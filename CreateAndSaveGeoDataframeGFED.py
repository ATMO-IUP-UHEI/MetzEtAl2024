#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""
from calendar import monthcalendar
import h5py 
import pandas as pd
import numpy as np
import datetime
import geopandas
import os

def createAnnualDataFrameGFED(filepath,datayear,datatype):
    '''
    # function to create a pandas dataframe out of one GFED4 file
    # input:
    #           filepath: file of GFED4 data
    #           datayear: year of the datafile
    #           datatype: 0=carbon emissions, 1=CO2 emissions, 2=CO emissions
    #
    # output:
    #           df: pandas DataFrame containing the variables
    #           -- 'emission',
    #           -- 'Month' integer
    #           -- 'Date' datetime object with day = 15
    #           -- 'Year' integer
    # Documentation GFED see https://www.geo.vu.nl/~gwerf/GFED/GFED4/Readme.pdf
    '''
    DS = h5py.File(filepath,'r')
    typeL = ['SAVA','BORF','TEMF','DEFO','PEAT','AGRI']
    CO2L = [1686,1489,1647,1643,1703,1585] #emission factor from https://www.geo.vu.nl/~gwerf/GFED/GFED4/ancill/GFED4_Emission_Factors.txt
    COL = [63,127,88,93,210,102] # emission factors see above
    
    #get grid cell areas
    dfy_gridarea = pd.DataFrame(data=DS['ancill/grid_cell_area'][:],index = DS['lat'][:,0],
                                columns = DS['lon'][0,:], dtype='float' )
    dfy_gridarea = pd.melt(dfy_gridarea.reset_index(), id_vars='index',
                                value_name ='Grid_area',var_name = 'Long')
    dfy_gridarea["Lat"] = dfy_gridarea['index']

    #store the values for each month
    for m in range(1,13):
        mo = str(m)
        print(m)
        if datatype == 0:
            dfy = pd.DataFrame(data=DS['emissions/'+mo.zfill(2)+'/C'][:],index = DS['lat'][:,0],
                                  columns = DS['lon'][0,:], dtype='float' )
            if datayear < 2017:
                dfy_area = pd.DataFrame(data=DS['burned_area/'+mo.zfill(2)+'/burned_fraction'][:],index = DS['lat'][:,0],
                                  columns = DS['lon'][0,:], dtype='float' )
            else: #cburned fraction only given in this format until 2016, add nans for later years
                arr = np.empty([720,1440])
                arr[:] = np.nan
                dfy_area = pd.DataFrame(data=arr,index = DS['lat'][:,0],
                                  columns = DS['lon'][0,:], dtype='float' )
            dfy_area = pd.melt(dfy_area.reset_index(), id_vars='index',
                                  value_name ='burned_fraction',var_name = 'Long')
            dfy_area["Lat"] = dfy_area['index']

        elif datatype == 1:
            datae = DS['emissions/'+mo.zfill(2)+'/DM'][:]
            CO2emission = np.zeros((720, 1440))
            for t in range(6):#for each of the eco systems, add up emissions
                contribution = DS['emissions/'+mo.zfill(2)+'/partitioning/DM_'+typeL[t]][:]
                CO2emission += datae * contribution * CO2L[t]
            #Create Pandas DataFrame
            dfy = pd.DataFrame(data=CO2emission,index = DS['lat'][:,0],
                                  columns = DS['lon'][0,:], dtype='float' )
       
        elif datatype == 2:
            datae = DS['emissions/'+mo.zfill(2)+'/DM'][:]
            CO2emission = np.zeros((720, 1440))
            for t in range(6):#for each of the eco systems, add up emissions
                contribution = DS['emissions/'+mo.zfill(2)+'/partitioning/DM_'+typeL[t]][:]
                CO2emission += datae * contribution * COL[t]
            #Create Pandas DataFrame
            dfy = pd.DataFrame(data=CO2emission,index = DS['lat'][:,0],
                                  columns = DS['lon'][0,:], dtype='float' )
        #unpivot dataframe (2D array -> columns: Lat, Long, emission)
        dfy = pd.melt(dfy.reset_index(), id_vars='index',
                                  value_name ='emission',var_name = 'Long')
        dfy["Lat"] = dfy['index']

        #merge burned area if applicable
        if datatype == 0:
            dfy = pd.merge(dfy, dfy_area, on=['Lat','Long'], how = 'outer')

        #merge grid area into df and calculate total emissions
        dfy = pd.merge(dfy, dfy_gridarea, on=['Lat','Long'], how = 'outer')
        dfy.insert(loc=1, column='total_emission', value = dfy.Grid_area*dfy['emission']/(1000000000000))

        dfy.insert(loc=1,column='Month',value= int(m)*np.ones(len(dfy['emission'])))
        dfy.insert(loc=1,column='Year',value= int(datayear)*np.ones(len(dfy['emission'])))
        MonthDate = dfy.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
        dfy.insert(loc=1,column='Date',value=MonthDate)

        if m ==1:
            df = dfy.copy()
        else:
            df = df.append(dfy)#merge dataframes
    

    return df


def createDataFrameGFED(datapath,savepath,datatype,Long_min,Long_max,Lat_min,Lat_max,RegionName,Num, Transcom):
    typelist = ['C','CO2','CO']
    if not os.path.isfile(savepath+"DFglobalnew_"+typelist[datatype]+".pkl"):
        print('start creating Dataframe')
        print('process 2009')
        df = createAnnualDataFrameGFED(datapath+'GFED4.1s_2009.hdf5',2009,datatype)
        for j in range(2010,2024):
            print('process '+str(j))
            if j < 2017:
                df2 = createAnnualDataFrameGFED(datapath+'GFED4.1s_'+str(j)+'.hdf5',j,datatype)
            else:
                df2 = createAnnualDataFrameGFED(datapath+'GFED4.1s_'+str(j)+'_beta.hdf5',j,datatype)
            df = df.append(df2)
        df.to_pickle(savepath+"DFglobalnew_"+typelist[datatype]+".pkl")
    else:
        print('reading dataframe, this may take some minutes')
        df = pd.read_pickle(savepath+"DFglobalnew_"+typelist[datatype]+".pkl")
    df = df[(df.Long >= Long_min) & (df.Long <= Long_max) & (df.Lat >= Lat_min) & (df.Lat <= Lat_max)]
    print('create GeoDataFrame')
    gdf = geopandas.GeoDataFrame(
                    df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 'EPSG:4326')

    if UseTranscom:
        print('cut Transcom Region, this too needs a cuple of minutes')
        if not os.path.isfile(savepath+"GDFnew_"+typelist[datatype]+"_"+RegionName+".pkl"):
            if Num >= 700: 
                #geodataframe (gdf) with borders of the Transcom regions
                if Num >= 900:
                    Transcom = pd.read_pickle("/Transcom_Regions.pkl")
                elif Num >= 700:
                    Transcom = pd.read_pickle("/CTTranscomMask1x1Borders.pkl")
                
                #not needed for the interpolated gdf as gdfGrid already cut by Transcom border
                igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
                gdf = gdf.loc[igdf]
    gdf.to_pickle(savepath+"GDFnew_"+typelist[datatype]+"_"+RegionName+".pkl")

    return gdf

def getMonthSum(gdf,date_min,date_max,ValueName):
    # function to create monthly sum of indicated column inside given borders
    # input:
    #       gdf: Geodataframe containing the data
    #       year_max, month_max, day_max: enddate
    #       year_min,month_min, day_min: startdate
    #       Long_min, Long_max: Longitude range
    #       Lat_min, Lat_max: Latitude range
    #       ValueName: Name of column to sum up
    #
    # output:
    #       df with monthly summed values
    print('calculate monthly means')
    gdf2 = gdf[(gdf.Date >= date_min) & (gdf.Date <= date_max)]

    output = gdf2.groupby(['Year','Month'])[ValueName].sum().reset_index()
    MonthDate = output.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    output.insert(loc=1,column='MonthDate',value=MonthDate)

    return output



#main
if __name__ == "__main__": #if this script is executed, not executed if function of this script is imported in another script
    datapath = "/mnt/data/users/eschoema/GFED4/"
    savepath = "/mnt/data/users/eschoema/GFED4/"
    datatype = 2 #0=C, 1=CO2, 2=CO
    UseTranscom = False
    Numm = 756
    RegionName, Long_min, Long_max, Lat_min, Lat_max = 'AU',112,180,-50,-10, #getRegion(Numm)

    typelist = ['C','CO2','CO']
    
    if not os.path.isfile(savepath+"GDFnew_"+typelist[datatype]+"_"+RegionName+".pkl"):
        gdf = createDataFrameGFED(savepath,savepath,datatype,Long_min,Long_max,Lat_min,Lat_max,RegionName,Numm, UseTranscom)
    else:
        print('reading GeoDataFrame, this will take some time')
        gdf = pd.read_pickle(savepath+"GDFnew_"+typelist[datatype]+"_"+RegionName+".pkl")
