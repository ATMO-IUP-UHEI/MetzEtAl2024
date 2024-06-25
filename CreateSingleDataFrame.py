#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Function to create a pandas dataframe from a single file, collection of different treatments
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz


import read_remotec_out
import numpy as np
import pandas as pd
import sys
from RegionParam import getRegion
from RegridXarray import RegridXarray
import xarray as xr
import datetime, os
from functions import CreateMask, get_date_Components, SelectTimeAndRegion, getNumDayOfMonth
import glob
import h5py


#main
def CreateDF(FilePath,Numm,DataName,DataRes):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           Numm: Number of selected region
    #           DataName: Name of the dataset to choose routine to extract variables from file
    #           DataRes: Resolution of the dataset to choose routine to extract variables from file
    #
    # returns:
    #           df: pandas DataFrame containing varaibles from the .out file
    '''
    if DataName == 'RT240' and DataRes == 'sounding':
        df = CreateDFRT240sounding(FilePath,Numm)
    elif DataName == 'ACOS' and DataRes == 'sounding':
        df = CreateDFACOSsounding(FilePath,Numm)
    elif DataName == 'OCO2v11' and DataRes == 'sounding':
        df = CreateDFOCO2v11sounding(FilePath,Numm)
    elif DataName == 'CT2022' and DataRes == 'grid':
        df = CreateDFCT2019Bgrid(FilePath,Numm)
    elif ('CT2022' in DataName or 'CAMS-IS2' in DataName) and 'cs' in DataRes:
        #get dataset to sample on
        COsampleOn = DataRes[2:]
        df = CreateDFCT2022cs(FilePath,DataName,Numm,COsampleOn)
    elif 'TM54DVar' in DataName and 'cs' in DataRes:
        #get dataset to sample on
        COsampleOn = DataRes[2:]
        df = CreateDFTM54DVarcs(FilePath,Numm,COsampleOn)
    elif (DataName == 'TM54DVarIS' or DataName == 'TM54DVarRemoTeC' or DataName == 'TM54DVarRemoTeC+ISloc') and DataRes == 'grid':
        df = CreateDFTM54DVargrid(FilePath,Numm)
    elif 'COCCON' in DataName and DataRes == 'sounding':
        df = CreateCOCCON(FilePath,Numm,DataName)
    elif 'CAMS' in DataName and DataRes == 'grid':
        df = CreateDFCAMSgrid(FilePath,Numm)
    elif 'CAMS' in DataName and 'cs' in DataRes:
        COsampleOn = DataRes[2:]
        df = CreateDFCAMScs(FilePath,DataName,Numm,COsampleOn)
    elif 'MIP' in DataName and 'cs' in DataRes:
        df = CreateMIPcsOCO2(FilePath,Numm, DataName)
    
    else:
        sys.exit('There is no function found in function CreateDF for the given dataset '+DataName+' and resolution '+DataRes)
        
    return df

#-----------------------
# individual functions
#-----------------------


def CreateDFRT240sounding(FilePath,Numm):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    # returns:
    #           df: pandas DataFrame containing varaibles from the .out file
    '''
    dic = read_remotec_out.read_output_file(FilePath)
            
    d = {'CO2': dic['REMOTEC-OUTPUT_X_COLUMN_CORR'].T[0], 
        'CH4': dic['REMOTEC-OUTPUT_X_COLUMN_CORR'].T[1],
        'CO2_uncorr': dic['REMOTEC-OUTPUT_X_COLUMN'].T[0],
        'CH4_uncorr': dic['REMOTEC-OUTPUT_X_COLUMN'].T[1],
        'Sec': dic['GOSAT_TIME_SEC'],
        'Min': dic['GOSAT_TIME_MIN'], 
        'Hour': dic['GOSAT_TIME_HOUR'],
        'Day': dic['GOSAT_TIME_DAY'], 
        'Month': dic['GOSAT_TIME_MONTH'],
        'Year': dic['GOSAT_TIME_YEAR'],
        'Lat':dic['GOSAT_LATITUDE'],
        'Long':dic['GOSAT_LONGITUDE'],
        'ElevationMin':dic['METEO_SURFACE_ELEVATION_MIN'],
        'ElevationMax':dic['METEO_SURFACE_ELEVATION_MAX'],
        'Height':dic['METEO_HEIGHT'][:,12],
        'sfp':dic['METEO_PRESSURE'][:,12],
        'CT_CO2': dic['REMOTEC-OUTPUT_X_APR_COLUMN'].T[0],
        'CT_CH4': dic['REMOTEC-OUTPUT_X_APR_COLUMN'].T[1],
        'CT_error': dic['REMOTEC-OUTPUT_X_APR_COLUMN_ERR'].T[0],
        'CO2_error': dic['REMOTEC-OUTPUT_X_COLUMN_ERR'].T[0],
        'meas_geom0' : dic['(0=NADIR,1=GLINT)'],
        'quality' : dic['REMOTEC-OUTPUT_FLAG_QUALITY'], #NICHT VERWENDEN, DA NICHT ZUVERLÄSSIG
        'gain': dic["GOSAT_GAIN"]}
    df = pd.DataFrame(data=d)
    
    df["Lat_round"] = np.floor(np.array(df["Lat"])) 
    df["Long_round"] = np.floor(np.array(df["Long"]))
    df.insert(loc=1,column='meas_geom',value=df.apply(lambda x: int(x.meas_geom0),axis=1))

    return df

def CreateDFACOSsounding(FilePath,Numm):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    # returns:
    #           df: pandas DataFrame containing varaibles from the .nc file
    '''
    
    DS = xr.open_mfdataset(FilePath,combine='by_coords',concat_dim='None',use_cftime=None)
    DS2 = xr.open_mfdataset(FilePath,group = 'Retrieval',combine='by_coords',concat_dim='None',use_cftime=None)
    DS3 = xr.open_mfdataset(FilePath,group = 'Sounding',combine='by_coords',concat_dim='None',use_cftime=None)
    year = int(str(DS.time.values[0])[0:4])
    month =int(str(DS.time.values[0])[5:7])
    day_min =int(str(DS.time.values[0])[8:10])
    day_max =int(str(DS.time.values[-1])[8:10])
    
    d = {'CO2': DS.xco2.values, 
         'CO2_error': DS.xco2_uncertainty.values, 
         'Lat':DS.latitude.values,
         'Long':DS.longitude.values,
         'Year':DS.date.values[:,0],
         'Month':DS.date.values[:,1], 
         'Day':DS.date.values[:,2],
         'Hour':DS.date.values[:,3], 
         'Min':DS.date.values[:,4],
         'Sec':DS.date.values[:,5],
         'CO2_uncorr': DS2.xco2_raw.values,
         'quality': DS.xco2_quality_flag.values,
         'glint': DS3.glint_angle.values,
         'land_frac': DS3.land_fraction,
         'DWS':DS2.dws.values,
         'delGrad':DS2.co2_grad_del.values,
         'psurf':DS2.psurf.values,
         'dpfrac':DS2.dpfrac.values,
         'albedo_sco2':DS2.albedo_sco2.values, 
         'gain':DS3.gain.values}
    df = pd.DataFrame(data=d)

    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    #IDENTIFY ACOS measurements also included in RemoTeCv2.4.0
    #get newest version of RT240_sounding
    Datainfo = pd.read_csv("/DataFrames/Data_Overview.csv")
    DFversionRT = int(Datainfo[(Datainfo.DataName == 'RT240')&(Datainfo.DataRes =='sounding')].DFversion.values[0])
    if not os.path.isfile("/DataFrames/DFv"+str(DFversionRT)+"_RT240_sounding_"+RegionName+str(Numm)+".pkl"):
        sys.exit('Please create a Dataframe for RT240 for the selected Region '+str(Numm))    
    gdf_R = pd.read_pickle("/DataFrames/DFv"+str(DFversionRT)+"_RT240_sounding_"+RegionName+str(Numm)+".pkl")
    gdf_R = SelectTimeAndRegion(gdf_R, year, month, day_min, year, month, day_max, Long_min, Long_max, Lat_min, Lat_max)
    
    if len(gdf_R) > 0:
        ACOSRemotecDate = df.apply(lambda x: datetime.datetime(int(x.Year),int(x.Month),int(x.Day),int(x.Hour),int(x.Min),int(x.Sec)),axis=1)
        df.insert(loc=1,column="RemotecDate",value=ACOSRemotecDate)
        #RT seconds are shifted by 2 seconds (no values > 2) and contain the value 60  no values >2
        modSec = np.round(gdf_R.Sec.values)-2 #correct seconds' shift 
        Iwrongsec = np.where(modSec == 60)
        modSec[Iwrongsec] = 0
        gdf_R.insert(loc=1,column="modSec",value=modSec)
        RemotecDate = gdf_R.apply(lambda x: datetime.datetime(int(x.Year),int(x.Month),int(x.Day),int(x.Hour),int(x.Min),int(x.modSec)),axis=1)
        RemotecDate.iloc[Iwrongsec] = RemotecDate.iloc[Iwrongsec].apply(lambda x: x+datetime.timedelta(0,0,0,0,1,0)) #add one minute to date former containing sec=60
        gdf_R.insert(loc=1,column="RemotecDate",value=RemotecDate)
        
        inRemotec = np.isin(df.RemotecDate,gdf_R.RemotecDate)
        gdf_R.rename(columns = {'CO2':'RT_CO2'}, inplace = True)
        df = pd.merge(df, gdf_R[['RT_CO2','RemotecDate']], on = 'RemotecDate', how = 'left' )
        df.drop(columns = 'RemotecDate',inplace=True)
    else:
        inRemotec = np.full(len(df),False)
        df.insert(loc=1,column="RT_CO2",value=np.full(len(df),np.nan))
        print('no RT data at '+str(DS.time.values[0]))

    df.insert(loc=1,column="in_Remotec_data",value=inRemotec)
    
    return df

def CreateDFOCO2v11sounding(FilePath,Numm):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    # returns:
    #           df: pandas DataFrame containing varaibles from the .nc file
    '''
    DS = xr.open_mfdataset(FilePath,combine='by_coords',concat_dim='None',use_cftime=None)
    DS2 = xr.open_mfdataset(FilePath,group = 'Retrieval', combine='by_coords',concat_dim='None')
    DS3 = xr.open_mfdataset(FilePath,group = 'Sounding', combine='by_coords',concat_dim='None')
        
    d = {'CO2': DS.xco2.values, 
         'CO2x2019': DS.xco2_x2019.values, 
         'CO2_uncorr':DS2.xco2_raw.values,
         'CO2_error': DS.xco2_uncertainty.values, 
         'Lat':DS.latitude.values,
         'Long':DS.longitude.values,
         'Year':DS.date.values[:,0],
         'Month':DS.date.values[:,1], 
         'Day':DS.date.values[:,2],
         'Hour':DS.date.values[:,3],
         'Minute':DS.date.values[:,4],
         'quality':DS.xco2_quality_flag.values, 
         'glint':DS3.glint_angle.values,
         'land':DS3.land_water_indicator.values,
         'sounding_id':DS.sounding_id.values}
         #'surface':np.abs(DS2.surface_type_flipped.values * DS3.land_water_indicator.values)}#for old 9r version
    df = pd.DataFrame(data=d)

    return df

def CreateDFCT2019Bgrid(FilePath,Numm):
    '''
    # function to create a xarray dataset for a certain region out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    # returns:
    #           xaray dataset: dataset containing varaibles from the .nc file
    '''
    DS = xr.open_dataset(FilePath)
    
    Molfrac = DS.co2[0] * DS.air_mass[0] 
    XCO2_0 = Molfrac.sum(dim = 'level')/DS.air_mass[0].sum(dim = 'level')
    
    DS_reg = InterpolateOnMaskXarray(XCO2_0, Numm,2, 3, 'CO2', 'nearest')
    
    return DS_reg

def CreateDFCT2022cs(FilePath,DataName,Numm,COsampleOn):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data cosampled on a second dataset
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    #           COsampleOn: dataset on which to be cosampled
    # returns:
    #           df: pandas DataFrame containing varaibles from the .nc file
    '''
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    
    if (DataName == 'CT2022' or DataName == 'CT2019B2' or DataName == 'CAMS-IS2') and COsampleOn == 'RT240':
        print('reading CT2022 cs RT240')
        if os.path.isfile(FilePath):# + RegionName + str(Numm)+'.pkl'):
            df0 = pd.read_pickle(FilePath)# + RegionName + str(Numm)+'.pkl')
        else:
            sys.exit('Cosampling of CT2022 on RT240 for region ' + str(Numm) +' not done yet. Do so by using SampleGridOnSatellite.py')
            
        df = df0.reset_index()
        df.insert(1, column = 'meas_geom', value = np.array(list(map(int,df.RT_mode.values))))
        df.drop(columns = ['index'],inplace = True)
        df["Lat_round"] = np.floor(np.array(df["Lat"])) 
        df["Long_round"] = np.floor(np.array(df["Long"]))

    
    else:
        sys.exit('Not implemented yet')
    return df

def CreateDFTM54DVargrid(FilePath,Numm):
    '''
    # function to create a xarray dataset for a certain region out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    # returns:
    #           xaray dataset: dataset containing varaibles from the .nc file
    '''
    #one day per file
    Datestr = FilePath.split('/')[-1].split('_')[1].split('.')[0]
    Year = int(Datestr[0:4]) 
    Month = int(Datestr[4:6]) 
    Day = int(Datestr[6:8]) 

    DS = xr.open_dataset(FilePath)
    DS2 = xr.open_dataset(FilePath,group='glb300x200')
        
    da = DS.at[:-1]-DS.at[1:]
    db = DS.bt[:-1]-DS.bt[1:]

    Airmass = DS2.pressure*db+((DS2.pressure/DS2.pressure)*da)#[times, latitude, longitude, levels]
    Airmass = Airmass.rename({'boundaries':'levels'})
    Airmass = Airmass.transpose('times', 'levels', 'latitude', 'longitude')
    Molfrac = DS2.mix[0] * Airmass #DS2.mix dims: tracers: 1, times: 8, levels: 25, latitude: 90, longitude: 120
    XCO2_0 = Molfrac.sum(dim = 'levels')/Airmass.sum(dim = 'levels')

    #mean over whole day
    XCO2_0meant = XCO2_0.mean(dim = 'times')
    XCO2_0meant = XCO2_0meant.expand_dims(dim={'time': [str(Year)+str(Month).zfill(2)+str(Day).zfill(2)]})#[datetime.date(Year,Month,Day)]})
    
    DS_reg = InterpolateOnMaskXarray(XCO2_0meant, Numm,2, 3, 'CO2', 'nearest')

    return DS_reg

def CreateDFTM54DVarcs(FilePath,Numm,COsampleOn):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data cosampled on a second dataset
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    #           COsampleOn: dataset on which to be cosampled
    # returns:
    #           df: pandas DataFrame containing varaibles from the .nc file
    '''
    if COsampleOn == 'ACOS' or COsampleOn == 'RT240':
        try: 
            if COsampleOn == 'ACOS':
                DS = xr.open_dataset(FilePath,group = 'CO2/GOSAT_ACOS',decode_times=False)
            else:
                DS = xr.open_dataset(FilePath,group = 'CO2/GOSAT_RemoTeC',decode_times=False)
        
            d = {'CO2': DS.modeled_XCO2.values, 
                'CO2_error': DS.sigma_XCO2.values, 
                'Lat':DS.latitude.values,
                'Long':DS.longitude.values,
                'decimal_date':DS.decimal_date.values,
                'CO2_uncorr': DS.modeled_raw.values, 
                'gain':DS.gain.values, 
                'quality':DS.xco2_quality_flag.values}

        except:
            d = {'CO2': [], 
                'CO2_error': [], 
                'Lat':[],
                'Long':[],
                'decimal_date':[],
                'CO2_uncorr': [], 
                'gain':[], 
                'quality':[]}
        df = pd.DataFrame(data=d)

        date = []
        Minute = []
        Month = []
        Day = []
        Hour = []
        Sec = []
        Year = []
    
        for i in df.decimal_date:
            #convert decimal year in datetime object
            yea = int(i)
            rem = i - yea
            
            datesec = datetime.datetime(yea,1,1) + datetime.timedelta(seconds=(datetime.datetime(yea+ 1,1,1)-datetime.datetime(yea,1,1)).total_seconds() * rem)        

            # get compounds of the datetime object
            Year.append(datesec.year)
            Month.append(datesec.month)
            Day.append(datesec.day)
            Hour.append(datesec.hour)
            Minute.append(datesec.minute)
            Sec.append(datesec.second)           
        df.insert(loc=1,column='Month',value=Month)
        df.insert(loc=1,column='Day',value=Day)
        df.insert(loc=1,column='Hour',value=Hour)
        df.insert(loc=1,column='Min',value=Minute)
        df.insert(loc=1,column='Year',value=Year)
        df.insert(loc=1,column='Sec',value=Sec)

    else:
        print('cosampling not implemented yet')
    return df

def CreateCOCCON(FilePath,Numm,DataName):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    # returns:
    #           df: pandas DataFrame containing varaibles from the .out file
    '''
    DS = h5py.File(FilePath, 'r')
    d = {'CO2': DS['/CO2.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR'][:], 
        'CO2_error': DS['/CO2.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR_UNCERTAINTY.RANDOM.STANDARD'][:],
        'CO': DS['/CO.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR'][:], 
        'CO_error': DS['/CO.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR_UNCERTAINTY.RANDOM.STANDARD'][:],
        'CH4': DS['/CH4.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR'][:], 
        'CH4_error': DS['/CH4.COLUMN.MIXING.RATIO.VOLUME.DRY_ABSORPTION.SOLAR_UNCERTAINTY.RANDOM.STANDARD'][:],
        'time':DS['/DATETIME'][:]}
    timeformat = 'fracDays2000'
    
    df = pd.DataFrame(data=d)
    
    if 'Gobabeb' in DataName:
        LongCOCCON, LatCOCCON = 15.0414, -23.5611
    elif 'Jinja' in DataName:
        LongCOCCON, LatCOCCON = 33.211, 0.421

    df.insert(loc = 1, column = 'Long',value = LongCOCCON)
    df.insert(loc = 1, column = 'Lat',value = LatCOCCON)
    df = get_date_Components(df,'time',True,timeformat)

    return df

def CreateDFCAMSgrid(FilePath,Numm):
    '''
    # function to create a xarray dataset for a certain region out of a single file containing spatial data
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    # returns:
    #           xaray dataset: dataset containing varaibles from the .nc file
    '''
    DS = xr.open_dataset(FilePath)
    year_str = str(DS.time.values[0]).split('-')[0] #extract year from first timestep
    month_str = str(DS.time.values[0]).split('-')[1] #extract month from first timestep
    
    for day in range(1,getNumDayOfMonth(int(year_str),int(month_str))+1):
        DS2 = DS.XCO2.loc[year_str+"-"+month_str+"-"+str(day).zfill(2):year_str+"-"+month_str+"-"+str(day).zfill(2)]
        if len(DS2) > 0:
            #mean over whole day
            XCO2_0meant = DS2.mean(dim = 'time')
            XCO2_0meant = XCO2_0meant.expand_dims(dim={'time': [year_str+month_str+str(day).zfill(2)]})

            if day == 1:        
                DSgrid = XCO2_0meant.copy(deep = True)            
            else:            
                DSgrid = xr.concat([DSgrid,XCO2_0meant], dim = 'time') 
    DSgridlat = DSgrid.reindex(latitude=DSgrid.latitude[::-1]) #reverse latitude axis from 90 - -90 to -90 - 90
    DSgridppm = DSgridlat * 1000000 # XCO2 in ppm
    DS_reg = InterpolateOnMaskXarray(DSgridppm, Numm,2, 4, 'CO2', 'linear')
    
    return DS_reg

def CreateDFCAMScs(FilePath,DataName,Numm,COsampleOn):
    '''
    # function to create a pandas dataframe out of a single file containing spatial data cosampled on a second dataset
    # arguments:
    #           FilePath: Filepath to single file with data
    #           NUmm: Region number
    #           COsampleOn: dataset on which to be cosampled
    # returns:
    #           df: pandas DataFrame containing varaibles from the .nc file
    '''
    if DataName == 'CAMS-IS' and COsampleOn == 'RT240':
        df = pd.read_pickle(FilePath)
    return df

def CreateMIPcsOCO2(FilePath, Numm, DataName):
    DS = xr.open_dataset(FilePath)
    experiment = DataName[3:]
    df0 = DS.to_dataframe().reset_index()
    df0.drop(columns = ['record'],inplace = True)
    df0.rename(columns = {'AMES':'AMES_CO2', 'BAKER':'BAKER_CO2', 'CAMS':'CAMS_CO2','CMS-FLUX':'CMS-FLUX_CO2', 'COLA':'COLA_CO2', 'CT':'CT_CO2', 'NIES':'NIES_CO2', 'OU':'OU_CO2', 'TM5-4DVAR':'TM5-4DVAR_CO2', 'UT':'UT_CO2', 'WOMBAT':'WOMBAT_CO2','WEIR':'WEIR_CO2',
                         'latitude':'Lat', 'longitude':'Long','obs_xco2':'OCO2_CO2'},inplace = True)
    RegionName = getRegion(Numm)[0]
    
    #read OCO-2 dataframe
    #DataNameOCO, DataResOCO, datafilesOCO, DataDirOCO, DFversionOCO, DatainfoOCO = getDataInfo('OCO2_sounding')
    df0.insert(0,column = 'Year',value = df0.apply(lambda  x: int(str(x.sounding_id)[0:4]), axis = 1))
    df0.insert(0,column = 'Month',value = df0.apply(lambda  x: int(str(x.sounding_id)[4:6]), axis = 1))
    df0.insert(0,column = 'Day',value = df0.apply(lambda  x: int(str(x.sounding_id)[6:8]), axis = 1))
    df0.insert(0,column = 'Hour',value = df0.apply(lambda  x: int(str(x.sounding_id)[8:10]), axis = 1))
    
    #replace CMs-Flux by correct values, the pandas dataframe was created by using ReadMATLABfile.py
    dfCMS = pd.read_pickle('/OCO-2_v10_MIP/OCO-2_Cosamples_Abhishek/DataFrames/CMS-Flux_'+experiment+'.pkl')
    df = pd.merge(df0,dfCMS,on = ['Year','Month','Day','sounding_id'],how = 'left')
    df.drop(columns = ['CMS-FLUX_CO2','observed xco2'], inplace = True)
    df.rename(columns={'modeled xco2':'CMS-FLUX_CO2'},inplace = True)
    
    #add JHU, which was provided seperately
    DSJHU = xr.open_dataset("/OCO-2_v10_MIP/OCO-2_Cosamples_Abhishek/JHU CoSamples/OCO2_v10mip.OCO2residuals."+experiment+".JHU.nc")
    dfJHU0 = DSJHU.to_dataframe().reset_index()
    dfJHU0.drop(columns = ['observed_value', 'model-data_mismatch', 'assimilation_flag','sounding_ID', 'latitude', 'longitude'],inplace = True)
    dfJHU0.rename(columns = {'model_estimate':'JHU_CO2','sounding':'sounding_id'},inplace = True)
    df = pd.merge(df,dfJHU0,on = ['sounding_id'],how = 'left')
    
    df.reset_index(drop=True,inplace=True)
    
    return df


# ---------------
# Extra functions
# ---------------

def InterpolateOnMaskXarray(DS, Numm,resLat, resLong, datasetname, intMethod, addTotalFluxes = False, Path_matrix = '', timeres = ''):
    #rectengular selcetion
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    XCO2 = DS.sel(latitude=slice(Lat_min-resLat, Lat_max+resLat), longitude=slice(Long_min-resLong, Long_max+resLong)) # take grid resolution as buffer
    
    if not os.path.isfile("/masks/netcdf/Mask"+str(Numm)+".nc"):
        CreateMask(Numm)
    mask = xr.open_dataset("/masks/netcdf/Mask"+str(Numm)+".nc")
    if intMethod == 'weightMean':
        DS_reg = RegridXarray(XCO2, datasetname, Long_min-0.5,Long_max+0.5, 1, Lat_min-0.5, Lat_max+0.5, 1, Path_matrix, 'weightMean')
    elif intMethod == 'non':
        DS_reg = XCO2.sel(latitude=slice(Lat_min, Lat_max), longitude=slice(Long_min, Long_max)) # take grid resolution as buffer
    else:
        DS_reg = XCO2.interp_like(mask.transcom_regions, method=intMethod)
    if Numm >= 900:
        DS_reg = DS_reg*mask.transcom_regions
    
    DS_reg = DS_reg.to_dataset(name = datasetname)

    if addTotalFluxes:
        if timeres == '':
            sys.exit('please provide timeres in InterpolateOnMaskXarray in CreateSingleDataFrames function')
        elif timeres == 'Month':
            NumDayOfMonth = [getNumDayOfMonth(int(x.replace('-','')[0:4]),int(x.replace('-','')[4:6])) for x in DS_reg.time.values]
            DoM = xr.DataArray(NumDayOfMonth, dims=['time'],coords = [DS_reg.time])
            DS_reg[datasetname+'_total'] = (DS_reg*mask.Area_Mask)[datasetname]*24*60*60*DoM
        elif timeres == 'Day':
            DS_reg[datasetname+'_total'] = (DS_reg*mask.Area_Mask)[datasetname]*24*60*60

    return DS_reg

