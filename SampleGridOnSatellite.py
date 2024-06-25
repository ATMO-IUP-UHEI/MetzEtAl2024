#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample Gridded 3D data on satellite
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import numpy as np
import pandas as pd
import xarray as xr
import read_remotec_out
from RegionParam import getRegion
import datetime
import glob
import sys


Numm = 756
year_min, month_min, day_min = 2009,4,15
year_max, month_max, day_max = 2018,12,31
timeperiod = [datetime.date(year_min,month_min,day_min),datetime.date(year_max, month_max,day_max)]
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)


#for each year
for year in range(year_min,year_max+1):
    print(year)
    print('load GOSAT')
    if year != 2018:
        GOSAT_path = glob.glob("Path to GOSAT data")[0]
        # Function to read GOSAT data
        GOSATout = read_remotec_out.read_output_file(GOSAT_path)

        #get datetime from year, month ...
        #No nicer way found yet
        timedf = pd.DataFrame({'year': GOSATout['GOSAT_TIME_YEAR'],
                        'month': GOSATout['GOSAT_TIME_MONTH'],
                        'day': GOSATout['GOSAT_TIME_DAY'],
                        'hour': GOSATout['GOSAT_TIME_HOUR'],
                        'minute': GOSATout['GOSAT_TIME_MIN'],
                        'second': GOSATout['GOSAT_TIME_SEC'] - 2})#GOSAT seconds are shifted by two seconds
        timearray = pd.to_datetime(timedf)

        GOSAT = xr.Dataset({
                "indexes": range(len(GOSATout['GOSAT_LATITUDE'])),
                "level": range(len(GOSATout["METEO_PRESSURE"][0])),
                "layer": range(len(GOSATout['REMOTEC-OUTPUT_X_STATE'].T[0:12].T[0])),
                "time": ("indexes",timearray),
                "latitude" : ("indexes",GOSATout['GOSAT_LATITUDE']),
                "longitude" : ("indexes",GOSATout['GOSAT_LONGITUDE']),
                "pressure" : (["indexes","level"],GOSATout["METEO_PRESSURE"]), #should contain an array/list of increasing pressure levels for each measurement
                "state_c" : (["indexes","layer"],GOSATout['REMOTEC-OUTPUT_X_STATE'].T[0:12].T),
                "column_c" : ("indexes",GOSATout['REMOTEC-OUTPUT_X_COLUMN_CORR'].T[0]),
                "CTcolumn_c" : ("indexes",GOSATout['REMOTEC-OUTPUT_X_APR_COLUMN'].T[0]),
                "state_apr" : (["indexes","layer"],GOSATout['REMOTEC-OUTPUT_X_APR_STATE'].T[0:12].T),
                "airmass" : (["indexes","layer"],GOSATout["METEO_AIRMASS"]),
                "AK" : (["indexes","layer"],GOSATout['REMOTEC-OUTPUT_X_AK_COLUMN'][:,0,:]),
                "mode" : ("indexes",GOSATout['(0=NADIR,1=GLINT)']),
                })
    else:
        GOSAT_path = " fp_short_fil_corr_210210_2018.out"
        GOSATout = read_remotec_out.read_output_file(GOSAT_path)
        GOSAT_path2 = " fp_short_fil_corr_201202_2018.out"
        GOSATout2 =read_remotec_out.read_output_file(GOSAT_path2)

        #get datetime from year, month ...
        #No nicer way found yet
        timedf = pd.DataFrame({'year': np.append(GOSATout['GOSAT_TIME_YEAR'],GOSATout2['GOSAT_TIME_YEAR']),
                        'month': np.append(GOSATout['GOSAT_TIME_MONTH'],GOSATout2['GOSAT_TIME_MONTH']),
                        'day': np.append(GOSATout['GOSAT_TIME_DAY'],GOSATout2['GOSAT_TIME_DAY']),
                        'hour': np.append(GOSATout['GOSAT_TIME_HOUR'],GOSATout2['GOSAT_TIME_HOUR']),
                        'minute': np.append(GOSATout['GOSAT_TIME_MIN'],GOSATout2['GOSAT_TIME_MIN']),
                        'second': np.append(GOSATout['GOSAT_TIME_SEC'] - 2,GOSATout2['GOSAT_TIME_SEC'] - 2)})#GOSAT seconds are shifted by two seconds
        timearray = pd.to_datetime(timedf)

        GOSAT = xr.Dataset({
                "indexes": range(len(np.append(GOSATout['GOSAT_LATITUDE'],GOSATout2['GOSAT_LATITUDE']))),
                "level": range(len(GOSATout["METEO_PRESSURE"][0])),
                "layer": range(len(GOSATout['REMOTEC-OUTPUT_X_STATE'].T[0:12].T[0])),
                "time": ("indexes",timearray),
                "latitude" : ("indexes",np.append(GOSATout['GOSAT_LATITUDE'],GOSATout2['GOSAT_LATITUDE'])),
                "longitude" : ("indexes",np.append(GOSATout['GOSAT_LONGITUDE'],GOSATout2['GOSAT_LONGITUDE'])),
                "pressure" : (["indexes","level"],np.append(GOSATout["METEO_PRESSURE"],GOSATout2["METEO_PRESSURE"],axis=0)), #should contain an array/list of increasing pressure levels for each measurement
                "state_c" : (["indexes","layer"],np.append(GOSATout['REMOTEC-OUTPUT_X_STATE'].T[0:12].T,GOSATout2['REMOTEC-OUTPUT_X_STATE'].T[0:12].T,axis=0)),
                "column_c" : ("indexes",np.append(GOSATout['REMOTEC-OUTPUT_X_COLUMN_CORR'].T[0],GOSATout2['REMOTEC-OUTPUT_X_COLUMN_CORR'].T[0])),
                "CTcolumn_c" : ("indexes",np.append(GOSATout['REMOTEC-OUTPUT_X_APR_COLUMN'].T[0],GOSATout2['REMOTEC-OUTPUT_X_APR_COLUMN'].T[0])),
                "state_apr" : (["indexes","layer"],np.append(GOSATout['REMOTEC-OUTPUT_X_APR_STATE'].T[0:12].T,GOSATout2['REMOTEC-OUTPUT_X_APR_STATE'].T[0:12].T,axis=0)),
                "airmass" : (["indexes","layer"],np.append(GOSATout["METEO_AIRMASS"],GOSATout2["METEO_AIRMASS"],axis=0)),
                "AK" : (["indexes","layer"],np.append(GOSATout['REMOTEC-OUTPUT_X_AK_COLUMN'][:,0,:],GOSATout2['REMOTEC-OUTPUT_X_AK_COLUMN'][:,0,:],axis=0)),
                "mode" : ("indexes",np.append(GOSATout['(0=NADIR,1=GLINT)'],GOSATout2['(0=NADIR,1=GLINT)'])),
                })

    print('load CT2022 datasets')
    if year > 2009:
        Model = xr.open_mfdataset([" CT2022/molefractions_co2_total_monthly/CT2022.molefrac_glb3x2_"+str(year-1)+"-12.nc"]
                +glob.glob(" CT2022/molefractions_co2_total_monthly/CT2022.molefrac_glb3x2_"+str(year)+"-*.nc")
                +[" CT2022/molefractions_co2_total_monthly/CT2022.molefrac_glb3x2_"+str(year+1)+"-01.nc"],combine='by_coords')
    else:
        Model = xr.open_mfdataset(glob.glob(" CT2022/molefractions_co2_total_monthly/CT2022.molefrac_glb3x2_"+str(year)+"-*.nc")
                +[" CT2022/molefractions_co2_total_monthly/CT2022.molefrac_glb3x2_"+str(year+1)+"-01.nc"],combine='by_coords')
    Model = Model.drop_vars(['v','u','time_components','temperature','specific_humidity','orography','gph','decimal_date','blh','air_mass','calendar_components'])
    dlat_model = abs(Model.latitude.values[1]-Model.latitude.values[0])
    dlong_model = abs(Model.longitude.values[1]-Model.longitude.values[0])
    Model = Model.where(((Model.latitude >= max(Lat_min -2*dlat_model,-90))
                        &(Model.latitude <= min(Lat_max +2*dlat_model,90))
                        &(Model.longitude >= max(Long_min -2*dlong_model,-180))
                        &(Model.longitude <= min(Long_max +2*dlong_model,180))),drop = True)
    Model = Model.where(((Model.time >= np.datetime64(str(year-1)+'-12-01'))&(Model.time <= np.datetime64(str(year+1)+'-01-31'))),drop = True)
    Model = Model.compute()

    GOSAT.set_coords(['latitude','longitude','time'])

    #GOSAT = GOSAT.where(((GOSAT.time >= np.datetime64('2018-05-01'))&(GOSAT.time <= np.datetime64('2018-05-31'))),drop = True)
    GOSAT = GOSAT.where(((GOSAT.latitude >= Lat_min)
                        &(GOSAT.latitude <= Lat_max)
                        &(GOSAT.longitude >= Long_min)
                        &(GOSAT.longitude <= Long_max)),drop = True)


    Modelint = Model.interp(latitude = GOSAT.latitude, longitude = GOSAT.longitude, time = GOSAT.time)
    #flip order of level and layer in model -> from TOA to surface like GOSAT
    Modelint = Modelint.isel(level=slice(None, None, -1))
    Modelint = Modelint.isel(boundary=slice(None, None, -1))
    
    ModelXCO2 = []
    ModelXCO2_wAK = []
    GOSATXCO2 = []
    GOSATXCO2_wAK = []
    #pressure interpolation by looping over all soundings and perform matrix multiplicatation
    for snd in range(len(GOSAT.time)):
        # get pressure -> this are borders!
        #GOSAT pressure = target
        GOSAT_p = GOSAT.pressure.values[snd,:]
        DFGOSAT_p = pd.DataFrame(data={'Pressure':GOSAT_p,'index_new':range(1,len(GOSAT_p)+1)})
        DFGOSAT_p.insert(1,column = 'dp_GOSAT', value = np.append(DFGOSAT_p.Pressure.values[1:13],[np.nan]) - DFGOSAT_p.Pressure)
        # Model pressure
        model_p = Modelint.pressure.values[snd,:]/100 # Pa -> hPa
        if model_p[-1] < GOSAT_p[-1]:
            model_p[-1] = GOSAT_p[-1]
        DFmodel_p = pd.DataFrame(data={'Pressure':model_p,'index_old':range(1,len(model_p)+1)})
            
        DFpressure = pd.merge(DFGOSAT_p,DFmodel_p,on='Pressure',how='outer')
        DFpressure.sort_values(by=['Pressure'],ignore_index=True,inplace=True)
        DFpressure.index_old.fillna(method="ffill", inplace = True)
        DFpressure.dp_GOSAT.fillna(method="ffill", inplace = True)
        DFpressure.index_new.fillna(method="ffill", inplace = True)
        DFpressure.index_new.fillna(0, inplace = True)
    
        
        DFpressure.insert(loc = 1,value=np.append(DFpressure.Pressure[1:],np.nan),column= 'Pressure_up')
        DFpressure.insert(loc = 1,value=(DFpressure.Pressure_up-DFpressure.Pressure)/DFpressure.dp_GOSAT,column= 'Pressure_weight')
        
        #set weight of model layers ouside the range of GOSAT to 0 -> not contributing
        DFpressure.loc[(DFpressure.Pressure < GOSAT_p[0])|(DFpressure.Pressure_up > GOSAT_p[-1]),'Pressure_weight'] = 0
        if len(DFpressure.index_new) != len(DFpressure.drop_duplicates(['index_new','index_old']).index_new):
            sys.exit('double entries in matrix')
        DFpressure.loc[(DFpressure.index_new == 0),'index_new'] = 1 #add the only 0 contributing models levels above TOA GOSAT (model pressure < min GOSAT pressure)), to the lowest GOSAT level -> avoids additional level
        DFpressure.loc[(DFpressure.index_new == 13),'index_new'] = 12 #add the only 0 contributing models levels (model pressure > max GOSAT pressure/SFP)), to the highest GOSAT level -> avoids additional level
        DFpressure.loc[(DFpressure.duplicated(['index_new','index_old'],keep = False)&(DFpressure.Pressure_weight == 0)),'Pressure_up'] = np.nan
        
        DFpressure.dropna(axis = 0,subset=['Pressure_up'],inplace = True) #drop last row
        
        DFpressure.drop(columns= ['Pressure','Pressure_up'],inplace=True)
        dfPressure_Matrix = pd.pivot(DFpressure, values='Pressure_weight', index='index_old', columns='index_new')
        dfPressure_Matrix.fillna(0,inplace=True)
        PressureMatrix = dfPressure_Matrix.to_numpy(copy=True)
        
        csModelCO2 = np.matmul(Modelint.co2.values[snd,:],PressureMatrix)
        csModelCO2airmass = csModelCO2 * GOSAT.airmass.values[snd,:]
        GOSATairmass = GOSAT.airmass.values[snd,:]
        state_apr_airmass = GOSAT.state_apr.values[snd,:]
        
        ModelXCO2.append(round(sum(csModelCO2airmass)/sum(GOSATairmass),6))
        GOSATXCO2.append(round(1000000*sum(GOSAT.state_c.values[snd,:])/sum(GOSATairmass),6))
        ModelXCO2_wAK.append(round(1000000*(np.dot(GOSAT.AK.values[snd,:],(csModelCO2airmass/1000000-state_apr_airmass)) 
                                    + sum(state_apr_airmass))/sum(GOSATairmass),6))
        GOSATXCO2_wAK.append(round(1000000*(np.dot(GOSAT.AK.values[snd,:],(GOSAT.state_c.values[snd,:]-state_apr_airmass)) 
                                    + sum(state_apr_airmass))/sum(GOSATairmass),6))
    
    d = {'CO2': ModelXCO2, 
         'CO2_AK':ModelXCO2_wAK,
         'Lat':Modelint.latitude.values,
         'Long':Modelint.longitude.values,
         'Year':pd.DatetimeIndex(GOSAT.time.values).year,
         'Month':pd.DatetimeIndex(GOSAT.time.values).month, 
         'Day':pd.DatetimeIndex(GOSAT.time.values).day,
         'Hour':pd.DatetimeIndex(GOSAT.time.values).hour, 
         'Min':pd.DatetimeIndex(GOSAT.time.values).minute,
         'Sec':pd.DatetimeIndex(GOSAT.time.values).second, 
         'RT_CO2':GOSAT.column_c.values,
         'RT_mode':GOSAT.mode.values}
    df = pd.DataFrame(data=d)
    
    if year == year_min:        
        #create Dataframe
        dfres= df.copy(deep = True)
    else:      
        #create Dataframe
        dfres = dfres.append(df, ignore_index=True)  
        
        
dfres.to_pickle(' DataFrames/DF0v0_CT2022cs_sounding_'+RegionName+str(Numm)+'.pkl')
        
        
    
    
    
    
    
    