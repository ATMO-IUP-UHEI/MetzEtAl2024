#!/usr/bin/env python
# functions for the Scripts 
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

import numpy as np
import pandas as pd
import datetime
import geopandas
from shapely.geometry import Polygon
from RegionParam import getRegion


def getMonthMeans(gdf,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,ValueName,Error_name = '',UnCorr_name='',weighting = 'Non'):
    # function to create monthly means of CO2 for different paper
    # input:
    # -- gdf: Geodataframe containing the data
    # -- year_max, month_max, day_max: enddate
    # -- year_min,month_min, day_min: startdate
    # -- Long_min, Long_max: Longitude range
    # -- Lat_min, Lat_max: Latitude range
    #
    # output:
    # -- gdf with monthly CO2 values
    if weighting == 'Area':
        RegGDF = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                            & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                            & (gdf.Long >= Long_min)
                            & (gdf.Long <= Long_max)
                            & (gdf.Lat >= Lat_min)
                            & (gdf.Lat <= Lat_max)]
        RegGDF.insert(1, column = ValueName+'xArea', value = RegGDF[ValueName]*RegGDF.Area)
        output_mean = pd.merge(RegGDF.groupby(['Year','Month'])[ValueName+'xArea'].sum().reset_index(),
                               RegGDF.groupby(['Year','Month'])['Area'].sum().reset_index(),
                               on = ['Year','Month'])
        output_mean.insert(1,column=ValueName,value =output_mean[ValueName+'xArea']/output_mean.Area)
    else:
        output_mean = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                            & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                            & (gdf.Long >= Long_min)
                            & (gdf.Long <= Long_max)
                            & (gdf.Lat >= Lat_min)
                            & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].mean().reset_index()
    output_count = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].count().reset_index()
    output_count.rename(columns = {ValueName:'number'}, inplace = True)
    output_all = pd.merge(output_mean, output_count, on=['Year','Month'])
        

    output_std = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].std(ddof = 0).reset_index()
    output_std.rename(columns = {ValueName:'StDev'}, inplace = True)
    output_all = pd.merge(output_all, output_std, on=['Year','Month'])
    output_all.insert(loc=1,column='StError',value= output_all.StDev/np.sqrt(output_all.number))

    if Error_name:
        output_error_all = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)]
        output_error_all.insert(loc=1,column='error_sq',value= output_error_all[Error_name]*output_error_all[Error_name])
        output_error = output_error_all.groupby(['Year','Month'])['error_sq'].sum().reset_index()
        output_all = pd.merge(output_all, output_error, on=['Year','Month']) 
        output_all.insert(loc=1,column=Error_name,value=np.sqrt(output_all.error_sq)/output_all.number)

    if UnCorr_name:
        output_raw_mean = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[UnCorr_name].mean().reset_index()
        output_all = pd.merge(output_all, output_raw_mean, on=['Year','Month'])
    
    output = output_all.copy(); print('Caution, min number of mean value = 0')
    #output = output_all[(output_all.number >= 10)]
    
    return output

def getMonthSum(gdf,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,ValueName,Errorname = ''):
    # function to create monthly sum of indicated column inside given borders
    # input:
    # -- gdf: Geodataframe containing the data
    # -- year_max, month_max, day_max: enddate
    # -- year_min,month_min, day_min: startdate
    # -- Long_min, Long_max: Longitude range
    # -- Lat_min, Lat_max: Latitude range
    # -- ValueName: Name of column to sum up
    #
    # output:
    # -- gdf with monthly summed values
    output_mean = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].sum().reset_index()

    output_count = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].count().reset_index()
    output_count.rename(columns = {ValueName:'number'}, inplace = True)
    output_all = pd.merge(output_mean, output_count, on=['Year','Month'])
    if len(Errorname) > 1:
        output_error_all = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)]
        output_error_all.insert(loc=1,column='error_sq',value= output_error_all[Errorname]*output_error_all[Errorname])
        output_error = output_error_all.groupby(['Year','Month'])['error_sq'].sum().reset_index()
        output_all = pd.merge(output_all, output_error, on=['Year','Month']) 
        output_all.insert(loc=1,column=Errorname,value=np.sqrt(output_all.error_sq))
        output_error2 = output_error_all.groupby(['Year','Month'])[Errorname].mean().reset_index()
        output_error2 = output_error2.rename(columns={Errorname:Errorname+'mean'})
        output_all = pd.merge(output_all, output_error2, on=['Year','Month']) 
    
    output = output_all[(output_all.number > 0)]



    return output

def getMonthSumFCday(gdf,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,ValueName):
    # function to create monthly sum of indicated column inside given borders
    # input:
    # -- gdf: Geodataframe containing the data
    # -- year_max, month_max, day_max: enddate
    # -- year_min,month_min, day_min: startdate
    # -- Long_min, Long_max: Longitude range
    # -- Lat_min, Lat_max: Latitude range
    # -- ValueName: Name of column to sum up
    #
    # output:
    # -- gdf with monthly summed values
    output_mean2 = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month','Day'])[ValueName].mean().reset_index()

    output_mean = output_mean2.groupby(['Year','Month'])[ValueName].sum().reset_index()


    output_count = gdf[(gdf.Date >= datetime.date(year_min,month_min,day_min))
                        & (gdf.Date <= datetime.date(year_max,month_max,day_max))
                        & (gdf.Long >= Long_min)
                        & (gdf.Long <= Long_max)
                        & (gdf.Lat >= Lat_min)
                        & (gdf.Lat <= Lat_max)].groupby(['Year','Month'])[ValueName].count().reset_index()
    output_count.rename(columns = {ValueName:'number'}, inplace = True)
    output_all = pd.merge(output_mean, output_count, on=['Year','Month'])
    output = output_all[(output_all.number > 0)]



    return output

def addMonthDate(DF):
    # function to add "MonthDate" column to Dtaframe based on column Year and Month
    MonthDate = []
    for m in range(len(DF['Year'])):
        try:
            MonthDate.append(datetime.date(int(DF.Year[m]),int(DF.Month[m]),15))
        except:
            print(DF.Year[m])
            print(DF.Month[m])
            MonthDate.append(datetime.date(int(DF.Year[m]),(DF.Month[m]),15))
    DF.insert(loc = 1,column='MonthDate',value = MonthDate)
    return DF

def getMeanAmplitude(df,varName, yearName, minYear, maxYear):
    '''
    # calculates the mean amplitude and its standard deviation of a time series 
    # with customizable year coordinate
    # input
    #   dataframe containing:
    #     value variable
    #     year coordinate
    #   name of value coordinate
    #   name of year coordinate
    #   Year to start the averaging 
    #   Year to end the averaging
    # returns:
    #   mean Amplitude, stDev, minimum and maximum as integers'''
    dfAmpli = df.groupby([yearName])[varName].min().reset_index()
    dfAmpli.rename(columns = {varName:'AMin'},inplace=True)
    dfAmpli2 = df.groupby([yearName])[varName].max().reset_index()
    dfAmpli2.rename(columns = {varName:'AMax'},inplace=True)
    dfAmpli = dfAmpli.merge(dfAmpli2,on=yearName)
    del dfAmpli2
    dfAmpli.insert(loc = 1, column = 'AAmpli',value = dfAmpli.AMax-dfAmpli.AMin)
    amplitude = dfAmpli[(dfAmpli[yearName] >= minYear)&(dfAmpli[yearName] <= maxYear)]['AAmpli'].mean()
    minimum = dfAmpli[(dfAmpli[yearName] >= minYear)&(dfAmpli[yearName] <= maxYear)]['AMin'].mean()
    maximum = dfAmpli[(dfAmpli[yearName] >= minYear)&(dfAmpli[yearName] <= maxYear)]['AMax'].mean()
    StDevAmplitude = dfAmpli[(dfAmpli[yearName] >= minYear)&(dfAmpli[yearName] <= maxYear)]['AAmpli'].std(ddof = 0)
    
    return amplitude, StDevAmplitude, minimum, maximum

def getReferenceDate(year_min,year_max,month_min,month_max):
    yea = []
    mon = []
    for k in range(year_min,year_max +1):
        if k == year_min:
            for p in range(month_min,12+1):
                yea.append(k)
                mon.append(p)
        elif k == year_max:
            for p in range(1,month_max+1):
                yea.append(k)
                mon.append(p)
        else:
            for p in range(1,13):
                yea.append(k)
                mon.append(p)

    dateData = {"Year":yea,"Month":mon}
    DateRef = pd.DataFrame(data=dateData)
    MonthDate2 = []
    for j in range(len(DateRef.Year)):
        MonthDate2.append(datetime.date(DateRef.Year[j],DateRef.Month[j],15))
    DateRef['MonthDate'] = MonthDate2

    #year from july to june index column
    jjYear = []
    for i in range(len(DateRef.Year)):
        if DateRef.Month.values[i] <= 6:
            jjYear.append(DateRef.Year.values[i]-1)
        elif DateRef.Month.values[i] >= 7:
            jjYear.append(DateRef.Year.values[i])
        else:
            print('Error')
    DateRef.insert(loc = 1, column= 'JJYear',value = jjYear)

    return DateRef

def getNumDayOfMonth(year,month):
    listDays = [31,28,31,30,31,30,31,31,30,31,30,31]
    listDaysl = [31,29,31,30,31,30,31,31,30,31,30,31]
    if year < 1900:
        print('year is out of implemented range, check code')
    elif year in list(range(1904,2100,4)):
        days = listDaysl[month-1]
    else:
        days = listDays[month-1]

    return days

def getReferencesDateDay(year_min,year_max,month_min,month_max,day_min,day_max):
    
    yea = []
    mon = []
    day = []
    for k in range(year_min,year_max +1):
        if k == year_min:
            for p in range(month_min,12+1):
                if p == month_min:
                    for d in range(day_min, getNumDayOfMonth(k,p)+1):
                        day.append(d)
                        yea.append(k)
                        mon.append(p)
                else:
                    for d in range(1, getNumDayOfMonth(k,p)+1):
                        day.append(d)
                        yea.append(k)
                        mon.append(p)
        elif k == year_max:
            for p in range(1,month_max+1):
                if p == month_max:
                    for d in range(1, day_max):
                        day.append(d)
                        yea.append(k)
                        mon.append(p)
                else:
                    for d in range(1, getNumDayOfMonth(k,p)+1):
                        day.append(d)
                        yea.append(k)
                        mon.append(p)
        else:
            for p in range(1,13):
                for d in range(1, getNumDayOfMonth(k,p)+1):
                    day.append(d)
                    yea.append(k)
                    mon.append(p)

    dateData = {"Year":yea,"Month":mon,"Day":day}
    DateRef = pd.DataFrame(data=dateData)
    MonthDate2 = []
    for j in range(len(DateRef.Year)):
        MonthDate2.append(datetime.date(DateRef.Year[j],DateRef.Month[j],DateRef.Day[j]))
    DateRef['MonthDate'] = MonthDate2

    return DateRef

def getAreaOfGrid():
    """Get Dataframe with Area dependend on Latitude for a 1°x1° grid"""
    AreaLat = []
    Lat = range(895,-905,-10)
    for i in Lat:
        geom = [Polygon(zip([100,100,101,101],[i/10-0.5,i/10+0.5,i/10+0.5,i/10-0.5]))]
        GData = geopandas.GeoDataFrame({'data':[1]}, geometry=geom)
        GData.crs = 'epsg:4326'
        GData = GData.to_crs({'proj':'cea'})
        AreaLat.append(GData.geometry.area[0])
    dfArea = pd.DataFrame({'Area':AreaLat,'Lat':np.array(Lat)/10})
    
    return(dfArea)

def getGdfOfGrid(Num):
    """Get Geodatframe of 1°x1° Grid for a certain transcom region"""
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Num)
    if Num >= 900:
        Transcom = pd.read_pickle("/Transcom_Regions.pkl")
        first = True
        AreaLat = []
        Lat = range(Lat_max*10+5, Lat_min*10-5,-10)
        if Long_max >= 179.5:
            Long_max = Long_max -1
        Long = range(Long_max*10+5, Long_min*10-5,-10)
        for i in Lat:
            for j in Long:
                geom = [Polygon(zip([j/10-0.5,j/10-0.5,j/10+0.5,j/10+0.5],[i/10-0.5,i/10+0.5,i/10+0.5,i/10-0.5]))]
                GData = geopandas.GeoDataFrame({'data':[1]}, geometry=geom)
                GData.crs = 'epsg:4326'
                GData.insert(loc=1,column = 'Lat',value = [i/10])
                GData.insert(loc=1,column = 'Long',value = [j/10])
                if first:
                    gdf = GData.copy()
                    first = False
                else:
                    gdf = gdf.append(GData, ignore_index=True)
                GData = GData.to_crs({'proj':'cea'})
                AreaLat.append(GData.geometry.area[0])
        
        gdf.insert(loc=1, value = AreaLat,column = 'AreaGrid') 
        gdf = gdf.drop(columns = 'data')
        gdf.insert(loc=1, value = gdf.geometry,column = 'geomPoly') 
        gdf.insert(loc=1, column = 'geomPoint', value = geopandas.points_from_xy(gdf.Long, gdf.Lat))
        gdf = gdf.set_geometry(gdf.geomPoint)
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
        
        gdf = gdf.set_geometry(gdf.geomPoly)
        
    elif Num >= 700:
        Transcom = pd.read_pickle("/masks/CTTranscomMask1x1Borders.pkl")
        gdf = pd.read_pickle("/masks/CTTranscomMask1x1Grid.pkl")
        
        gdf = gdf.set_geometry(gdf.geomPoint)
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
        gdf = gdf.set_geometry(gdf.geomPoly)
    else:
        raise Exception('Not yet implemented for Num < 700')
    gdf = gdf.reset_index()
    gdf.insert(loc=1,value=gdf.index,column='GridID')
    gdf.to_pickle("/Grid1_1"+RegionName+".pkl")
    gdf = gdf[(gdf.Long >= Long_min)
         & (gdf.Long <= Long_max)
         & (gdf.Lat >= Lat_min)
         & (gdf.Lat <= Lat_max)]
    gdf.to_pickle("/Grid1_1"+str(Num)+".pkl")
    
    return gdf

def getEnsembleFluxMIP(PMIPOCOModel, PMIPpriorModel, MIPModelNames, MIPModelNamesPrior,MIPEnsNames):
    '''
    # calculates the mean carbon flux of posteriori and prior  MIP fluxes 
    # input
    #   PMIPOCOModel: list of MIP Model dataframes containing:
    #     Landtot variable POSTERIORI fluxes
    #     monthDate, Year, month variable
    #   PMIPpriorModel: list of MIP Model dataframes containing:
    #     Landtot variable Prior Fluxes
    #     monthDate, Year, month variable
    #   MIPModelNames: List of MIP model names matching order in dataframe list
    #   MIPModelNamesPriorname: List of MIP model names matching order in prior dataframe list
    #   MIPEnsNames: list of MIP model names which should go into the ensemble
    # returns:
    #   EnsembleFluxMIP: dataframe with mean of posteriori and prior fluxes of models in MIPEnsNames and post fluxes values for each moel
    '''
    EnsembleFluxMIP = PMIPOCOModel[0][['MonthDate','Month','Year']].copy()

    MIPEnsNamesCorr = []
    for Nam0 in MIPEnsNames:
        if Nam0 in MIPModelNamesPrior: #Lofi is not
            EnsembleFluxMIP = pd.merge(EnsembleFluxMIP,
                                        PMIPpriorModel[np.where(np.array(MIPModelNamesPrior) == Nam0)[0][0]]
                                            .rename(columns = {'Landtot':Nam0})[[Nam0,'MonthDate','Month','Year']],
                                        on = ['MonthDate','Month','Year'])
            MIPEnsNamesCorr.append(Nam0)
        else:
            print(Nam0 + ' has no prior, skipt in Ensemble Prior calculation')
    EnsembleFluxMIP.insert(0,column = 'meanPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].mean(axis=1))
    EnsembleFluxMIP.insert(0,column = 'StdPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].std(ddof = 0,axis=1))
    EnsembleFluxMIP.insert(0,column = 'minPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].min(axis=1))
    EnsembleFluxMIP.insert(0,column = 'maxPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].max(axis=1))
    EnsembleFluxMIP.insert(0,column = 'low50nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.25,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up50nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.75,axis=1))
    EnsembleFluxMIP.insert(0,column = 'low75nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.875,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up75nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.125,axis=1))
    EnsembleFluxMIP.insert(0,column = 'low80nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.9,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up80nbpPrior', value = EnsembleFluxMIP[MIPEnsNamesCorr].quantile(0.1,axis=1))
    EnsembleFluxMIP.drop(columns = MIPEnsNamesCorr, inplace = True)

    for Nam0 in MIPEnsNames:
        EnsembleFluxMIP = pd.merge(EnsembleFluxMIP,
                                    PMIPOCOModel[np.where(np.array(MIPModelNames) == Nam0)[0][0]]
                                        .rename(columns = {'Landtot':Nam0})[[Nam0,'MonthDate','Month','Year']],
                                    on = ['MonthDate','Month','Year'])
    EnsembleFluxMIP.insert(0,column = 'meannbp', value = EnsembleFluxMIP[MIPEnsNames].mean(axis=1))
    EnsembleFluxMIP.insert(0,column = 'Stdnbp', value = EnsembleFluxMIP[MIPEnsNames].std(ddof = 0,axis=1))
    EnsembleFluxMIP.insert(0,column = 'minnbp', value = EnsembleFluxMIP[MIPEnsNames].min(axis=1))
    EnsembleFluxMIP.insert(0,column = 'maxnbp', value = EnsembleFluxMIP[MIPEnsNames].max(axis=1))
    EnsembleFluxMIP.insert(0,column = 'low50nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.25,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up50nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.75,axis=1))
    EnsembleFluxMIP.insert(0,column = 'low75nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.875,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up75nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.125,axis=1))
    EnsembleFluxMIP.insert(0,column = 'low80nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.9,axis=1))
    EnsembleFluxMIP.insert(0,column = 'up80nbp', value = EnsembleFluxMIP[MIPEnsNames].quantile(0.1,axis=1))
    #EnsembleFluxMIP.drop(columns = MIPEnsNames, inplace = True)
    
    return EnsembleFluxMIP

def getMSC(df,datavar,Year_min=2009,Year_max=2023):
    '''
    # Function to get mean seasonal cycle of a monthly Dataset
    # arguments:
    #           df: dataframe containing at least Month and datavar column
    #           datavar: name of column containing data
    # returns:
    #           df with Mean seasonal cycle
    '''
    if 'Month' not in df.keys():
        df.insert(1, column = 'Month',value = df.apply(lambda x: x.MonthDate.month,axis = 1))
    result = df[(df.Year >= Year_min)&(df.Year <= Year_max)].groupby(['Month'])[datavar].mean().reset_index()

    return result
    
def GoodModelsNames(Datatype, Num):
    # define the subselections of TRENDY , dependend on region
    if Num in [756, 766, 769, 770, 771, 772, 784]:
        return ['ORCHIDEE'+Datatype,
                'ORCHIDEEv3'+Datatype,
                'CABLE-POP'+Datatype] 
    else:
        print('TRENDY subselections not yet defined for region '+str(Num))           
def EarlyGPPModelNames(Datatype, Num):
    if Num in [756, 766, 769, 770, 771, 772, 784]:
        #before 19/2 'JSBACH','LPX-Bern','YIBs','SDGVM','CLASSIC','DLEM'
        return ['ISBA-CTRIP'+Datatype,
                'OCN'+Datatype] 
    else:
        print('TRENDY subselections not yet defined for region '+str(Num)) 
def LateRespNames(Datatype, Num): 
    if Num in [756, 766, 769, 770, 771, 772, 784]:
        return ['ORCHIDEE-CNP'+Datatype,
             'CLM5.0'+Datatype,#
             'VISIT'+Datatype,
             'ISAM'+Datatype,
             'IBIS'+Datatype,
             'JULES-ES-1p0'+Datatype,
             'LPJ'+Datatype,
             'JSBACH'+Datatype,
             'LPX-Bern'+Datatype,
             'YIBs'+Datatype,
             'SDGVM'+Datatype,
             'CLASSIC'+Datatype,
             'DLEM'+Datatype,] 
    else:
        print('TRENDY subselections not yet defined for region '+str(Num)) 
def BadModelNames(Datatype, Num):
    return EarlyGPPModelNames(Datatype, Num) + LateRespNames(Datatype, Num)

def NormalizeMeanCycle(MeanCycle, Var='', returnAmplFactor = False):
    '''
    # function to normalize mean seasonal cycle 
    # arguments:
    #           MeanCycle: Pandas Dataframe with variable column or numpy array
    #           Var: variable column name or '' for numpy array
    # returns:
    #           NMeanCycle: normalized mean seasonal cycle values
    '''
    if Var == '':
        NMeanCycle = ((MeanCycle-MeanCycle.min()-0.5*(MeanCycle.max()-MeanCycle.min()))
                      /abs((MeanCycle-MeanCycle.min()-0.5*(MeanCycle.max()-MeanCycle.min()))).max())
    else:
        NMeanCycle = ((MeanCycle[Var]-MeanCycle[Var].min()-0.5*(MeanCycle[Var].max()-MeanCycle[Var].min()))
                      /abs((MeanCycle[Var]-MeanCycle[Var].min()-0.5*(MeanCycle[Var].max()-MeanCycle[Var].min()))).max())
    AmplFactor = abs((MeanCycle[Var]-MeanCycle[Var].min()-0.5*(MeanCycle[Var].max()-MeanCycle[Var].min()))).max()
    
    if returnAmplFactor:
        return NMeanCycle, AmplFactor
    else:
        return NMeanCycle
