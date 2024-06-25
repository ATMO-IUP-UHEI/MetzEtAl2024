#!/usr/bin/env python
# -*- coding: utf-8 -*-
# functions #2 
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz


import numpy as np
import pandas as pd
import sys
import datetime
import xarray as xr
from RegionParam import getRegion


def getTimeMeans(gdf,ValueName,timeInterval,minNum,Error_name = '',UnCorr_name=''):
    '''
    # function to create means over a specified time period
    # arguments:
    #           gdf: Geodataframe containing the data and the columns Year and month
    #           ValueName: Name of column to be aggregated temporally
    #           timeInterval: 'Month' or 'Day' or 'Year'
    #           minNum: threshold: min number of measurements per month
    #           Error_name: name of the column containing the value uncertainties
    #           Uncorr_name: name of the column containing the values uncorrected data
    # output:
    #           gdf with monthly CO2 values, values in month with to less measurements are set nan
    '''
    # get timecomponents to aggregate on
    timeList = ['Year', 'Month', 'Day']
    timeList = timeList[:timeList.index(timeInterval)+1]
    if timeList[-1] not in gdf.keys():
        sys.exit("GDF does not contain column '"+timeList[-1]+"' necessary to calculate "+timeList[-1]+"ly means")

    # monthly mean
    output_mean = gdf.groupby(timeList)[ValueName].mean().reset_index()
    output_mean.rename(columns = {ValueName:ValueName+'_'+str(minNum)}, inplace = True)

    # number of measurements per month
    output_count = gdf.groupby(timeList)[ValueName].count().reset_index()
    output_count.rename(columns = {ValueName:'number'}, inplace = True)
    output_all = pd.merge(output_mean, output_count, on=timeList)

    #standard deviation within month
    output_std = gdf.groupby(timeList)[ValueName].std(ddof=0).reset_index()
    output_std.rename(columns = {ValueName:'StDev'+'_'+str(minNum)}, inplace = True)
    output_all = pd.merge(output_all, output_std, on=timeList)

    # standard error
    output_all.insert(loc=1,column='StError'+'_'+str(minNum),value= output_all['StDev'+'_'+str(minNum)]/np.sqrt(output_all.number))
    columns = [ValueName+'_'+str(minNum),'StDev'+'_'+str(minNum),'StError'+'_'+str(minNum)]

    #error propagation
    if UnCorr_name:
        output_raw_mean = gdf.groupby(timeList)[UnCorr_name].mean().reset_index()
        output_raw_mean.rename(columns = {UnCorr_name:UnCorr_name+'_'+str(minNum)}, inplace = True)
        output_all = pd.merge(output_all, output_raw_mean, on=timeList)
        columns.append(UnCorr_name+'_'+str(minNum))
    if Error_name:
        gdf.insert(loc=1,column='error_sq'+'_'+str(minNum),value= gdf[Error_name]*gdf[Error_name])
        output_error = gdf.groupby(timeList)['error_sq'+'_'+str(minNum)].sum().reset_index()
        gdf.drop(columns = 'error_sq'+'_'+str(minNum),inplace = True)
        output_all = pd.merge(output_all, output_error, on=timeList) 
        #error propagation error = sqrt(sum(error^2))/number
        output_all.insert(loc=1,column=Error_name+'_'+str(minNum),value=np.sqrt(output_all['error_sq'+'_'+str(minNum)])/output_all.number)
        del output_all['error_sq'+'_'+str(minNum)]
        columns.append(Error_name+'_'+str(minNum))

    #set all values in month with less measurements than threshold to nan
    output_all.loc[(output_all.number < minNum),columns] = np.nan

    return output_all

def getTimeSums(gdf,ValueName,timeInterval):
    '''
    # function to create sums over a specified time period
    # arguments:
    #           gdf: Geodataframe containing the data and the columns Year and month
    #           ValueName: Name of column to be aggregated temporally
    #           timeInterval: 'Month' or 'Day' or 'Year'
    # output:
    #           gdf with monthly flux sum
    '''
    # get timecomponents to aggregate on
    timeList = ['Year', 'Month', 'Day']
    timeList = timeList[:timeList.index(timeInterval)+1]
    if timeList[-1] not in gdf.keys():
        sys.exit("GDF does not contain column '"+timeList[-1]+"' necessary to calculate "+timeList[-1]+"ly means")

    # monthly mean
    print(ValueName)
    output = gdf.groupby(timeList)[ValueName].sum().reset_index()
    
    return output


def DetrendMonthdate(Month_means_CO2, offset, minNum, ValueName, rate = ''):
    ''' 
    # function to detrend the monthly data
    # arguments:
    #           Month_means_CO2: DataFrame containing a "Date" colume with a datetime date and a value column
    #           offsets: offset for April 2009
    #           minNum: min number of values per month
    #           ValueName: name of the value column
    #           rate: optional individual rate
    # returns:
    #           Month_means_CO2: Dataframe with the Detrended timeseries          
    '''
    if offset == 0:
        SetMeanZero = True # calculate detrended timeseries with mean = 0
    else:
        SetMeanZero = False 
    if not rate or 'Noaa' in rate: 
        # create background dataset
        rate = [1.58,2.42,1.68,2.41,2.45,2.04,2.95,2.83,2.14,2.39,2.52,2.36,2.45,2.17,2.17]#2023 rate = 2022 as it is not published yet! updated on 08.01.2024
        if 'CH4' in ValueName:
            rate = [0.0047,0.00519,0.00483,0.00501,0.0057,0.01277,0.01002,0.00709,0.00685,0.00867,0.00989,0.01527,0.01699]#updated on 08.04.2022
            #rate = [0.00467,0.00521,0.00485,0.00502,0.00569,0.01276,0.01003,0.00708,0.00687,0.00854,0.01034]#old before 08.04.2022
        ratename = 'rateNoaa'
    elif rate == 'NIESprior': #rate NOAA - 4ppm/13yrs
        rate = np.array([1.58,2.41,1.69,2.41,2.43,2.05,2.96,2.83,2.16,2.30,2.57,2.34,2.66])
        rate = list(rate - np.ones(13)*3/(13)) #ppm/year
        ratename = 'rateNIESprior'
    else:
        ratename = 'rate'+str(rate)
    detrend_all_month = []
    detrend_year = []
    detrend_month = []
    
    for y in range(2009,2024):
        if y == 2009:
            for m in range(4,13):
                detrend_all_month.append(offset+rate[y-2009]/12*(m-4))
                detrend_month.append(m)
                detrend_year.append(y)
        else:
            for m in range(1,13):
                detrend_all_month.append(detrend_all_month[-1]+rate[y-2009]/12)
                detrend_month.append(m)
                detrend_year.append(y)
    data_month = {'Background' : detrend_all_month,'Year':detrend_year,'Month': detrend_month}
    BG = pd.DataFrame(data = data_month)
    
    Month_means_CO2 = pd.merge(Month_means_CO2, BG, on = ['Year','Month']).reset_index()

    Month_means_CO2.drop(columns = ['index'],inplace = True)
    DetrendColumnName = ValueName+'_Detrend'+str(offset)+'_'+str(minNum)+'_'+ratename
    Month_means_CO2.insert(loc=1,column =DetrendColumnName,value = Month_means_CO2[ValueName+'_'+str(minNum)] - Month_means_CO2['Background'])
    MonthDate = Month_means_CO2.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    Month_means_CO2.insert(loc = 1, column = 'MonthDate', value = MonthDate)
    if SetMeanZero:
        meanMeans = Month_means_CO2[DetrendColumnName].mean()
        Month_means_CO2.insert(loc=1,column=ValueName+'_DetrendMean_'+str(minNum)+'_'+ratename,value=(Month_means_CO2[DetrendColumnName] - meanMeans))         
    return Month_means_CO2


def getReferenceDate(year_min,year_max,month_min,month_max):
    '''
    # function returning a pandas dataframe with year, month and datetime.date (day is set to 15th) column for the given time periode
    # arguments:
    #           year_min: start year of the time period
    #           year_max: end year of the time period
    #           month_min: first month in the start year of the time period  
    #           month_max: last month in the end year of the time period
    # returns:
    #           df: containing Year, Month and MonthDate column for all month in the given time period
    '''
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
    MonthDate2 = DateRef.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    DateRef.insert(loc=1,column='MonthDate',value=MonthDate2)

    return DateRef

def getReferenceDateDay(year_min,year_max,month_min,month_max):
    yea = []
    mon = []
    day = []
    numday = []
    for k in range(year_min,year_max +1):
        len0 = len(yea)
        if k == year_min:
            for p in range(month_min,12+1):
                for d in range(1,getNumDayOfMonth(k,p)+1):
                    yea.append(k)
                    mon.append(p)
                    day.append(d)
        elif k == year_max:
            for p in range(1,month_max+1):
                for d in range(1,getNumDayOfMonth(k,p)+1):
                    yea.append(k)
                    mon.append(p)
                    day.append(d)
        else:
            for p in range(1,13):
                for d in range(1,getNumDayOfMonth(k,p)+1):
                    yea.append(k)
                    mon.append(p)
                    day.append(d)
        numday = numday + list(range(len(yea)-len0))
        

    dateData = {"Year":yea,"Month":mon,"Day":day,"NumDayInYear":numday}
    DateRef = pd.DataFrame(data=dateData)
    MonthDate2 = []
    for j in range(len(DateRef.Year)):
        MonthDate2.append(datetime.date(DateRef.Year[j],DateRef.Month[j],DateRef.Day[j]))
    DateRef['Date'] = MonthDate2

    return DateRef

def CreateMask(Numm):
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    
    #Transcom name to number
    NameNr = pd.DataFrame(data = {'Name':['AU','SAT','SA','SATr','NA','Af'],'Num':[[10],[4],[6],[3],[5],[5,6]]})
    TranscomNr = NameNr[(NameNr.Name == RegionName)].Num.values[0]
    
    DS = xr.open_dataset("/regions.nc")
    DSt = DS.transcom_regions
    DS_clip = DSt.sel(latitude=slice(Lat_min, Lat_max), longitude=slice(Long_min, Long_max))
    if Numm >= 900:
        if len(TranscomNr) == 1:
            DS_clip = DS_clip.where(DS_clip == TranscomNr[0],np.nan)
            DS_clip = DS_clip.where(DS_clip != TranscomNr[0],1)
        elif len(TranscomNr) == 2:
            DS_clip = DS_clip.where((DS_clip == TranscomNr[0]) | (DS_clip == TranscomNr[1]),np.nan)
            DS_clip = DS_clip.where((DS_clip != TranscomNr[0]) & (DS_clip != TranscomNr[1]),1)
    else:
        DS_clip = DS_clip.where(DS_clip < 0,1)
    
    #get areas
    dfArea = getAreaOfGrid()
    dfArea2 = dfArea[(dfArea.Lat >= Lat_min)&(dfArea.Lat <= Lat_max)].sort_values('Lat')
    AreasA = np.tile(dfArea2.Area,(len(DS_clip.longitude),1))
    DSArea = xr.DataArray(AreasA.transpose(),dims={'latitude':dfArea2.Lat,'longitude': DS_clip.longitude.values})

    DDS_clip = DS_clip.to_dataset()
    DDS_clip['Area'] = DSArea
    DDS_clip['Area_Mask'] = DDS_clip.Area.where(DDS_clip.transcom_regions == 1,np.nan)
    DDS_clip.to_netcdf("/masks/netcdf/Mask"+str(Numm)+".nc")

def get_date_Components(df, timecolumn, insert=bool, dateformat=str):
    '''
    # Function to get the time components of a datetime64 column in a dataframe. Returns the components as Year, Month, Day, Hour, Minute, Seconds (insert = False) or inserts the components into the dataframe.
    # arguments:
    #           df: dataframe containing a column with datetime64 values
    #           timecolumn: name of the column containing the datetime64 dates
    #           insert: boolean, whether to insert the date components and datetime date into the dataframe (=True) or whether to return them
    #           dateformat: choose from 'datetime64','fracDays', 'fracDays2000'
    '''
    if dateformat == 'fracDays2000':
        datetime64time = df.apply(lambda x: datetime.datetime(2000,1,1,0,0,0)+datetime.timedelta(days=x[timecolumn]),axis=1)
        df.drop(columns = timecolumn,inplace=True)
        df.insert(loc=1, column =timecolumn, value = datetime64time)
    elif dateformat == 'fracDays':
        datetime64time = df.apply(lambda x: datetime.datetime(1970,1,1,0,0,0)+datetime.timedelta(days=x[timecolumn]),axis=1)
        df.drop(columns = timecolumn,inplace=True)
        df.insert(loc=1, column =timecolumn, value = datetime64time)
    elif dateformat == 'biteStr':
        s = df.apply(lambda x: int(round(float(str(x.time.decode("utf-8"))[17:]))),axis=1).values
        Iwrongsec = np.where(s >= 60)
        s[Iwrongsec] = 0
        df.insert(1,column = 'secround',value = s)
        datetime64time = df.apply(lambda x: datetime.datetime(int(str(x.time.decode("utf-8"))[0:4]),
                                                              int(str(x.time.decode("utf-8"))[5:7]),
                                                              int(str(x.time.decode("utf-8"))[8:10]),
                                                              int(str(x.time.decode("utf-8"))[11:13]),
                                                              int(str(x.time.decode("utf-8"))[14:16]),
                                                              x.secround),axis=1)
        datetime64time.iloc[Iwrongsec] = datetime64time.iloc[Iwrongsec].apply(lambda x: x+datetime.timedelta(0,0,0,0,1,0)) #add one minute to date former containing sec=60
        df.drop(columns = [timecolumn,'secround'],inplace=True)
        df.insert(loc=1, column =timecolumn, value = datetime64time)
    if dateformat == 'datetime64' or 'fracDays' in dateformat or dateformat == 'biteStr':
        Y = df.apply(lambda x: int(str(x[timecolumn])[0:4]),axis=1)
        M = df.apply(lambda x: int(str(x[timecolumn])[5:7]),axis=1)
        D = df.apply(lambda x: int(str(x[timecolumn])[8:10]),axis=1)
        h = df.apply(lambda x: int(str(x[timecolumn])[11:13]),axis=1)
        m = df.apply(lambda x: int(str(x[timecolumn])[14:16]),axis=1)
        s = df.apply(lambda x: int(str(x[timecolumn])[17:19]),axis=1)
        
    if insert: 
        df.insert(loc=1,column = 'Year',value = Y)
        df.insert(loc=1,column = 'Month',value = M)
        df.insert(loc=1,column = 'Day',value = D)
        df.insert(loc=1,column = 'Hour',value = h)
        df.insert(loc=1,column = 'Min',value = m)
        df.insert(loc=1,column = 'Sec',value = s)
        return df
    else:
        return Y,M,D,h,m,s

def SelectTimeAndRegion(df, year_min, month_min, day_min, year_max, month_max, day_max, Long_min, Long_max, Lat_min, Lat_max):
    '''
    # Function to select rows in column within a time period and rectengular region
    # arguments:
    #           df: dataframe containing Year, Month, Day, Lat and Long column
    #           year_min, month_min, day_min: begin of time period (included)
    #           year_max, month_max, day_max: max of time period (included)
    #           Long_min, Long_max, Lat_min, Lat_max: region boundaries
    # returns:
    #           df with temporal and spatial subselection of rows
    '''
    df = df[(df.Year >= year_min)&(df.Year >= year_max)&(df.Month >= month_min)&(df.Month <= month_max)&(df.Day >= day_min)&(df.Day <= day_max)
                    & (df.Long >= Long_min)
                    & (df.Long <= Long_max)
                    & (df.Lat >= Lat_min)
                    & (df.Lat <= Lat_max)]
    return df    

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
    if 'Year' not in df.keys():
        df.insert(1, column = 'Year',value = df.apply(lambda x: x.MonthDate.year,axis = 1))
    result = df[(df.Year >= Year_min)&(df.Year <= Year_max)].groupby(['Month'])[datavar].mean().reset_index()

    return result

def getNumDayOfMonth(year,month):
    """returns list of number of days within given month in given year"""
    listDays = [31,28,31,30,31,30,31,31,30,31,30,31]
    listDaysl = [31,29,31,30,31,30,31,31,30,31,30,31]
    if year < 1900:
        print('year is out of implemented range, check code')
    elif year in list(range(1904,2100,4)):
        days = listDaysl[month-1]
    else:
        days = listDays[month-1]

    return days

def CombineTimeseries(dfs: list,dfnames: list, dataVar, Year_min = 2009, Year_max = 2023):
    '''
    # Function to combine multiple timeseries dataframes and calculate mean, min, max and std of given variable dataVar
    # arguments:
    #           dfs: list of datframes to be combined, need column 'MonthDate'
    #           dfnames: list of name of dataframes dataset
    #           dataVar: name of variable to calculate min, max, mean, std
    # returns:
    #           dfComb: combined dataframe with MonthDate, Month and mean, min, max and std of given variable
    #           dfCombYmean: Mean seasonal cycle mean
    #           dfCombYstd: Std. over time of Mean seasonal cycle
    '''
    for i, df in enumerate(dfs):
        if i == 0:
            dfComb = df[['MonthDate',dataVar]].copy()
            MonthlistM = dfComb.apply(lambda x: x.MonthDate.month,axis=1)
            YearlistY = dfComb.apply(lambda x: x.MonthDate.year,axis=1)
            dfComb.insert(loc = 1, column = 'Month', value = MonthlistM)
            dfComb.insert(loc = 1, column = 'Year', value = YearlistY)
            
        else:
            dfComb = pd.merge(dfComb,df[['MonthDate',dataVar]],on='MonthDate',how='inner')
        dfComb = dfComb.rename(columns = {dataVar:dfnames[i]})
   
    dfComb.insert(loc = 1, column = 'meanValue', value = dfComb[dfnames].mean(axis = 1,skipna = False))
    dfComb.insert(loc = 1, column = 'stdValue', value = dfComb[dfnames].std(axis = 1,ddof = 0,skipna = False))
    dfComb.insert(loc = 1, column = 'minValue', value = dfComb[dfnames].min(axis = 1,skipna = False))
    dfComb.insert(loc = 1, column = 'maxValue', value = dfComb[dfnames].max(axis = 1,skipna = False))
    #skipna = False: if one dataset has nan, result also nan
    dfCombYmean = dfComb[(dfComb.Year >= Year_min)&(dfComb.Year <= Year_max)].groupby('Month')['meanValue'].mean().reset_index()
    dfCombYstd = dfComb[(dfComb.Year >= Year_min)&(dfComb.Year <= Year_max)].groupby('Month')['meanValue'].std(ddof = 0).reset_index()
    #calculate Mean Seasonal Cycle of individual Datasets and get min/max 
    dfCombYDataSetMean = dfComb[(dfComb.Year >= Year_min)&(dfComb.Year <= Year_max)].groupby('Month')[dfnames].mean().reset_index()
    dfCombYDataSetMean.insert(1, column = 'minValue', value = dfCombYDataSetMean[dfnames].min(axis = 1,skipna = False))
    dfCombYDataSetMean.insert(1, column = 'maxValue', value = dfCombYDataSetMean[dfnames].max(axis = 1,skipna = False))
    
    return dfComb, dfCombYmean, dfCombYstd, dfCombYDataSetMean

def addMonthDate(df):
    '''
    # Function to to insert MonthDate column with datetime(year,month,15), day = 15
    # arguments:
    #           df: dataframe, need column 'Month' and 'Year'
    # returns:
    #           df with additional MotnhDateColumn
    '''
    MonthDate = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    df.insert(loc = 1, column = 'MonthDate', value = MonthDate)
    return df

def getDataInfo(Dataset):
    '''
    # function to get informations about dataset stored in Data_Overview.csv
    # arguments:
    #           Dataset: 'DataName_resolution' as given in /DataFrames/Data_Overview.csv
    # returns:
    #           DataName, DataRes, datafiles, DataDir, DFversion, Datainfo as defined in csv file
    '''
    Datainfo = pd.read_csv("/DataFrames/Data_Overview.csv", index_col=[0])
    DataName, DataRes = Dataset.split('_')
    Datainfo_sel = Datainfo[(Datainfo.DataName == DataName)&(Datainfo.DataRes ==DataRes)].copy()
    if len(Datainfo_sel) == 0:
        print('Please enter the new dataset "'+Dataset+'" in the list: /DataFrames/Data_Overview.csv')
        sys.exit()
    DFversion = int(float(Datainfo_sel.DFversion.values[0]))
    DataDir = Datainfo_sel.DataDir.values[0]
    datafiles = Datainfo_sel.datafiles.values[0]
    print('get '+ Dataset )

    return DataName, DataRes, datafiles, DataDir, DFversion, Datainfo  
