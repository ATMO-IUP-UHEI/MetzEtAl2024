#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot Suppfigure 5 for Metz et al 2024
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""
import numpy as np
import pandas as pd

from functions2 import getReferenceDate, getMSC, DetrendMonthdate, getTimeMeans
from RegionParam import getRegion
from ReadDataFrames import getData
import datetime
import matplotlib.pyplot as plt


#settings
year_min ,month_min, day_min = 2015, 1, 1
year_max, month_max, day_max = 2018, 12, 31
minyearmean = 2010
timeperiod = [datetime.date(year_min,month_min,day_min),datetime.date(year_max, month_max,day_max)]
Numm = 756
offset,minnum = 384, 0

RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

Datainfo = pd.read_csv("Data_Overview.csv")
DateRef = getReferenceDate(year_min,year_max,month_min,month_max)


colors = ['grey','orange','blue','grey','grey','grey','grey','grey','red','grey','grey','grey','darkblue','darkgrey','moccasin','tan','azure','teal','crimson','olive']
        
def getMonthSum(gdf,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,ValueName,Errorname = ''):
   
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
    output = output_all[(output_all.number > 0)]

    return output

# MIP cs
MIPmodels = ['AMES','BAKER','CAMS','CMS-FLUX', 'COLA', 'CT','NIES','OU','TM5-4DVAR','UT','WOMBAT','JHU','OCO2']
MMMIPLNLGISs = []
for MIPmod in MIPmodels:
    MMMIPLNLGISs.append(getData('MIPLNLGIS_csOCO2all',MIPmod+'_CO2',timeperiod, RegionName, Numm,'MM',offset,minnum))#, data_type = 1))
#all MIP cs soundings
DS = pd.read_pickle("/mnt/data/users/eschoema/DataFrames/DFv4_MIPLNLGIS_csOCO2all_Af756.pkl")
DSdiff = DS.copy(deep = True)
DSdiff[[x + '_CO2' for x in MIPmodels]] = DSdiff[[x + '_CO2' for x in MIPmodels]]  - DS.OCO2_CO2.values[:,None]
DSsquare = DSdiff.copy(deep=True)
DSsquare[[x + '_CO2' for x in MIPmodels]+['obs_unc']] = np.square(DSsquare[[x + '_CO2' for x in MIPmodels]+['obs_unc']])


#MIP fluxes
MIPModelNames = ['Ames', 'Baker', 'CAMS','CMS-Flux', 'COLA', 'CT', 'NIES', 'OU', 'TM5-4DVAR', 'UT', 'WOMBAT','JHU']
PMIPOCOModel = []
PMIPOCOModelMeanSais = []
for Mname in MIPModelNames:
    MIPOCOModelmean = pd.read_pickle("/SavePath/MonthMeansMIP_"+Mname+"_"+RegionName+str(Numm)+".pkl")
    
    PMIPOCOModel.append(pd.merge(DateRef, MIPOCOModelmean, on=['MonthDate','Year','Month'], how = 'left'))
    PMIPOCOModelMeanSais.append(PMIPOCOModel[-1][(PMIPOCOModel[-1].Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index())
    del MIPOCOModelmean        


#TM5 fluxes
TM51x1_regionl = []
TM51x1dataset = ["flux_1x1_RemoTeC_2.4.0+IS","flux_1x1_IS","flux_1x1_ACOS+IS","flux_1x1_prior"]

for TM51x1name in TM51x1dataset:
    TM51x1 = pd.read_pickle("/SavePaths/GDF2FLUXmonthly1x1_"+TM51x1name+RegionName+str(Numm)+"V2.pkl")
    TM51x1_region_nee = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_nee_total')
    TM51x1_region_fire = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_fire_total')
    TM51x1_region = pd.merge(TM51x1_region_nee,TM51x1_region_fire, on =['Year','Month'])      
    date = TM51x1_region.apply(lambda x: datetime.date(int(x.Year),
                                                    int(x.Month),
                                                    15),axis=1)
    TM51x1_region.insert(loc=1,column='MonthDate',value=date)
    
    PTM51x1_region = pd.merge(DateRef, TM51x1_region, on=['MonthDate'], how = 'left') 
    TM51x1_regionl.append(PTM51x1_region)

allSatModel = TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_ACOS+IS')[0][0]][['MonthDate','CO2_flux_nee_total','CO2_flux_fire_total']]
allSatModel.insert(loc = 1, column = 'TM5ACOSIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
allSatModel = pd.merge(allSatModel,TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_RemoTeC_2.4.0+IS')[0][0]][['MonthDate','CO2_flux_nee_total','CO2_flux_fire_total']])
allSatModel.insert(loc = 1, column = 'TM5RTIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
allSatModel = allSatModel.merge(DateRef, on='MonthDate')
ParamList = ['TM5RTIS_nbp','TM5ACOSIS_nbp']
allSatModel.insert(0, column = 'mean_nbp', value = allSatModel[ParamList].mean(axis = 1))    
AMeanSatModel = allSatModel.groupby(['Month'])['mean_nbp'].mean().reset_index()
    
PlotVar = 'CO2_Detrend384_'+str(minnum)+'_rateNoaa'#'number'#

#only sept-nov
fig, ax = plt.subplots()    
plotMIP = MMMIPLNLGISs
markerstyles = ['v','^','<','>','1','2','3','4','.','x','+','s']
for MIPnum, MIPmod in enumerate(MIPmodels):
    if MIPmod != 'OCO2':
        MIP = getMSC(plotMIP[MIPnum],MIPmod+'_'+PlotVar)
        OCO2 = getMSC(plotMIP[np.where(np.array(MIPmodels) == 'OCO2')[0][0]],'OCO2_'+PlotVar)
        for month in [1,2,3,4]:
            #Mean Error
            #csDiff = MIP[MIP.Month == month][MIPmod+'_'+PlotVar] - OCO2[OCO2.Month == month]['OCO2_'+PlotVar]
            #RMSE
            #csDiff = np.sqrt(DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].sum()/DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].count())
            #RMSE OCO-2 weighted
            csDiff = np.sqrt((DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2']/(DSsquare[(DSsquare.Month == month)]['obs_unc'])).sum()/DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].count())
            #FluxDiff
            fluxDiff = 12/44*(AMeanSatModel[AMeanSatModel.Month == month]['mean_nbp']) - PMIPOCOModelMeanSais[MIPnum][(PMIPOCOModelMeanSais[MIPnum].Month == month)]['Landtot']
            #Abs Flux MIP
            #fluxDiff = PMIPOCOModelMeanSais[MIPnum][(PMIPOCOModelMeanSais[MIPnum].Month == month)]['Landtot']
            if colors[MIPnum] != 'grey' and month == 2:
                plt.plot(csDiff,fluxDiff.values, marker = markerstyles[month-1], ls = '',color = colors[MIPnum], label = MIPmod)
            elif MIPnum == 0:
                plt.plot(csDiff,fluxDiff.values, marker = markerstyles[month-1],color = colors[MIPnum],ls = '', label = 'Month '+str(month))
            else:
                plt.plot(csDiff,fluxDiff.values, marker = markerstyles[month-1],color = colors[MIPnum])
ax.axhline(color = 'grey',ls=':')
ax.axvline(color = 'grey',ls=':')
ax.set_ylabel('TM5-4DVar/GOSAT+IS flux - MIP/LGLNIS \n [TgC/month]')
ax.set_xlabel('RMSE errorweight MIP/LGLNIS - OCO-2 XCO2')
ax.legend()
plt.savefig('FigureS5.png', dpi = 300)

#all month
fig, ax = plt.subplots()    
plotMIP = MMMIPLNLGISs
for MIPnum, MIPmod in enumerate(MIPmodels):
    if MIPmod != 'OCO2':
        markers = '.'#markerstyles[month-1]
        MIP = getMSC(plotMIP[MIPnum],MIPmod+'_'+PlotVar)
        OCO2 = getMSC(plotMIP[np.where(np.array(MIPmodels) == 'OCO2')[0][0]],'OCO2_'+PlotVar)
        for month in range(1,13):
            #Mean Error
            #csDiff = MIP[MIP.Month == month][MIPmod+'_'+PlotVar] - OCO2[OCO2.Month == month]['OCO2_'+PlotVar]
            #RMSE
            #csDiff = np.sqrt(DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].sum()/DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].count())
            #RMSE OCO-2 weighted
            csDiff = np.sqrt((DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2']/(DSsquare[(DSsquare.Month == month)]['obs_unc'])).sum()/DSsquare[(DSsquare.Month == month)][MIPmod+'_CO2'].count())
            #flux Diff
            fluxDiff = 12/44*(AMeanSatModel[AMeanSatModel.Month == month]['mean_nbp']) - PMIPOCOModelMeanSais[MIPnum][(PMIPOCOModelMeanSais[MIPnum].Month == month)]['Landtot']
            #Abs flux MIP
            #fluxDiff = PMIPOCOModelMeanSais[MIPnum][(PMIPOCOModelMeanSais[MIPnum].Month == month)]['Landtot']
            if colors[MIPnum] != 'grey' and month == 9:
                plt.plot(csDiff,fluxDiff.values, marker = markers, ls = '',color = colors[MIPnum], label = MIPmod)
            elif MIPnum == 0:
                plt.plot(csDiff,fluxDiff.values, marker = markers,color = colors[MIPnum],ls = '')#, label = 'Month '+str(month))
            else:
                plt.plot(csDiff,fluxDiff.values, marker = markers,color = colors[MIPnum])
ax.axhline(color = 'grey',ls=':')
ax.axvline(color = 'grey',ls=':')
#ax.set_ylabel('MIP/LGLNIS flux [TgC/month]')
ax.set_ylabel('TM5-4DVar/GOSAT+IS flux - MIP/LGLNIS \n [TgC/month]')
#ax.set_xlabel('MIP/LGLNIS - OCO-2 XCO2')
ax.set_xlabel('RMSE errorweight MIP/LGLNIS - OCO-2 XCO2')
ax.legend()
#plt.savefig('FigureS5.png', dpi = 300,bbox_inches='tight')
