#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
# Script to plot timeseries of CO2 concentrations for Metz et al 2024
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz

print('get modules')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas
import argparse
print('and functions')
from functions2 import getReferenceDate, getMSC, CombineTimeseries,addMonthDate,DetrendMonthdate, getTimeMeans
from RegionParam import getRegion
from ReadDataFrames import getData
import datetime

### start with e.g.
# python Plot_TimeSeries.py 756 '2009-01-01' '2020-12-31' 0

# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to Plot a TimeSeries for a chosen region ")
    parser.add_argument("Number", type=int, help="Enter Region ID Number")
    parser.add_argument("startdate", type=str, help="enter startdate as string in form 'yyyy-mm-dd'")
    parser.add_argument("enddate", type=str, help="enter enddate as string in form 'yyyy-mm-dd'")    
    parser.add_argument("plot_type", type=int, help="type of plot which should be created")
    args = parser.parse_args()
    return args

# SETTINGS:
print('get settings')
args = parse_arguments()
year_min = int(args.startdate[0:4])#2009
month_min = int(args.startdate[5:7])#4
day_min = int(args.startdate[8:10])#1
year_max = int(args.enddate[0:4])#2019
month_max = int(args.enddate[5:7])#12
day_max = int(args.enddate[8:10])#31
timeperiod = [datetime.date(year_min,month_min,day_min),datetime.date(year_max, month_max,day_max)]
Numm = args.Number
Plot_type = args.plot_type

RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
print(getRegion(Numm))

Datainfo = pd.read_csv("Data_Overview.csv")
DateRef = getReferenceDate(year_min,year_max,month_min,month_max)
minnum = 10
offset = 384
Var = 'CO2_DetrendMean_'+str(minnum)+'_rateNoaa'
Var2 = 'CO2_DetrendMean_0_rateNoaa'

# getData(Dataset ,ValueName, timeperiod,RegionName, Numm,FrameType, offset, minNum,rate='Noaa', **selection)
# selection: meas_geom=0, quality = 0, newDFversion = 2, in_Remotec_data, ,DFversion = 3
# see "Data_Overview.csv" for possible datasets
MMRT240 = getData('RT240_sounding','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum)#, meas_geom = 1)#,DFversion = 16)
MMACOS = getData('ACOS_sounding','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum, quality = 0)#, land_frac = 0)#,DFversion = 12)
MMOCOv11 = getData('OCO2v11_sounding','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum, quality = 0)#, land_frac = 0)#,DFversion = 12)

MMCT2022csRT240 = getData('CT2022_csRT240','CO2',timeperiod, RegionName, Numm,'MM',0,minnum,rate='Noaa')#,newDFversion=2)
MMCAMS_IScsRT240 = getData('CAMS-IS_csRT240','CO2',timeperiod, RegionName, Numm,'MM',0,minnum,rate='Noaa')
MMTM54DVar_IS3x2csACOS = getData('TM54DVarIS3x2_csACOS','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum)#minnum makes sense as cosampled
MMTM54DVar_IS3x2csRT240 = getData('TM54DVarIS3x2_csRT240','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum)
MMCOCCONG = getData('COCCONGobabeb_sounding','CO2',timeperiod, RegionName, Numm,'MM',offset,minnum,rate='Noaa')    

MMGOSAT, GOSATMSC, GOSATMSCstd, GOSATMSCDS  = CombineTimeseries([MMRT240,MMACOS],['RT240','ACOS'],Var)
MMTM5IS, TM5ISMSC, TM5ISMSCstd, TM5ISMSCDS  = CombineTimeseries([MMTM54DVar_IS3x2csRT240,MMTM54DVar_IS3x2csACOS],['TM5IScsRT','TM5IScsACOS'],Var)
MMModelcs, ModelcsMSC, ModelMSCcsstd, ModelcsMSCDS = CombineTimeseries([MMCAMS_IScsRT240,MMTM54DVar_IS3x2csRT240, MMCT2022csRT240],['CAMS-IScsRT240','TM5-4DVarcsRT240','CT2022csRT240'],Var)





if Plot_type == 4: # paper figure 1
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8#0.8
    cs = 1.3
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
    
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.98*cm),
                            gridspec_kw={'width_ratios': [20, 4],'height_ratios':[1]})
    
    PlotVar = 'CO2_Detrend384_'+str(minnum)+'_rateNoaa'#'number'#
    ylabel = r'$\rm Detrended~CO_{2}~[ppm]$'
    
    #Panel A
    #Ensembles
    ax2[0].fill_between(MMGOSAT.MonthDate, MMGOSAT.minValue, MMGOSAT.maxValue,color = 'red',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[0].plot(MMGOSAT.MonthDate, MMGOSAT.meanValue,color = 'red',ls = '-',linewidth= lw,marker = '',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    ax2[0].fill_between(MMModelcs.MonthDate, MMModelcs.minValue, MMModelcs.maxValue,color = 'blue',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[0].plot(MMModelcs.MonthDate, MMModelcs.meanValue,color = 'blue',ls = '-',linewidth= lw,marker = '',markersize = 2,zorder = 2, label = r'$\rm Inverse~Models_{IS,cs}$')# $\rm CO_2$ ')
    
    #Panel B
    ax2[1].plot(range(1,13), GOSATMSC.meanValue.iloc[indexrange],color = 'red',linewidth = lw,ls = '-',marker = '', label = r'GOSAT')# $\rm CO_2$ ')
    ax2[1].fill_between(range(1,13), GOSATMSC.meanValue.iloc[indexrange]-GOSATMSCstd.meanValue.iloc[indexrange],GOSATMSC.meanValue.iloc[indexrange]+GOSATMSCstd.meanValue.iloc[indexrange],color = 'red',alpha = 0.3)# $\rm CO_2$ ')
    ax2[1].plot(range(1,13), ModelcsMSC.meanValue.iloc[indexrange],color = 'blue',linewidth = lw,ls = '-',marker = '', label = r'$\rm Inverse~Models_{IS,cs}$')# $\rm CO_2$ ')
    ax2[1].fill_between(range(1,13), ModelcsMSC.meanValue.iloc[indexrange]-ModelMSCcsstd.meanValue.iloc[indexrange],ModelcsMSC.meanValue.iloc[indexrange]+ModelMSCcsstd.meanValue.iloc[indexrange],linewidth= lw,color = 'blue',alpha = 0.3)# $\rm CO_2$ ')
    
    

     # SETTINGS
    ax2[0].text(0.03, 0.95, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=9)#, weight='bold')
    ax2[0].legend(fontsize=6,loc = 1,ncol = 2)#,loc=3)#, loc = 9)
    ax2[0].set_ylabel(ylabel,fontsize=7)
    ax2[0].tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2[0].tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2[0].grid(True,which = 'major', axis='x',zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    #ax2[0].set_ylim(-280,280)
    ax2[0].set_xlim(datetime.date(year_min-1,12,15),datetime.date(year_max+1,3,15))
    ax2[0].set_axisbelow(True)

    ax2[1].text(0.15, 0.95, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=9)#, weight='bold')      
    ax2[1].set_xticks(ticks = [3,6,9,12])
    #ax2[1].set_xticklabels(['3','6','9','12'])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1].grid(True,which = 'major', axis='y',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    ax2[1].set_ylim(ax2[0].get_ylim())
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6, zorder = 0)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 

    plt.savefig("Figure1.png", dpi=400, bbox_inches = "tight")#,format = 'pdf')

if Plot_type == 2: # paper appendix 1
    MMGOSAT, GOSATMSC, GOSATMSCstd, GOSATMSCDS  = CombineTimeseries([MMRT240,MMACOS],['RT240','ACOS'],Var,2015,2018)
    MMModelcs, ModelcsMSC, ModelMSCcsstd, ModelcsMSCDS = CombineTimeseries([MMCAMS_IScsRT240,MMTM54DVar_IS3x2csRT240, MMCT2019BcsRT240],['CAMS-IScsRT240','TM5-4DVarcsRT240','CT2019BcsRT240'],Var,2015,2018)
    
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8#0.8
    cs = 1.3
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.98*cm),
                            gridspec_kw={'width_ratios': [20, 4],'height_ratios':[1]})
    
    PlotVar = datavar = 'CO2_DetrendMean_'+str(minnum)+'_rateNoaa'#'number'#
    #PlotVar = datavar = 'CO2_Detrend384_'+str(minnum)+'_rateNoaa'#'number'#
    PlotVar2 = datavar2 = 'CO2_DetrendMean_0_rateNoaa'#'number'#
    ylabel = r'$\rm Detrended~CO_{2}~[ppm]$'
    #Panel A
    #Ensembles
    ax2[0].fill_between(MMGOSAT.MonthDate, MMGOSAT.minValue, MMGOSAT.maxValue,linewidth = lw,color = 'red',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[0].plot(MMGOSAT.MonthDate, MMGOSAT.meanValue,color = 'red',ls = '-',linewidth= lw,marker = '',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    
    ax2[0].fill_between(MMModelcs.MonthDate, MMModelcs.minValue, MMModelcs.maxValue,color = 'blue',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[0].plot(MMModelcs.MonthDate, MMModelcs.meanValue,color = 'blue',ls = '-',linewidth= lw,marker = '',markersize = 2,zorder = 2, label = r'$\rm Inverse~Models_{IS,cs}$')# $\rm CO_2$ ')
    
    ax2[0].plot(MMOCOv11.MonthDate, MMOCOv11[PlotVar],color = 'black',ls = '-', marker = '',linewidth = lw, label = r'OCO-2')
    #cs
    ax2[0].plot(MMTM54DVar_IS3x2csRT240.MonthDate, MMTM54DVar_IS3x2csRT240[PlotVar],color = 'darkblue',ls = ':',linewidth= lw, marker = '', label = r'$\rm TM5-4DVar_{IS,cs}$')
    ax2[0].plot(MMCT2022csRT240.MonthDate, MMCT2022csRT240[PlotVar],color = 'darkblue',linewidth= lw,ls = '--', marker = '', label = r'$\rm CT2019B_{IS,cs}$')
    ax2[0].plot(MMCAMS_IScsRT240.MonthDate, MMCAMS_IScsRT240[PlotVar],color = 'darkblue',linewidth= lw,ls = '-.', marker = '', label = r'$\rm CAMS_{IS,cs}$')
    
    ax2[1].plot(range(1,13), GOSATMSC.meanValue.iloc[indexrange],color = 'red',ls = '-',linewidth= lw,marker = '', label = r'GOSAT')# $\rm CO_2$ ')
    ax2[1].fill_between(range(1,13), GOSATMSC.meanValue.iloc[indexrange]-GOSATMSCstd.meanValue.iloc[indexrange],GOSATMSC.meanValue.iloc[indexrange]+GOSATMSCstd.meanValue.iloc[indexrange],linewidth= lw,color = 'red',alpha = 0.3)# $\rm CO_2$ ')
    
    ax2[1].plot(range(1,13), ModelcsMSC.meanValue.iloc[indexrange],color = 'blue',linewidth = lw,ls = '-',marker = '', label = r'$\rm Inverse~Models_{IS,cs}$')# $\rm CO_2$ ')
    ax2[1].fill_between(range(1,13), ModelcsMSC.meanValue.iloc[indexrange]-ModelMSCcsstd.meanValue.iloc[indexrange],ModelcsMSC.meanValue.iloc[indexrange]+ModelMSCcsstd.meanValue.iloc[indexrange],color = 'blue',alpha = 0.3)# $\rm CO_2$ ')
    
    ax2[1].plot(range(1,13), getMSC(MMOCOv11,datavar,2015,2018)[datavar].iloc[indexrange],color = 'black',ls = '-',linewidth = lw,marker = '', label = r'OCO-2')# $\rm CO_2$ ')
    ax2[1].plot(range(1,13), getMSC(MMTM54DVar_IS3x2csRT240,datavar,2015,2018)[datavar].iloc[indexrange],color = 'darkblue',ls = ':',linewidth = lw,marker = '')
    ax2[1].plot(range(1,13), getMSC(MMCT2022csRT240,datavar,2015,2018)[datavar].iloc[indexrange],color = 'darkblue',ls = '--',linewidth = lw,marker = '')
    ax2[1].plot(range(1,13), getMSC(MMCAMS_IScsRT240,datavar,2015,2018)[datavar].iloc[indexrange],color = 'darkblue',ls = '-.',linewidth = lw,marker = '')

     # SETTINGS
    ax2[0].text(0.03, 0.95, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=9)#, weight='bold')
    ax2[0].legend(fontsize=6,loc = 4,ncol = 4)#,loc=3)#, loc = 9)
    ax2[0].set_ylabel(ylabel,fontsize=7)
    ax2[0].tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2[0].tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2[0].grid(True,which = 'major', axis='x',zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    #ax2[0].set_ylim(-280,280)
    ax2[0].set_xlim(datetime.date(year_min-1,12,15),datetime.date(year_max+1,3,15))
    #ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2022,12,31))
    ax2[0].set_axisbelow(True)

    ax2[1].text(0.15, 0.95, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=9)#, weight='bold')      
    ax2[1].set_xticks(ticks = [3,6,9,12])
    #ax2[1].set_xticklabels(['3','6','9','12'])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1].grid(True,which = 'major', axis='y',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    ax2[1].set_ylim(ax2[0].get_ylim())
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6, zorder = 0)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 

    plt.savefig("FigureS1.png", dpi=400, bbox_inches = "tight")#,format = 'pdf')


if Plot_type == 3: #Supp Figure 3
    
    PlotVar = 'CO2_DetrendMean_'+str(minnum)+'_rateNoaa'
    PlotVarg = 'CO2_DetrendMean_0_rateNoaa'
    ylabel = r'$\rm Detrended~CO_2~(ppm)$'
    
    fig2, ax2 = plt.subplots(figsize=(12,5))
    #red, darkred
    ax2.fill_between(MMGOSAT.MonthDate, MMGOSAT.minValue, MMGOSAT.maxValue,color = 'red',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2.plot(MMGOSAT.MonthDate, MMGOSAT.meanValue,color = 'red',ls = '-',marker = '.',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    
    #blue
    ax2.fill_between(MMModelcs.MonthDate, MMModelcs.minValue, MMModelcs.maxValue,color = 'blue',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2.plot(MMModelcs.MonthDate, MMModelcs.meanValue,color = 'blue',ls = '-',linewidth= 1,marker = '.',markersize = 2,zorder = 2, label = r'$\rm Inverse~Models_{IS,cs}$')# $\rm CO_2$ ')
    
    ax2.plot(MMCOCCONG.MonthDate, MMCOCCONG[PlotVar],color = 'black',ls = ':', marker = '.', label = r'COCCON Gobabeb')
    
    ax2.legend()
    #ax2.set_yscale('log')
    ax2.set_ylabel(ylabel,fontsize=15)
    ax2.tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 15)
    ax2.tick_params(axis = 'y',labelsize = 15,zorder = 8)
    ax2.set_xlabel(r'Date',fontsize=15)
    #ax2.set_xlim(datetime.date(2009,9,15),datetime.date(2010,3,31))
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_axisbelow(True)
    plt.savefig("SuppFigure3.png", dpi=250, bbox_inches = "tight")

