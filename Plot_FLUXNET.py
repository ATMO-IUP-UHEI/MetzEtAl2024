#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

#use the following three statiosn for SA
#CG-Tch, ZM-Mon, ZA-Kru

import pandas as pd
import numpy as np
import glob
import datetime
import matplotlib.pyplot as plt

station = 'CG-Tch'#'ZA-Kru'#,'CG-Tch','ZM-Mon'
filepath = glob.glob("/FLUXNET/FLX_"+station+"_FLUXNET2015_SUBSET_DD_*.csv")
df = pd.read_csv(filepath[0])

# for variable description see: https://fluxnet.org/data/fluxnet2015-dataset/subset-data-product/
# P_F: Precipitation consolidated from P and P_ERA
# P_f_QC: Quality flag for P_F, fraction between 0-1, indicating percentage of measured data
# CO2_F_MDS: CO2 mole fraction, gapfilled with MDS
if station != 'CG-Tch':
    dfsel = df[['TIMESTAMP','SW_IN_POT', 'TA_F', 'TA_F_QC','VPD_F', 'VPD_F_QC','P_F',
           'P_F_QC','CO2_F_MDS', 'CO2_F_MDS_QC','TS_F_MDS_1','TS_F_MDS_1_QC',
           'SWC_F_MDS_1', 'SWC_F_MDS_2','SWC_F_MDS_1_QC', 'SWC_F_MDS_2_QC','NEE_VUT_REF',
           'NEE_VUT_REF_QC','RECO_NT_VUT_REF','GPP_NT_VUT_REF','RECO_DT_VUT_REF','GPP_DT_VUT_REF','RECO_SR']].copy()
    dfsel.loc[dfsel[(dfsel['SWC_F_MDS_1_QC'] > 1.5)|(dfsel['SWC_F_MDS_1_QC'] < 0)].index,
                                    'SWC_F_MDS_1'] = np.nan
    dfsel.loc[dfsel[(dfsel['SWC_F_MDS_1_QC'] > 1.5)|(dfsel['SWC_F_MDS_1_QC'] < 0)].index,
                                    'SWC_F_MDS_1_QC'] = np.nan

else:
    dfsel = df[['TIMESTAMP','SW_IN_POT', 'TA_F', 'TA_F_QC','VPD_F', 'VPD_F_QC','P_F',
           'P_F_QC','CO2_F_MDS', 'CO2_F_MDS_QC','TS_F_MDS_1','TS_F_MDS_1_QC',
           'NEE_VUT_REF',
           'NEE_VUT_REF_QC','RECO_NT_VUT_REF','GPP_NT_VUT_REF','RECO_DT_VUT_REF','GPP_DT_VUT_REF','RECO_SR']].copy()

del df

for var in ['NEE_VUT_REF']:
    try:
        if dfsel[var].min() == -9999:
            dfsel.loc[dfsel[(dfsel[var] == -9999)].index,
                                    var] = np.nan
    except:
        pass

for var in ['TS_F_MDS_1']:
    try:
        if dfsel[var].min() == -9999:
            dfsel.loc[dfsel[(dfsel[var] == -9999)].index,
                                    var] = np.nan
    except:
        pass


dfsel.insert(0,column = 'Year', value = np.floor(dfsel.TIMESTAMP.values/10000))
dfsel.insert(0,column = 'Month', value = np.floor(dfsel.TIMESTAMP.values/100)-dfsel.Year*100)
dfsel.insert(0,column = 'Day', value = np.floor(dfsel.TIMESTAMP.values)-(dfsel.Year*10000+dfsel.Month*100))
dfsel.insert(0,column = 'Date',value = dfsel.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1))

for y in range(int(dfsel.Year.min()),int(dfsel.Year.max())+1):
    
    dfp = dfsel[(dfsel.Year == y)].copy(deep=True)
    print(dfp.NEE_VUT_REF.max())
    fig, ax= plt.subplots(figsize=(10,6))
    secAx = ax.twinx()
    secAx2 = ax.twinx()
    #secAx3 = ax2.twinx()
    secAx2.spines["right"].set_position(("axes", 1.12))
    secAx2.set_frame_on(True)
    secAx2.patch.set_visible(False)
    for sp in secAx2.spines.values():
        sp.set_visible(False)
    secAx2.spines["right"].set_visible(True)
    
    if False: # quality filtering
        dfp.loc[dfp[(dfp['P_F_QC'] > 0)].index,
                  'P_F'] = np.nan
        
    
    #ax.plot([dfp.Date.values[0],dfp.Date.values[1]], [0,0], color = 'grey')
    #ax.plot(dfp.Date, dfp.NEE_VUT_REF, color = 'green')
    ax.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['NEE_VUT_REF_QC'] >=0.5)][['NEE_VUT_REF','Date']], on = 'Date', how = 'left').NEE_VUT_REF, color = 'green',label = 'all meas and good Q>50%')
    ax.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['NEE_VUT_REF_QC'] <0.5)][['NEE_VUT_REF','Date']], on = 'Date', how = 'left').NEE_VUT_REF, color = 'lightgreen' , label = 'partly not good')
    #RECO_NT_VUT_REF, GPP_NT_VUT_REF, SWC_F_MDS_1
    secAx2.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['P_F_QC'] == 1)][['P_F','Date']], on = 'Date', how = 'left').P_F, color = 'blue',label = 'all meas')
    secAx2.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['P_F_QC'] <1)][['P_F','Date']], on = 'Date', how = 'left').P_F, color = 'lightblue' , label = 'partly ERA')
    #secAx.plot(dfp.Date, dfp.GPP_NT_VUT_REF, color = 'black')
    #secAx2.plot(dfp.Date, dfp.RECO_NT_VUT_REF, color = 'brown')
    
    
    if station == 'CG-Tch':
        pass
    else:
        #dfp.loc[dfp[(dfp['SWC_F_MDS_1_QC'] > 0)].index,
        #          'SWC_F_MDS_1'] = np.nan
        secAx.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['SWC_F_MDS_1_QC'] == 1)][['SWC_F_MDS_1','Date']], on = 'Date', how = 'left').SWC_F_MDS_1, color = 'brown',label = 'all meas')
        secAx.plot(dfp.Date, pd.merge(dfp['Date'],dfp[(dfp['SWC_F_MDS_1_QC'] < 1)][['SWC_F_MDS_1','Date']], on = 'Date', how = 'left').SWC_F_MDS_1, color = 'peru',label = 'partly ERA')
    
    ax.set_xlabel(r'Date',fontsize=10)
    #ax2.set_ylabel(r'$\rm CO_2~flux~($' +VegVar+r'$\rm )~[umol/(m^2s)]$',fontsize=7)
    ax.set_ylabel(r'$\rm NEE~CO_2~flux~[gC/(m^2d)]$',fontsize=10)
    #ax2.set_ylabel(r'Daily mean Fc [umol/m2/s]',fontsize=15)
    ax.yaxis.label.set_color('green')
    ax.yaxis.label.set_fontsize(10)
    
    #ax2.spines['left'].set_color('green')
    ax.tick_params(axis='y', colors='green',labelsize = 10)
    ax.legend(loc = 2)
    secAx.set_ylabel(r'Soil Water Fraction [%]',fontsize=10)
    #secAx.set_ylabel(r'GPP [umol/(m^2s)]',fontsize=10)
    secAx.yaxis.label.set_color('brown')
    secAx.yaxis.label.set_fontsize(10)
    #secAx.spines['right'].set_color('brown')
    secAx.tick_params(axis='y', colors='brown',labelsize = 10)
    secAx.legend(loc = 9)
    secAx2.set_ylabel(r'Precipitation [mm/d]',fontsize=10)
    #secAx2.set_ylabel(r'RECO [umol/(m^2s)]',fontsize=10)
    secAx2.yaxis.label.set_color('blue')
    #secAx2.spines['right'].set_color('blue')
    secAx2.tick_params(axis='y', colors='blue',labelsize = 10)
    secAx2.legend()
    #secAx2.set_ylim(secAx.get_ylim())
    plt.title(str(y))
    plt.savefig("/"+station+str(y)+"_NEEwQualityBT05_PqualityFpartlyERA_SWCqualityFpartlyERA.png",bbox_inches='tight')
    
    