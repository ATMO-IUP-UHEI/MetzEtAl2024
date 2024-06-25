#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
Script to investigate Carbon Fluxes 

"""
import datetime
from tarfile import REGULAR_TYPES
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import getMeanAmplitude, getMonthMeans, getMonthSum, addMonthDate, EarlyGPPModelNames, GoodModelsNames, BadModelNames, LateRespNames,  getReferenceDate,  getMonthSumFCday,getNumDayOfMonth,getEnsembleFluxMIP, getMSC, NormalizeMeanCycle
import glob
import xarray as xr
import os.path
from RegionParam import getRegion
from CreateAndSaveGeoDataframeCTFlux1x1 import CreateDataFrameCTFluxMonthly1x1
from CreateAndSaveGeoDataframeFLUXCOM import CreateDataFrameSIFFC
from CreateAndSaveGeoDataframeCAMSFlux1x1 import CreateDataFrameCAMSflux, CombineGDFCAMS
from CreateAndSaveGeoDataframeMIP import CreateDataFrameMIPflux
from CreateAndSaveGeoDataFrameSIFGOMEv2 import CreateDataFrameSIFGOMEv2
from CreateAndSaveGeoDataframeGFAS import CreateDataframeGFAS
from CreateAndSaveGeoDataframeFINN import CreateDataFrameFINN

import matplotlib.font_manager as fm

font = {'family' : 'Arial'}
mpl.rc('font', **font)

WithOCO2 = False

#choose which datasets are getting loaded
CAMS_b = True###True###
CT_flux2022 = True#####
SIF_FC = False#True#
TM5Flux1x1 = True###True###
TRENDY = True
Fires = True#True###True###
MIPdata = False#True###True###
MIPISdata = False###
MIPTM54DVardata = False#True
MIPOCOModelsdata = False#True
MIPpriorModelsdata = False#True
MIPISModelsdata = False#True
GOMESIFv2 = False#True

#further specifications, which TRENDY models are taken into the analyses
TrendyModels = ['DLEM', 'IBIS', 'ISAM', 'JSBACH', 
                'LPX-Bern', 'OCN','ORCHIDEE', 'ORCHIDEE-CNP', 'ORCHIDEEv3', 
                'CLM5.0','ISBA-CTRIP','LPJ','SDGVM',
                'VISIT','YIBs','CABLE-POP', 'CLASSIC','JULES-ES-1p0'] #'JULES-ES-1p0'


#select the Inversion TM5-4DVar datasets
TM51x1dataset = ["flux_1x1_RemoTeC_2.4.0+IS","flux_1x1_IS","flux_1x1_ACOS+IS","flux_1x1_prior","flux_1x1_LNLGIS"]#,"flux_1x1_LNLGIS","flux_1x1_LNLG"]#,"flux_1x1_ACOS","flux_1x1_prior","flux_1x1_RemoTeC"]

#input (could be moved to arg parser)
# which region and which period
Numm = 756
#RegionArea used in TRENDY models and MeanCycleTM5
RegionAreaDic = {756: 6.3025380864*10**12,
                 758:3.910027*10**12, 
                 766: 1.0212565*10**13, 
                 769:3517620413886,
                 770:2784917672538, 
                 771: 2545922571910,
                 772: 3717600157633,
                 784: 2333898024102,
                 949: 7.965*10**12 ,
                 950: 5422341009727 ,
                 951 : 2525018121358 ,
                 952 : 4902147095848 ,
                 953 : 3045212035237}

RegionArea = RegionAreaDic[Numm]#Af 756 6.3025380864*10**12#SA 766 1.0212565*10**13## AU 949 7.965*10**12 ##DryERA5 (950):5422341009727 #949: 7.965*10**12 #wetERA5 (951): 2525018121358 #dryLC(952): 4902147095848 #wetLC (953): 3045212035237
startdate = '2009-04-15'#'2009-04-15'#'2009-05-01'#'2009-01-01'#'2015-01-01'#'2009-04-15'#'2014-08-15'#
enddate = '2018-12-31'
minyearmean = 2010#2015 #needs to be changed for annual flux calculation later starting than 2009
                   # exception: the june-july calculations are starting mid of 2009
    
#main settings
year_min = int(startdate[0:4])#2009
month_min = int(startdate[5:7])#4
day_min = int(startdate[8:10])#1
year_max = int(enddate[0:4])#2019
month_max = int(enddate[5:7])#12
day_max = int(enddate[8:10])#31

RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

Transcom = pd.read_pickle("CTTranscomMask1x1Borders.pkl")

DateRef = getReferenceDate(year_min,year_max,month_min,month_max)



# get the data sets
# --------------------------------------------------------
if Fires:
    #read GFED
    # if not available yet, create via: Plot_TimeSeries_Detmers_CompareData.py
    GFEDregion = pd.read_pickle("/GFED4/MonthMeansCO2_"+RegionName+str(Numm)+".pkl")
    PGFEDregion = pd.merge(DateRef, GFEDregion, on=['MonthDate'], how = 'left')
    PGFEDregion.rename(columns={'Month_x':'Month','Year_x':'Year'},inplace=True)
    # get annual mean and Stdev
    GFEDMonthlyAnnualMean = PGFEDregion.groupby('Month')['total_emission'].mean().reset_index()
    GFEDMonthlyAnnualStd = PGFEDregion.groupby('Month')['total_emission'].std(ddof=0).reset_index()
    #get mean Amplitude
    GFEDmeanAmplitude, GFEDmeanAmplitudeStDev = getMeanAmplitude(PGFEDregion,'total_emission','Year',2009,2018)[0:2]
    GFEDmeanAmplitude_no09, GFEDmeanAmplitudeStDev_no09 = getMeanAmplitude(PGFEDregion,'total_emission','Year',2010,2018)[0:2]
    GFEDmeanAmplitudejj, GFEDmeanAmplitudeStDevjj,GFEDmeanAmplitudeMin,GFEDmeanAmplitudeMax = getMeanAmplitude(PGFEDregion,'total_emission','JJYear',2009,2017)

    #read GFAS
    if not os.path.isfile("/GFAS/dataframes/GDF2_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameGFAS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Month')

    gdfGFAS = pd.read_pickle("/GFAS/dataframes/GDF2_"+RegionName+str(Numm)+".pkl")

    if os.path.isfile("/GFAS/dataframes/monthFrames/MonthMeans_GFAS_"+RegionName+ str(Numm)+".pkl") and os.path.isfile("/GFAS/dataframes/monthFrames/MonthMeans_GFASten_"+RegionName+ str(Numm)+".pkl"):
        Month_means_GFAS = pd.read_pickle("/GFAS/dataframes/monthFrames/MonthMeans_GFAS_"+RegionName+ str(Numm)+".pkl")

    else:
        Month_means_GFASCO21 = getMonthSum(gdfGFAS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2fireE')
        Month_means_GFAS = addMonthDate(Month_means_GFASCO21)
        
    PMonth_means_GFAS = pd.merge(DateRef, Month_means_GFAS, on=['MonthDate'], how = 'left')

    #read FINN
    if not os.path.isfile("/FINN/dataframes/GDF1_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameFINN(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Month')

    gdfFINN = pd.read_pickle("/FINN/dataframes/GDF1_"+RegionName+str(Numm)+".pkl")

    if os.path.isfile("/FINN/dataframes/monthFrames/MonthMeans_FINN_"+RegionName+ str(Numm)+".pkl"):
        Month_means_FINN = pd.read_pickle("/FINN/dataframes/monthFrames/MonthMeans_FINN_"+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_FINNCO21 = getMonthSum(gdfFINN,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2fireE')
        Month_means_FINN = addMonthDate(Month_means_FINNCO21)

    PMonth_means_FINN = pd.merge(DateRef, Month_means_FINN, on=['MonthDate'], how = 'left')



if MIPdata:
    if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl"):
        MIPmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl")
    else:
        if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+RegionName+str(Numm)+".pkl"):
            MIP = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,version = 'regular',ModelNames=[])
            MIP = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+RegionName+str(Numm)+".pkl")
    
        MIPmean = getMonthSum(MIP,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
        #gaussian errorpropagation for Std: reg Std = sqrt(sum(grid std^2))
        MIPStd = getMonthMeans(MIP,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot','LandStdtot')
        # * number of grid cell as getMonthMeans divides Std by that
        MIPStd.insert(loc = 1, column = 'LandStdtot2',value= MIPStd.LandStdtot*MIPStd.number)
        MIPmonth = pd.merge(MIPmean,MIPStd[['LandStdtot2','Month','Year']], on=['Month','Year'], how = 'left')
        # insert MontDate in datframe
        md = MIPmonth.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
        MIPmonth.insert(loc=1,column='MonthDate',value=md)
        MIPmonth.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl")
        del MIPStd,MIPmean,MIP,md

    PMIP = pd.merge(DateRef, MIPmonth, on=['MonthDate','Month','Year'], how = 'left')
    PMIPMeanSais = PMIP[(PMIP.MonthDate >= datetime.date(minyearmean,1,1))].groupby(['Month'])['Landtot'].mean().reset_index()
    PMIPMeanSaisStDev = PMIP[(PMIP.MonthDate >= datetime.date(minyearmean,1,1))].groupby(['Month'])['Landtot'].std(ddof = 0).reset_index()
    
    del MIPmonth



if MIPISdata:
    if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl"):
        MIPISmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_IS_"+RegionName+str(Numm)+".pkl")
    else:
        if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_IS_"+RegionName+str(Numm)+".pkl"):
            MIPIS = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_IS_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,version = 'regularIS')
            MIPIS = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_IS_"+RegionName+str(Numm)+".pkl")
    
        MIPISmean = getMonthSum(MIPIS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
        #gaussian errorpropagation for Std: reg Std = sqrt(sum(grid std^2))
        MIPISStd = getMonthMeans(MIPIS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot','LandStdtot')
        # * number of grid cell as getMonthMeans divides Std by that
        MIPISStd.insert(loc = 1, column = 'LandStdtot2',value= MIPISStd.LandStdtot*MIPISStd.number)
        
        MIPISmonth = pd.merge(MIPISmean,MIPISStd[['LandStdtot2','Month','Year']], on=['Month','Year'], how = 'left')
        md = MIPISmonth.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
        MIPISmonth.insert(loc=1,column='MonthDate',value=md)
        MIPISmonth.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_IS_"+RegionName+str(Numm)+".pkl")
        del MIPISStd,MIPISmean,MIPIS,md


    PMIPIS = pd.merge(DateRef, MIPISmonth, on=['MonthDate','Year','Month'], how = 'left')
    PMIPISMeanSais = PMIPIS[(PMIPIS.Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index()
    PMIPISMeanSaisStDev = PMIPIS[(PMIPIS.Year >= minyearmean)].groupby(['Month'])['Landtot'].std(ddof = 0).reset_index()
    del MIPISmonth

if MIPTM54DVardata:
    if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl"):
        MIPTM54DVarmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
    else:
        if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl"):
            MIPTM54DVar = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,version = 'TM5-4DVar')
            MIPTM54DVar = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
    
        MIPTM54DVarmean = getMonthSum(MIPTM54DVar,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
        
        md = MIPTM54DVarmean.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
        MIPTM54DVarmean.insert(loc=1,column='MonthDate',value=md)
        MIPTM54DVarmean.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
        del MIPTM54DVar,md


    PMIPTM54DVar = pd.merge(DateRef, MIPTM54DVarmean, on=['MonthDate','Year','Month'], how = 'left')
    PMIPTM54DVarMeanSais = PMIPTM54DVar[(PMIPTM54DVar.Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index()
    del MIPTM54DVarmean

if MIPOCOModelsdata:
    PMIPOCOModel = []
    PMIPOCOModelMeanSais = []
    MIPModelNames = ['UT','WOMBAT','OU','CSU','NIES','CAMS','LoFI','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames']
    for Mname in MIPModelNames:
        if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl"):
            MIPTM54DVarmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
        else:
            if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_"+RegionName+str(Numm)+".pkl"):
                MIPOCOModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_"+RegionName+str(Numm)+".pkl")
            else:
                CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, Mname,MIPModelNames)
                MIPOCOModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_"+RegionName+str(Numm)+".pkl")
        
            MIPOCOModelmean = getMonthSum(MIPOCOModel,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
            
            md = MIPOCOModelmean.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
            MIPOCOModelmean.insert(loc=1,column='MonthDate',value=md)
            MIPOCOModelmean.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+Mname+"_"+RegionName+str(Numm)+".pkl")
            del MIPOCOModel,md


        PMIPOCOModel.append(pd.merge(DateRef, MIPOCOModelmean, on=['MonthDate','Year','Month'], how = 'left'))
        PMIPOCOModelMeanSais.append(PMIPOCOModel[-1][(PMIPOCOModel[-1].Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index())
        del MIPOCOModelmean

if MIPISModelsdata:
    PMIPISModel = []
    PMIPISModelMeanSais = []
    MIPModelNamesIS = ['UT','WOMBAT','OU','CSU','NIES','CAMS','LoFI','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames']
    for Mname in MIPModelNamesIS:
        if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl"):
            MIPTM54DVarmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
        else:
            if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_IS_"+RegionName+str(Numm)+".pkl"):
                MIPISModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_IS_"+RegionName+str(Numm)+".pkl")
            else:
                CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, Mname+'_IS',MIPModelNamesIS)
                MIPISModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_IS_"+RegionName+str(Numm)+".pkl")
        
            MIPISModelmean = getMonthSum(MIPISModel,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
            
            md = MIPISModelmean.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
            MIPISModelmean.insert(loc=1,column='MonthDate',value=md)
            MIPISModelmean.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+Mname+"_IS_"+RegionName+str(Numm)+".pkl")
            del MIPISModel,md


        PMIPISModel.append(pd.merge(DateRef, MIPISModelmean, on=['MonthDate','Year','Month'], how = 'left'))
        PMIPISModelMeanSais.append(PMIPISModel[-1][(PMIPISModel[-1].Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index())
        del MIPISModelmean

if MIPpriorModelsdata:
    PMIPpriorModel = []
    PMIPpriorModelMeanSais = []
    MIPModelNamesPrior = ['UT','WOMBAT','OU','CSU','NIES','CAMS','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames']
    for Mname in MIPModelNamesPrior:
        if False:#os.path.isfile("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl"):
            MIPTM54DVarmean = pd.read_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_TM5-4DVar_"+RegionName+str(Numm)+".pkl")
        else:
            if os.path.isfile("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_prior_"+RegionName+str(Numm)+".pkl"):
                MIPpriorModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_prior_"+RegionName+str(Numm)+".pkl")
            else:
                CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, Mname+'_prior',MIPModelNamesPrior)
                MIPpriorModel = pd.read_pickle("/OCO-2_v10_MIP/dataframes/GDF2MIP_"+Mname+"_prior_"+RegionName+str(Numm)+".pkl")
        
            MIPpriorModelmean = getMonthSum(MIPpriorModel,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
            
            md = MIPpriorModelmean.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
            MIPpriorModelmean.insert(loc=1,column='MonthDate',value=md)
            MIPpriorModelmean.to_pickle("/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+Mname+"_prior_"+RegionName+str(Numm)+".pkl")
            del MIPpriorModel,md


        PMIPpriorModel.append(pd.merge(DateRef, MIPpriorModelmean, on=['MonthDate','Year','Month'], how = 'left'))
        PMIPpriorModelMeanSais.append(PMIPpriorModel[-1][(PMIPpriorModel[-1].Year >= minyearmean)].groupby(['Month'])['Landtot'].mean().reset_index())
        del MIPpriorModelmean



if CAMS_b:
    if False:#os.path.isfile("/CAMS/dataframes/monthFrames/MonthMeansFluxSat_"+RegionName+str(Numm)+".pkl"):
        CAMS_region = pd.read_pickle("/CAMS/dataframes/monthFrames/MonthMeansFluxSat_"+RegionName+str(Numm)+".pkl")
    else:
        if os.path.isfile("/CAMS/dataframes/GDF21x1FLUXSat_"+RegionName+str(Numm)+".pkl"):
            #CAMS = pd.read_pickle("/CAMS/dataframes/GDF3FLUXSat_"+RegionName+str(Numm)+".pkl")
            CAMS = pd.read_pickle("/CAMS/dataframes/GDF21x1FLUXSat_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Datatype='Sat')
            CombineGDFCAMS(RegionName, Numm,'Sat')
            #CAMS = pd.read_pickle("/CAMS/dataframes/GDF3FLUXSat_"+RegionName+str(Numm)+".pkl")
            CAMS = pd.read_pickle("/CAMS/dataframes/GDF21x1FLUXSat_"+RegionName+str(Numm)+".pkl")
    
        #CAMS_region = CAMS[(CAMS.Long <= Long_max) & (CAMS.Long > Long_min) & (CAMS.Lat <= Lat_max) & (CAMS.Lat > Lat_min)]
        CAMS_region_b = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_bio_tot')
        CAMS_region_o = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_ocean_tot')
        CAMS_region_ba = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_bio_tot')
        CAMS_region_oa = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_ocean_tot')
        CAMS_region_ff = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_foss_tot')

        CAMS_region = pd.merge(CAMS_region_b,CAMS_region_o, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ba, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_oa, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ff, on=['Month','Year'], how = 'left')
        
        
        md = []
        for num,row2 in CAMS_region.iterrows():
            md.append(datetime.date(int(row2.Year),int(row2.Month),15))
        CAMS_region.insert(loc=1,column='MonthDate',value= md)
        #CAMS_region.to_pickle("/CAMS/dataframes/monthFrames/MonthMeansFluxSat_"+RegionName+str(Numm)+".pkl")
        del CAMS_region_b,CAMS_region_ba,CAMS_region_o,CAMS_region_oa,CAMS_region_ff


    PCAMS_region_sat = pd.merge(DateRef, CAMS_region, on=['MonthDate'], how = 'left')
    del CAMS_region
    
    if False:#os.path.isfile("/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl"):
        CAMS_region = pd.read_pickle("/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl")
    else:
        CAMSGDFversion = '1x1'#'1x1'#'2'
        print("/CAMS/dataframes/GDF2"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl")
        if os.path.isfile("/CAMS/dataframes/GDF2"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl"):
            CAMS = pd.read_pickle("/CAMS/dataframes/GDF2"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Datatype='Sur')
            CombineGDFCAMS(RegionName, Numm,'Sur')
            CAMS = pd.read_pickle("/CAMS/dataframes/GDF2"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl")
    
        #CAMS_region = CAMS[(CAMS.Long <= Long_max) & (CAMS.Long > Long_min) & (CAMS.Lat <= Lat_max) & (CAMS.Lat > Lat_min)]
        CAMS_region_b = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_bio_tot')
        CAMS_region_o = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_ocean_tot')
        CAMS_region_ba = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_bio_tot')
        CAMS_region_oa = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_ocean_tot')
        CAMS_region_ff = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_foss_tot')

        CAMS_region = pd.merge(CAMS_region_b,CAMS_region_o, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ba, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_oa, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ff, on=['Month','Year'], how = 'left')
        
        md = []
        for num,row2 in CAMS_region.iterrows():
            md.append(datetime.date(int(row2.Year),int(row2.Month),15))
        CAMS_region.insert(loc=1,column='MonthDate',value= md)
        CAMS_region.to_pickle("/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl")

    PCAMS_region_Sur = pd.merge(DateRef, CAMS_region, on=['MonthDate'], how = 'left')

if TM5Flux1x1:
    TM51x1_regionl = []
    for TM51x1name in TM51x1dataset:
        try:
            TM51x1 = pd.read_pickle("/TM5Inversion/glb3x2_20220413/new_gridded_flux/dataframes/GDF2FLUXmonthly1x1_"+TM51x1name+RegionName+str(Numm)+"V2.pkl")
        except:
            try:
                TM51x1 = pd.read_pickle("/TM5Inversion/glb3x2_20220413/gridded_flux/dataframes/GDF2FLUXmonthly1x1_"+TM51x1name+RegionName+str(Numm)+".pkl")
            except:
                print('dataframes not created yet, please do so by using CreateAndSaveGeoDataframeTM5Flux1x1.py')
        try:
            TM51x1_region_nee = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_nee_total')
            TM51x1_region_fire = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_fire_total')
            TM51x1_region = pd.merge(TM51x1_region_nee,TM51x1_region_fire, on =['Year','Month'])      
        except:
            TM51x1_region = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_nbe_total')
            
        date = TM51x1_region.apply(lambda x: datetime.date(int(x.Year),
                                                        int(x.Month),
                                                        15),axis=1)
        TM51x1_region.insert(loc=1,column='MonthDate',value=date)
        
        PTM51x1_region = pd.merge(DateRef, TM51x1_region, on=['MonthDate'], how = 'left') 
        TM51x1_regionl.append(PTM51x1_region)


if CT_flux2022:
    if True:#new
        if os.path.isfile("/CT2022/dataframes/GDF2FLUXmonthly1x1_"+RegionName+str(Numm)+".pkl"):
            CTF22 = pd.read_pickle("/CT2022/dataframes/GDF2FLUXmonthly1x1_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameCTFluxMonthly1x1(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Version='2022')
            CTF22 = pd.read_pickle("/CT2022/dataframes/GDF2FLUXmonthly1x1_"+RegionName+str(Numm)+".pkl")
                
        CTF22_region_b = getMonthSum(CTF22,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Biotot')
        CTF22_region_o = getMonthSum(CTF22,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Ocntot')
        CTF22_region_ff = getMonthSum(CTF22,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'fftot')
        CTF22_region_f = getMonthSum(CTF22,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'firetot')

        CTF22_region = pd.merge(CTF22_region_b,CTF22_region_o, on=['Month','Year'], how = 'left')
        CTF22_region = pd.merge(CTF22_region,CTF22_region_ff, on=['Month','Year'], how = 'left')
        CTF22_region = pd.merge(CTF22_region,CTF22_region_f, on=['Month','Year'], how = 'left')
        date = CTF22_region.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     15),axis=1)
        CTF22_region.insert(loc=1,column='MonthDate',value=date)

    PCTF22_region = pd.merge(DateRef, CTF22_region, on=['MonthDate'], how = 'left')   

if SIF_FC:
    PSIFFC_region = []
    PSIFFC_region2 = []
    for Datatype in ['NEE','GPP','TER']:#,'GPP']:#['GPP','NEE']:#,'TER']:
        print("/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
        if os.path.isfile("/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl"):
            SIFFC_region = pd.read_pickle("/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
        else:
            print('get dataframe')
            if os.path.isfile("/SIF/FLUXCOM/GDF3_"+Datatype+RegionName+str(Numm)+".pkl"):
                print('yes')
                SIFFC = pd.read_pickle("/SIF/FLUXCOM/GDF3_"+Datatype+RegionName+str(Numm)+".pkl")
            else:
                CreateDataFrameSIFFC(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Datatype)
                SIFFC = pd.read_pickle("/SIF/FLUXCOM/GDF3_"+Datatype+RegionName+str(Numm)+".pkl")
        
            SIFFC_region = getMonthSum(SIFFC,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,Datatype+'tot',Datatype+'_madtot')
            SIFFC_region2 = getMonthSumFCday(SIFFC,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,Datatype)
            #print(SIFFC_region.keys())
            md = []
            md = SIFFC_region.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)  
            SIFFC_region.insert(loc=1,column='MonthDate',value= md)
            SIFFC_region2.insert(loc=1,column='MonthDate',value= md)
            SIFFC_region.to_pickle("/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
            SIFFC_region2.to_pickle("/SIF/FLUXCOM/MonthMeans2_"+Datatype+RegionName+str(Numm)+".pkl")        
        print('got FC' + Datatype)
        PSIFFC_region.append(pd.merge(DateRef, SIFFC_region, on=['MonthDate','Year','Month'], how = 'left'))   
        #PSIFFC_region2.append(pd.merge(DateRef, SIFFC_region2, on=['MonthDate'], how = 'left'))   
    FCGFED = pd.merge(PSIFFC_region[0],PGFEDregion,on=['MonthDate','Year','Month'])
    FCGFED.insert(loc=1, column = 'nbp',value=FCGFED.NEEtot/10**12*(44/12)+FCGFED.total_emission)
    FCGFED = pd.merge(DateRef, FCGFED, on=['MonthDate','Year','Month'], how = 'left')
    FCGFED.to_pickle("/SIF/FLUXCOM/MonthMeans_FCGFED_"+RegionName+str(Numm)+".pkl")#_newForSuppPaper.pkl")
    #get mean Amplitude
    FCGFEDmeanAmplitude, FCGFEDmeanAmplitudeStDev = getMeanAmplitude(FCGFED,'nbp','Year',2009,2018)[0:2]
    FCGFEDmeanAmplitude_no09, FCGFEDmeanAmplitudeStDev_no09 = getMeanAmplitude(FCGFED,'nbp','Year',2010,2018)[0:2]
    FCGFEDmeanAmplitudejj, FCGFEDmeanAmplitudeStDevjj,FCGFEDmeanAmplitudeMin,FCGFEDmeanAmplitudeMax = getMeanAmplitude(FCGFED,'nbp','JJYear',2009,2017)
    #print('FCGFED')
    #print(str(FCGFEDmeanAmplitudejj/44*12) +' +- '+ str(FCGFEDmeanAmplitudeStDevjj/44*12))
        
    #print('This are the mean annual FCGFED fluxes')
    FCannualMean = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Year'])['nbp'].sum().mean()
    FCannualStDev = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Year'])['nbp'].sum().std(ddof = 0)
    FCMonthlyAnnualMean = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Month'])['nbp'].mean().reset_index()
    FCMonthlyAnnualStDev = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Month'])['nbp'].std(ddof = 0).reset_index()
    FCMonthlyAnnualMeanNEE = PSIFFC_region[0][(PSIFFC_region[0].Year >= minyearmean)].groupby(['Month'])['NEEtot'].mean().reset_index()
    FCMonthlyAnnualStDevNEE = PSIFFC_region[0][(PSIFFC_region[0].Year >= minyearmean)].groupby(['Month'])['NEEtot'].std(ddof = 0).reset_index()
    FCMonthlyAnnualMeanGPP = PSIFFC_region[1][(PSIFFC_region[1].Year >= minyearmean)].groupby(['Month'])['GPPtot'].mean().reset_index()
    FCMonthlyAnnualStDevGPP = PSIFFC_region[1][(PSIFFC_region[1].Year >= minyearmean)].groupby(['Month'])['GPPtot'].std(ddof = 0).reset_index()
    FCMonthlyAnnualMeanTER = PSIFFC_region[2][(PSIFFC_region[2].Year >= minyearmean)].groupby(['Month'])['TERtot'].mean().reset_index()
    FCMonthlyAnnualStDevTER = PSIFFC_region[2][(PSIFFC_region[2].Year >= minyearmean)].groupby(['Month'])['TERtot'].std(ddof = 0).reset_index()
  
    
if TRENDY:
    PMonthSumTrendy = []
    PMonthMeansTrendy = []
    AllMMTrendyNBP = DateRef.copy()
    MeanCycle = pd.DataFrame(data = {'Month':np.array(range(1,13))})
    firstmodel = True
    for Model in TrendyModels:
        DataList = ['nbp','gpp','ra','rh','fFire','npp']# ['nbp','gpp','ra','rh']#['gpp','nbp','npp','ra','rh','fFire','fLuc']
        for Datatype in DataList:
            if os.path.isfile("/TRENDY/dataframes/MonthFrames/MonthMeansV2_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl"):
                MonthMeansTrendy = pd.read_pickle("/TRENDY/dataframes/MonthFrames/MonthMeansV2_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")
                MonthSumTrendy = pd.read_pickle("/TRENDY/dataframes/MonthFrames/MonthSumV2_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")
            else:
                print('get dataframe')
                if Model in ['CABLE-POP','DLEM'] and Datatype == 'nbp':
                    print(Model +' does not provide nbp, nbp will be calculated as npp - rh')
                    gdfTrendyNPP = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"npp_"+RegionName+str(Numm)+".pkl")
                    gdfTrendyRH = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"rh_"+RegionName+str(Numm)+".pkl")
                  
                    MonthSumTrendyNPP = getMonthSum(gdfTrendyNPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'npptot')
                    MonthMeansTrendyNPP = getMonthMeans(gdfTrendyNPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'npp',weighting = 'Area')
                    MonthSumTrendyRH = getMonthSum(gdfTrendyRH,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'rhtot')
                    MonthMeansTrendyRH = getMonthMeans(gdfTrendyRH,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'rh',weighting = 'Area')
                    
                    MonthMeansTrendy = pd.merge(MonthMeansTrendyNPP,MonthMeansTrendyRH,on=['Month','Year','Area'])
                    if not len(MonthMeansTrendy) == len(MonthMeansTrendyRH) or (not MonthMeansTrendy.Area.min() == MonthMeansTrendy.Area.max()):
                        raise Exception('check Region Area in TRENDY mean fluxes')
                    #convert units from kgC/m2s to TgCO2 in Australia
                    SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                    MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                    RegionArea = MonthMeansTrendy.Area.min()
                    
                    # value * SecPerMonth * g/kg * Tg/g * A_Region * gCO2/gC
                    rhTg = MonthMeansTrendy.rh*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12 * RegionArea
                    nppTg = MonthMeansTrendy.npp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12 * RegionArea 
                    MonthMeansTrendy.insert(loc = 1, column = 'rhTg', value=rhTg)
                    MonthMeansTrendy.insert(loc = 1, column = 'nppTg', value=nppTg)
                    MonthMeansTrendy.drop(['npp','rh'],axis = 1, inplace = True)
                    MonthMeansTrendy.rename(columns={'nppTg':'npp','rhTg':'rh'}, inplace = True)
                    MonthSumTrendy = pd.merge(MonthSumTrendyNPP,MonthSumTrendyRH,on=['Month','Year'])
                    MonthMeansTrendy.insert(loc = 1, value = MonthMeansTrendy.npp-MonthMeansTrendy.rh,column = 'nbp')
                    MonthSumTrendy.insert(loc = 1, value = MonthSumTrendy.npptot-MonthSumTrendy.rhtot,column = 'nbptot')
                elif Model in ['JULES-ES-1p0'] and Datatype == 'ra':
                    print(Model +' does not provide ra, calculated as gpp - npp as npp=gpp-ra')
                    gdfTrendyNPP = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"npp_"+RegionName+str(Numm)+".pkl")
                    gdfTrendyGPP = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"gpp_"+RegionName+str(Numm)+".pkl")
                  
                    MonthSumTrendyNPP = getMonthSum(gdfTrendyNPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'npptot')
                    MonthMeansTrendyNPP = getMonthMeans(gdfTrendyNPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'npp',weighting = 'Area')
                    MonthSumTrendyGPP = getMonthSum(gdfTrendyGPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'gpptot')
                    MonthMeansTrendyGPP = getMonthMeans(gdfTrendyGPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'gpp',weighting = 'Area')
                    MonthMeansTrendy = pd.merge(MonthMeansTrendyNPP,MonthMeansTrendyGPP,on=['Month','Year','Area'])
                    if not len(MonthMeansTrendy) == len(MonthMeansTrendyGPP) or (not MonthMeansTrendy.Area.min() == MonthMeansTrendy.Area.max()):
                        raise Exception('check Region Area in TRENDY mean fluxes')
                    RegionArea = MonthMeansTrendy.Area.min()
                    #convert units from kgC/m2s to TgCO2 in Australia
                    SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                    MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                    # value * SecPerMonth * g/kg * Tg/g * A_Australia * gCo2/gC 
                    gppTg = MonthMeansTrendy.gpp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12* RegionArea 
                    nppTg = MonthMeansTrendy.npp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12* RegionArea 
                    MonthMeansTrendy.insert(loc = 1, column = 'gppTg', value=gppTg)
                    MonthMeansTrendy.insert(loc = 1, column = 'nppTg', value=nppTg)
                    MonthMeansTrendy.drop(['npp','gpp'],axis = 1, inplace = True)
                    MonthMeansTrendy.rename(columns={'nppTg':'npp','gppTg':'gpp'}, inplace = True)
                    MonthSumTrendy = pd.merge(MonthSumTrendyNPP,MonthSumTrendyGPP,on=['Month','Year'])
                    MonthMeansTrendy.insert(loc = 1, value = MonthMeansTrendy.gpp-MonthMeansTrendy.npp,column = 'ra')
                    MonthSumTrendy.insert(loc = 1, value = MonthSumTrendy.gpptot-MonthSumTrendy.npptot,column = 'ratot')
                elif Model in ['VISIT'] and Datatype == 'npp':
                    print(Model +' does not provide npp, calculated as gpp - ra as npp=gpp-ra')
                    gdfTrendyra = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"ra_"+RegionName+str(Numm)+".pkl")
                    gdfTrendyGPP = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+"gpp_"+RegionName+str(Numm)+".pkl")
                  
                    MonthSumTrendyra = getMonthSum(gdfTrendyra,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'ratot')
                    MonthMeansTrendyra = getMonthMeans(gdfTrendyra,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'ra',weighting = 'Area')
                    MonthSumTrendyGPP = getMonthSum(gdfTrendyGPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'gpptot')
                    MonthMeansTrendyGPP = getMonthMeans(gdfTrendyGPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'gpp',weighting = 'Area')
                    MonthMeansTrendy = pd.merge(MonthMeansTrendyra,MonthMeansTrendyGPP,on=['Month','Year','Area'])
                    if not len(MonthMeansTrendy) == len(MonthMeansTrendyGPP) or (not MonthMeansTrendy.Area.min() == MonthMeansTrendy.Area.max()):
                        raise Exception('check Region Area in TRENDY mean fluxes')
                    RegionArea = MonthMeansTrendy.Area.min()
                    #convert units from kgC/m2s to TgCO2 in Australia
                    SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                    MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                    # value * SecPerMonth * g/kg * Tg/g * A_Australia * gCo2/gC 
                    gppTg = MonthMeansTrendy.gpp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12* RegionArea 
                    raTg = MonthMeansTrendy.ra*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12* RegionArea 
                    MonthMeansTrendy.insert(loc = 1, column = 'gppTg', value=gppTg)
                    MonthMeansTrendy.insert(loc = 1, column = 'raTg', value=raTg)
                    MonthMeansTrendy.drop(['ra','gpp'],axis = 1, inplace = True)
                    MonthMeansTrendy.rename(columns={'raTg':'ra','gppTg':'gpp'}, inplace = True)
                    MonthSumTrendy = pd.merge(MonthSumTrendyra,MonthSumTrendyGPP,on=['Month','Year'])
                    MonthMeansTrendy.insert(loc = 1, value = MonthMeansTrendy.gpp-MonthMeansTrendy.ra,column = 'npp')
                    MonthSumTrendy.insert(loc = 1, value = MonthSumTrendy.gpptot-MonthSumTrendy.ratot,column = 'npptot')
                else:
                    if Model in ['CABLE-POP','DLEM','IBIS','ISAM','JULES-ES-1p0','OCN','ORCHIDEEv3','YIBs'] and Datatype == 'fFire':#models without fire emissions
                        # empty fire variable
                        #print(len(MonthSumTrendy))
                        MonthSumTrendy = DateRef.copy(deep = True)
                        MonthMeansTrendy = DateRef.copy(deep = True)
                        #print(len(MonthSumTrendy))
                        nanarray = np.empty(len(MonthSumTrendy.MonthDate))
                        nanarray = np.nan
                        onesarray = np.ones(len(MonthSumTrendy.MonthDate))
                        MonthSumTrendy.insert(loc = 1, column = 'fFire', value=nanarray)
                        MonthMeansTrendy.insert(loc = 1, column = 'fFire', value=nanarray)
                        MonthSumTrendy.insert(loc = 1, column = 'Area', value=onesarray)
                        MonthMeansTrendy.insert(loc = 1, column = 'Area', value=onesarray)
                        MonthSumTrendy = MonthSumTrendy[['fFire','Year','Month','Area']]
                        MonthMeansTrendy = MonthMeansTrendy[['fFire','Year','Month','Area']]
                    else:
                        gdfTrendy = pd.read_pickle("/TRENDY/dataframes/GDF2"+Model+Datatype+"_"+RegionName+str(Numm)+".pkl")
                    
                        MonthSumTrendy = getMonthSum(gdfTrendy,year_min,month_min,day_min,
                                                    year_max,month_max,day_max,
                                                    Long_min,Long_max,Lat_min,Lat_max,
                                                    Datatype+'tot')
                        MonthMeansTrendy = getMonthMeans(gdfTrendy,year_min,month_min,day_min,
                                                        year_max,month_max,day_max,
                                                        Long_min,Long_max,Lat_min,Lat_max,
                                                        Datatype,weighting = 'Area')
                        
                    if Datatype is not 'mrso':
                        #print('convert units from mugCo2/m2s to TgCO2 in Australia')
                        #print('convert units from kgC/m2s to TgCO2/month in Australia')
                        SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                        MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                        if not MonthMeansTrendy.Area.min() == MonthMeansTrendy.Area.max():
                            print(MonthMeansTrendy.Area.min())
                            print(MonthMeansTrendy.Area.max())
                            raise Exception('check Region Area in TRENDY mean fluxes')
                        RegionArea = MonthMeansTrendy.Area.min()
                        # value * SecPerMonth * g/kg * Tg/g * A_Region *gCO2/gC
                        varTg = MonthMeansTrendy[Datatype]*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * 44/12* RegionArea 
                        MonthMeansTrendy.insert(loc = 1, column = 'varTg', value=varTg)
                        MonthMeansTrendy.drop([Datatype],axis = 1, inplace = True)
                        MonthMeansTrendy.rename(columns={'varTg':Datatype}, inplace = True)
                Monthdate = MonthSumTrendy.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
                MonthSumTrendy.insert(loc=1,column='MonthDate',value=Monthdate)
                Monthdate = MonthMeansTrendy.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
                MonthMeansTrendy.insert(loc=1,column='MonthDate',value=Monthdate)
                #MonthMeansTrendy.to_pickle("/TRENDY/dataframes/MonthFrames/MonthMeans_"+Model+"_"+Datatype+RegionName+str(Numm)+"_newForSuppPaper.pkl")
                MonthMeansTrendy.to_pickle("/TRENDY/dataframes/MonthFrames/MonthMeansV2_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")
                MonthSumTrendy.to_pickle("/TRENDY/dataframes/MonthFrames/MonthSumV2_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")                
                print('got Trendy ' + Datatype)
            MeanCycle = MeanCycle.merge(MonthMeansTrendy[(MonthMeansTrendy.Year >= year_min)&(MonthMeansTrendy.Year <= year_max)].groupby(['Month'])[Datatype].mean(), on = 'Month')
            MeanCycle = MeanCycle.rename(columns = {Datatype: Model+Datatype})
            PMonthSumTrendy.append(pd.merge(DateRef, MonthSumTrendy, on=['MonthDate'], how = 'left'))   
            PMonthMeansTrendy.append(pd.merge(DateRef, MonthMeansTrendy, on=['MonthDate'], how = 'left'))   
            AllMMTrendyNBP = pd.merge(AllMMTrendyNBP,MonthMeansTrendy[[Datatype,'MonthDate']],on='MonthDate', how = 'left')
            AllMMTrendyNBP = AllMMTrendyNBP.rename(columns={Datatype:Model+Datatype})
        if True: #set false if using precipitation
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'ra']+MeanCycle[Model+'rh'],column = Model+'rarh')
            maxFlux = max([MeanCycle[Model+'rarh'].max(),-MeanCycle[Model+'gpp'].max()])
            MeanCycle.insert(loc=1,
                            value= MeanCycle[Model+'rarh']/maxFlux,column = Model+'rarhNorm')
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'gpp']/maxFlux,column = Model+'gppNorm')
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'nbp']/maxFlux,column = Model+'nbpNorm')

    if len(TrendyModels) >= 16:        
            #get mean and stdev for all models
        meanModel = AllMMTrendyNBP[BadModelNames('nbp',Numm)+GoodModelsNames('nbp',Numm)].mean(axis=1)
        stDevModel = AllMMTrendyNBP[BadModelNames('nbp',Numm)+GoodModelsNames('nbp',Numm)].std(axis=1,ddof = 0)
        minModel = AllMMTrendyNBP[BadModelNames('nbp',Numm)+GoodModelsNames('nbp',Numm)].min(axis=1)
        maxModel = AllMMTrendyNBP[BadModelNames('nbp',Numm)+GoodModelsNames('nbp',Numm)].max(axis=1)
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBP', value = minModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBP', value = maxModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBP', value = meanModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBP', value = stDevModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp',Numm)].min(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp',Numm)].max(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp',Numm)].std(axis=1,ddof = 0))
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp',Numm)].min(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp',Numm)].max(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp',Numm)].std(axis=1,ddof = 0))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('nbp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPlateModel', value = AllMMTrendyNBP[LateRespNames('nbp',Numm)].mean(axis=1))
        
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppbadModel', value = AllMMTrendyNBP[BadModelNames('gpp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppgoodModel', value = AllMMTrendyNBP[GoodModelsNames('gpp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('gpp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpplateModel', value = AllMMTrendyNBP[LateRespNames('gpp',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhrabadModel', value = AllMMTrendyNBP[BadModelNames('ra',Numm)].mean(axis=1)+AllMMTrendyNBP[BadModelNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhragoodModel', value = AllMMTrendyNBP[GoodModelsNames('ra',Numm)].mean(axis=1)+AllMMTrendyNBP[GoodModelsNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhraearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('ra',Numm)].mean(axis=1)+AllMMTrendyNBP[EarlyGPPModelNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhralateModel', value = AllMMTrendyNBP[LateRespNames('ra',Numm)].mean(axis=1)+AllMMTrendyNBP[LateRespNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhbadModel', value = AllMMTrendyNBP[BadModelNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhgoodModel', value = AllMMTrendyNBP[GoodModelsNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhlateModel', value = AllMMTrendyNBP[LateRespNames('rh',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrabadModel', value = AllMMTrendyNBP[BadModelNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanragoodModel', value = AllMMTrendyNBP[GoodModelsNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanraearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanralateModel', value = AllMMTrendyNBP[LateRespNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpp-rabadModel', value = AllMMTrendyNBP[BadModelNames('gpp',Numm)].mean(axis=1)-AllMMTrendyNBP[BadModelNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpp-ragoodModel', value = AllMMTrendyNBP[GoodModelsNames('gpp',Numm)].mean(axis=1)-AllMMTrendyNBP[GoodModelsNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpp-raearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('gpp',Numm)].mean(axis=1)-AllMMTrendyNBP[EarlyGPPModelNames('ra',Numm)].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpp-ralateModel', value = AllMMTrendyNBP[LateRespNames('gpp',Numm)].mean(axis=1)-AllMMTrendyNBP[LateRespNames('ra',Numm)].mean(axis=1))
        
        #get mean amplitude
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPgoodModelSwitch',value = -1*AllMMTrendyNBP.meanNBPgoodModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPbadModelSwitch',value = -1*AllMMTrendyNBP.meanNBPbadModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPSwitch',value = -1*AllMMTrendyNBP.meanNBP)
        AllMMTrendyNBP.to_pickle("/TRENDY/dataframes/MonthFrames/AllMMTrendyNBP"+RegionName+str(Numm)+"_newforSuppPaper.pkl")   
        goodTRENDYmeanAmplitude, goodTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','Year',2009,2018)[0:2]
        goodTRENDYmeanAmplitude_no09, goodTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','Year',2010,2018)[0:2]
        goodTRENDYmeanAmplitudejj, goodTRENDYmeanAmplitudeStDevjj,goodTRENDYmeanAmplitudeMin,goodTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','JJYear',2009,2017)
        badTRENDYmeanAmplitude, badTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','Year',2009,2018)[0:2]
        badTRENDYmeanAmplitude_no09, badTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','Year',2010,2018)[0:2]
        badTRENDYmeanAmplitudejj, badTRENDYmeanAmplitudeStDevjj,badTRENDYmeanAmplitudeMin,badTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','JJYear',2009,2017)
        allTRENDYmeanAmplitude, allTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','Year',2009,2018)[0:2]
        allTRENDYmeanAmplitude_no09, allTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','Year',2010,2018)[0:2]
        allTRENDYmeanAmplitudejj, allTRENDYmeanAmplitudeStDevjj,allTRENDYmeanAmplitudeMin,allTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','JJYear',2009,2017)
        
        TRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBP'].sum().mean()
        TRENDYannualStDev = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBP'].sum().std(ddof = 0)
        TRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp',Numm)+BadModelNames('nbp',Numm)].sum().mean(axis=0).min()
        TRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp',Numm)+BadModelNames('nbp',Numm)].sum().mean(axis=0).max()
        #print('mean: '+str(TRENDYannualMean))
        #print('std: '+str(TRENDYannualStDev))
        goodTRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPgoodModel'].sum().mean()
        badTRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPbadModel'].sum().mean()
        goodTRENDYannualStd = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPgoodModel'].sum().std(ddof = 0)
        badTRENDYannualStd = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPbadModel'].sum().std(ddof = 0)
        goodTRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp',Numm)].sum().mean(axis=0).min()
        #print(AllMMTrendyNBP.groupby(['Year'])[GoodModelsNames('nbp')].sum())
        #print(AllMMTrendyNBP.groupby(['Year'])[GoodModelsNames('nbp')].sum().mean(axis = 0))
        badTRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[BadModelNames('nbp',Numm)].sum().mean(axis=0).min()
        goodTRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp',Numm)].sum().mean(axis=0).max()
        badTRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[BadModelNames('nbp',Numm)].sum().mean(axis=0).max()
        
        MeanGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanAll = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevAll = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        
    if len(TrendyModels) >= 16:
        for Datatype in DataList+['rarh','rarhNorm','gppNorm','nbpNorm']:
            GoodModel = MeanCycle[GoodModelsNames(Datatype,Numm)].copy()
            AllModel = MeanCycle[GoodModelsNames(Datatype,Numm)+BadModelNames(Datatype,Numm)].copy()
            EarlyGPPModel = MeanCycle[EarlyGPPModelNames(Datatype,Numm)].copy()
            LateResp = MeanCycle[LateRespNames(Datatype,Numm)].copy()
            BadModel = MeanCycle[BadModelNames(Datatype,Numm)].copy()
            MeanGood.insert(loc=1,
                            value=GoodModel.mean(axis=1), column = Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.mean(axis=1), column = Datatype)
            MeanAll.insert(loc=1,
                            value=AllModel.mean(axis=1), column = Datatype)
            MeanGood.insert(loc=1,
                            value=GoodModel.min(axis=1), column = 'min_'+Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.min(axis=1), column = 'min_'+Datatype)
            MeanGood.insert(loc=1,
                            value=GoodModel.max(axis=1), column = 'max_'+Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.max(axis=1), column = 'max_'+Datatype)
            MeanEarlyGPP.insert(loc=1, value=EarlyGPPModel.mean(axis=1),
                                column = Datatype)
            MeanLateResp.insert(loc=1, value=LateResp.mean(axis=1),
                                column = Datatype)
            StErrorGood.insert(loc=1,value=GoodModel.sem(axis=1,ddof = 0), 
                               column = Datatype)
            StErrorBad.insert(loc=1,value=BadModel.sem(axis=1,ddof = 0), 
                               column = Datatype)
            StErrorEarlyGPP.insert(loc=1, value=EarlyGPPModel.sem(axis=1,ddof = 0),
                                column = Datatype)
            StErrorLateResp.insert(loc=1, value=LateResp.sem(axis=1,ddof = 0),
                                column = Datatype)
            StDevGood.insert(loc=1,value=GoodModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevBad.insert(loc=1,value=BadModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevAll.insert(loc=1,value=AllModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevEarlyGPP.insert(loc=1, value=EarlyGPPModel.std(axis=1,ddof = 0),
                                column = Datatype)
            StDevLateResp.insert(loc=1, value=LateResp.std(axis=1,ddof = 0),
                                column = Datatype)
            AllModel.insert(1,column = 'Month', value = range(1,13))
            if 'AllModelall' in locals():
                AllModelall = pd.merge(AllModelall,AllModel, on=['Month'])    
            else:
                AllModelall = AllModel.copy(deep = True)
    

if GOMESIFv2:
    GOMEvar = 'SIF_740'
    if False:
        pass
    else:
        if not os.path.isfile("/SIF/GOME2_v2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl"):
            print('prepare GOME2 v2')
            CreateDataFrameSIFGOMEv2(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
        print('prepare Monthmeans')
        gdfGOMESIFv2 = pd.read_pickle("/SIF/GOME2_v2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl")

        Month_means_GOMESIFv2 = getMonthMeans(gdfGOMESIFv2,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,GOMEvar)
        MonthD = Month_means_GOMESIFv2.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
        Month_means_GOMESIFv2.insert(loc = 1, column = 'MonthDate', value = MonthD)
        Month_means_GOMESIFv2.to_pickle("/SIF/GOME2_v2/DataFrames/monthFrames/MonthMeans_GOMESIFv2_"+RegionName+ str(Numm)+".pkl")
    
        gdfGOMESIFv2grid = gdfGOMESIFv2.groupby(['Lat','Long','Year','Month'])[GOMEvar].mean().reset_index()
        MonthDAllg = gdfGOMESIFv2grid.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
        gdfGOMESIFv2grid.insert(loc = 1, column = 'Date', value = MonthDAllg)
        Month_means_GOMESIFv2grid = getMonthMeans(gdfGOMESIFv2grid,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,GOMEvar)
        MonthDg = Month_means_GOMESIFv2grid.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
        Month_means_GOMESIFv2grid.insert(loc = 1, column = 'MonthDate', value = MonthDg)
        Month_means_GOMESIFv2grid.to_pickle("/SIF/GOME2_v2/DataFrames/monthFrames/MonthMeans_GOMESIFv2grid_"+RegionName+ str(Numm)+".pkl")
    
        gdfGOMESIFv2days = gdfGOMESIFv2.groupby(['Year','Month','Day','Date'])[GOMEvar].mean().reset_index()
        LatsSIF = gdfGOMESIFv2days.apply(lambda x: Lat_max,axis=1)
        gdfGOMESIFv2days.insert(loc = 1, column = 'Lat', value = LatsSIF)
        LongsSIF = gdfGOMESIFv2days.apply(lambda x: Long_max,axis=1)
        gdfGOMESIFv2days.insert(loc = 1, column = 'Long', value = LongsSIF)
        Month_means_GOMESIFv2days = getMonthMeans(gdfGOMESIFv2days,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,GOMEvar)
        MonthDg = Month_means_GOMESIFv2days.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
        Month_means_GOMESIFv2days.insert(loc = 1, column = 'MonthDate', value = MonthDg)
        Month_means_GOMESIFv2days.to_pickle("/SIF/GOME2_v2/DataFrames/monthFrames/MonthMeans_GOMESIFv2days_"+RegionName+ str(Numm)+".pkl")
    

    PMonth_means_GOMESIFv2 = pd.merge(DateRef, Month_means_GOMESIFv2, on=['MonthDate','Year'], how = 'left')
    PMonth_means_GOMESIFv2grid = pd.merge(DateRef, Month_means_GOMESIFv2grid, on=['MonthDate','Year'], how = 'left')
    PMonth_means_GOMESIFv2days = pd.merge(DateRef, Month_means_GOMESIFv2days, on=['MonthDate','Year'], how = 'left')


        
#calculate mean seasonal cycle and key numbers
if TM5Flux1x1:#1x1:#mean of the satellite based TM5
    #1x1 flux
    
    allSatModel = TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_ACOS+IS')[0][0]][['MonthDate','JJYear','CO2_flux_nee_total','CO2_flux_fire_total']]
    allSatModel.insert(loc = 1, column = 'TM5ACOSIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
    allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
    allSatModel = pd.merge(allSatModel,TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_RemoTeC_2.4.0+IS')[0][0]][['MonthDate','CO2_flux_nee_total','CO2_flux_fire_total']])
    allSatModel.insert(loc = 1, column = 'TM5RTIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
    allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
    
    
    ParamList = ['TM5RTIS_nbp','TM5ACOSIS_nbp']#,'TM5RT238IS_nbp']
    if TM5Flux1x1 and WithOCO2:
        print(np.where(np.array(TM51x1dataset) == 'flux_1x1_LNLGIS')[0][0])
        print(len(TM51x1dataset))
        allSatModel = pd.merge(allSatModel,TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_LNLGIS')[0][0]][['MonthDate','CO2_flux_nbe_total']])
        allSatModel.insert(loc = 1, column = 'TM5OCO2IS_nbp', value = allSatModel.CO2_flux_nbe_total*10**(-12))
        print(allSatModel.TM5OCO2IS_nbp.max())
        print(allSatModel.TM5RTIS_nbp.max())
        print(allSatModel.TM5ACOSIS_nbp.max())
        allSatModel = allSatModel.drop(columns = ['CO2_flux_nbe_total'])
        ParamList = ['TM5RTIS_nbp','TM5ACOSIS_nbp','TM5OCO2IS_nbp']
    allSatModel.insert(loc = 1, column = 'mean_nbp', value = allSatModel[ParamList].mean(axis = 1))
    allSatModel.insert(loc = 1, column = 'std_nbp', value = allSatModel[ParamList].std(axis = 1,ddof = 0))
    allSatModel.insert(loc = 1, column = 'min_nbp', value = allSatModel[ParamList].min(axis = 1))
    allSatModel.insert(loc = 1, column = 'max_nbp', value = allSatModel[ParamList].max(axis = 1))
    Monthlist = allSatModel.apply(lambda x: x.MonthDate.month,axis=1)
    allSatModel.insert(loc = 1, column = 'Month', value = Monthlist)
    Yearlist = allSatModel.apply(lambda x: x.MonthDate.year,axis=1)
    allSatModel.insert(loc = 1, column = 'Year', value = Yearlist) 
    #allSatModel.to_pickle("/TM5Inversion/dataframes_flux/MonthFrames/MonthMeansAllSatModels_newForSuppPaper_"+RegionName+str(Numm)+".pkl")          
    AMeanSatModel = allSatModel.groupby(['Month'])['mean_nbp'].mean().reset_index()
    
    #get mean amplitude
    allSatModelYmeanAmplitude, allSatModelmeanAmplitudeStDev,allSatModelmeanAmplitudeMin,allSatModelmeanAmplitudeMax = getMeanAmplitude(allSatModel,'mean_nbp','Year',2009,2018)
    allSatModelYmeanAmplitude_no09, allSatModelmeanAmplitudeStDev_no09 = getMeanAmplitude(allSatModel,'mean_nbp','Year',2010,2018)[0:2]
    allSatModelYmeanAmplitudejj, allSatModelmeanAmplitudeStDevjj,allSatModelmeanAmplitudeMinjj,allSatModelmeanAmplitudeMaxjj = getMeanAmplitude(allSatModel,'mean_nbp','JJYear',2009,2017)
    
    allSatModelannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().mean()
    allSatModelannualStDev = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().std(ddof = 0)
    allSatModelMeanSais = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].mean().reset_index()
    allSatModelMeanSaisStDev = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].std(ddof = 0).reset_index()
    
    RTannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5RTIS_nbp'].sum().mean()
    ACOSannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5ACOSIS_nbp'].sum().mean()
    if WithOCO2:
        OCO2annualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5OCO2IS_nbp'].sum().mean()
    
    if WithOCO2:
        satmeans = [ACOSannualMean,RTannualMean,OCO2annualMean]
    else:
        satmeans = [ACOSannualMean,RTannualMean]

if CAMS_b and CT_flux2022 and TM5Flux1x1:#mean of the inverse models
    
    if TM5Flux1x1:
        allModel = TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_IS')[0][0]][['MonthDate','JJYear','CO2_flux_nee_total','CO2_flux_fire_total']].copy()
        allModel.insert(loc = 1, column = 'TM5_nbp', value = (allModel.CO2_flux_nee_total+allModel.CO2_flux_fire_total)*10**(-12))
        allModel = allModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
    allModel = pd.merge(allModel,PCTF22_region[['MonthDate','Biotot','firetot']])
    allModel.insert(loc = 1, column = 'CT_nbp', value = allModel.Biotot*10**-12+allModel.firetot*10**-12)
    allModel = allModel.drop(columns = ['Biotot','firetot'])
    allModel = pd.merge(allModel,PCAMS_region_Sur[['MonthDate','flux_apos_bio_tot']])
    allModel.insert(loc = 1, column = 'CAMSsur_nbp', value = allModel.flux_apos_bio_tot*10**-12)
    allModel = allModel.drop(columns = ['flux_apos_bio_tot'])
    allModel.insert(loc = 1, column = 'mean_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].mean(axis = 1))
    allModel.insert(loc = 1, column = 'std_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].std(axis = 1,ddof = 0))
    allModel.insert(loc = 1, column = 'min_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].min(axis = 1))
    allModel.insert(loc = 1, column = 'max_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].max(axis = 1))
    Monthlist = allModel.apply(lambda x: x.MonthDate.month,axis=1)
    allModel.insert(loc = 1, column = 'Month', value = Monthlist)
    Yearlist = allModel.apply(lambda x: x.MonthDate.year,axis=1)
    allModel.insert(loc = 1, column = 'Year', value = Yearlist)     
    
        
    #get mean Amplitude
    allModelYmeanAmplitude, allModelmeanAmplitudeStDev = getMeanAmplitude(allModel,'mean_nbp','Year',2009,2018)[0:2]
    allModelYmeanAmplitude_no09, allModelmeanAmplitudeStDev_no09 = getMeanAmplitude(allModel,'mean_nbp','Year',2010,2018)[0:2]
    allModelYmeanAmplitudejj, allModelmeanAmplitudeStDevjj,allModelmeanAmplitudeMin,allModelmeanAmplitudeMax = getMeanAmplitude(allModel,'mean_nbp','JJYear',2009,2017)
    
    allModelannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().mean()
    allModelannualStDev = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().std(ddof = 0)
    allModelMeanSais = allModel[(allModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].mean().reset_index()
    allModelMeanSaisStDev = allModel[(allModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].std(ddof = 0).reset_index()
    
    CTannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['CT_nbp'].sum().mean()
    CAMSannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['CAMSsur_nbp'].sum().mean()
    TM5annualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['TM5_nbp'].sum().mean()
    

###################
##Africa Paper Plots
###################
# Figure 1 is study area plot

# Figure 2 is created with Plot_TimeSeries.py

# Figure 3 
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.1*cm),
                            gridspec_kw={'width_ratios': [4, 1]})

    #get TRENDY std
    StdCycnbp = AllMMTrendyNBP.groupby('Month')['meanNBP'].std(ddof=0).reset_index()
    StdCycgoodnbp = AllMMTrendyNBP.groupby('Month')['meanNBPgoodModel'].std(ddof=0).reset_index()
    StdCycbadnbp = AllMMTrendyNBP.groupby('Month')['meanNBPbadModel'].std(ddof=0).reset_index()
        

    # Panel a)
    ax2[0].text(0.03, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10, weight='bold')
    #Shading
    ## stdDev for TRENDY
    
    ax2[0].fill_between(AllMMTrendyNBP.MonthDate,
                        -(AllMMTrendyNBP.meanNBP-AllMMTrendyNBP.stdNBP)*factor,
                        -(AllMMTrendyNBP.meanNBP+AllMMTrendyNBP.stdNBP)*factor,
                        color= 'grey', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    
    ##Min Max for Sat and Models
    '''
    ax2[0].fill_between(allModel.MonthDate,allModel.max_nbp*factor,
                          allModel.min_nbp*factor
                          ,color= 'blue', alpha=0.3,zorder = 2)
    '''
    ax2[0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    #ax2[0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPgoodModel-AllMMTrendyNBP.stdNBPgoodModel)*factor,-(AllMMTrendyNBP.meanNBPgoodModel+AllMMTrendyNBP.stdNBPgoodModel)*factor,color= 'black', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    
    #monthly means
    ax2[0].plot([allModel.MonthDate.values[0],allModel.MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    #ax2[0].plot(allModel.MonthDate,allModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label='In-situ-only inversions',zorder = 4)
    ax2[0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/IS+GOSAT',zorder = 4)
    ##ax2[0].plot(PSIFFC_region[0].MonthDate,(PSIFFC_region[0].NEEtot/10**12*(44/12)+PGFEDregion.total_emission)*factor,ls='-',linewidth = lw,marker= '',color= 'gold',label='FLUXCOM+GFED',zorder = 4)
    ax2[0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBP*factor,ls='-',linewidth = lw,marker= '',color= 'dimgrey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ##ax2[0].plot(PSIFFC_region[0].MonthDate,(PGFEDregion.total_emission)*factor,ls='-',linewidth = lw,marker= '',color= 'orange',label=r'$\rm GFED~fire~CO_2$',zorder = 4)
    #ax2[0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBPgoodModel*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{selection}$',zorder = 4)
    
    # Panel b)
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]

    ax2[1].text(0.15, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10, weight='bold')
    ax2[1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1].plot(range(1,13),-MeanAll.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'dimgrey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ##ax2[1].plot(range(1,13),(FCMonthlyAnnualMean.nbp.iloc[indexrange])*factor,ls='-',linewidth = lw,marker= '',color= 'gold',label='FLUXCOM+GFED',zorder = 4)
    #ax2[1].plot(range(1,13),allModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm Inverse~model_{in-situ}$',zorder = 4)
    ax2[1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ##ax2[1].plot(range(1,13),GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'orange',label=r'$\rm GFED~fire~CO_{2}$',zorder = 4)    
    # SHADING
    ##ax2[1].fill_between(range(1,13),(allModelMeanSais.mean_nbp.iloc[indexrange]-allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allModelMeanSais.mean_nbp.iloc[indexrange]+allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'blue', alpha=0.3,zorder = 2)
    ax2[1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    #ax2[1].fill_between(range(1,13),-(MeanAll.nbp.iloc[indexrange]-StDevAll.nbp.iloc[indexrange])*factor,-(MeanAll.nbp.iloc[indexrange]+StDevAll.nbp.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[1].fill_between(range(1,13),-(MeanAll.nbp.iloc[indexrange]-StdCycnbp.meanNBP.iloc[indexrange])*factor,-(MeanAll.nbp.iloc[indexrange]+StdCycnbp.meanNBP.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    #ax2[1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    ##ax2[1].fill_between(range(1,13),(FCMonthlyAnnualMean.nbp.iloc[indexrange]-FCMonthlyAnnualStDev.nbp.iloc[indexrange])*factor,(FCMonthlyAnnualMean.nbp.iloc[indexrange]+FCMonthlyAnnualStDev.nbp.iloc[indexrange])*factor,color= 'gold', alpha=0.3,zorder = 2)
    ##ax2[1].fill_between(range(1,13),(GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]-GFEDMonthlyAnnualStd.total_emission.iloc[indexrange])*factor,(GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]+GFEDMonthlyAnnualStd.total_emission.iloc[indexrange])*factor,color= 'orange', alpha=0.3,zorder = 2)
    #ax2[1].plot(range(1,13),-MeanGood.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{all}$',zorder = 4)
    
    # SETTINGS
    
    #Panel a)
    
    ax2[0].legend(fontsize=6,loc = 4,ncol = 3)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    ax2[0].set_ylim(-350,350)
    #ax2[0].set_ylim(-280,280)
    #ax2[0].set_xticklabels(ax2[1,0].get_xticks(),fontsize=7)
    #ax2[0].set_yticklabels(ax2[1,0].get_yticks(),fontsize=7)
    ax2[0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[0].set_axisbelow(True)

    #Panel b)  
          
    secaxB = ax2[0].twinx()     
    secaxB.set_ylim(-350,350)
    #secaxB.set_ylim(-280,280)
    secaxB.set_yticklabels([])
    secaxB.tick_params(axis = 'y',direction = 'in')
    ax2[1].set_xticks(ticks = [3,6,9,12])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    #ax2[1].set_ylim(-280,280)
    ax2[1].set_ylim(-350,350)
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    print('done')
    plt.savefig("/Results/Plots/Africa/Paper2024/Fig3_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')

# Figure 4 TM5 prior with MSC
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.1*cm),
                            gridspec_kw={'width_ratios': [4, 1]})

    #get TM5 MSC std
    TM51x1_regionl[1].insert(0,column = 'NBP', value = TM51x1_regionl[1].CO2_flux_nee_total+TM51x1_regionl[1].CO2_flux_fire_total)
    TM51x1_regionl[3].insert(0,column = 'NBP', value = TM51x1_regionl[3].CO2_flux_nee_total+TM51x1_regionl[3].CO2_flux_fire_total)
    TM51x1_regionl[1] = pd.merge(DateRef, TM51x1_regionl[1], on=['MonthDate'], how = 'left')
    TM51x1_regionl[3] = pd.merge(DateRef, TM51x1_regionl[3], on=['MonthDate'], how = 'left')
    StdCycTM5IS = TM51x1_regionl[1].groupby('Month')['NBP'].std(ddof=0).reset_index()
    StdCycTM5prior = TM51x1_regionl[3].groupby('Month')['NBP'].std(ddof=0).reset_index()
    MeanCycTM5IS = TM51x1_regionl[1].groupby('Month')['NBP'].mean().reset_index()
    MeanCycTM5prior = TM51x1_regionl[3].groupby('Month')['NBP'].mean().reset_index()

    # Panel a)
    ax2[0].text(0.03, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10, weight='bold')
    ax2[0].plot([DateRef.MonthDate[0],DateRef.MonthDate[len(DateRef.MonthDate)-1]], [0,0], ls='-',linewidth = lw,marker= '',color= 'lightgrey',zorder = 4)
    
    ax2[0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    #monthly means
    ax2[0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'TM5-4DVar/GOSAT+IS',zorder = 4)
    ax2[0].plot(TM51x1_regionl[1].MonthDate, (TM51x1_regionl[1].CO2_flux_nee_total+TM51x1_regionl[1].CO2_flux_fire_total)*factor*10**(-12), ls='-',linewidth = lw,marker= '',color= 'grey',label=r'TM5-4DVar/IS',zorder = 4)
    ax2[0].plot(TM51x1_regionl[3].MonthDate, (TM51x1_regionl[3].CO2_flux_nee_total+TM51x1_regionl[3].CO2_flux_fire_total)*factor*10**(-12), ls=':',linewidth = lw,marker= '',color= 'grey',label='TM5-4DVa'+r'$\rm r_{prior}$',zorder = 4)
     
    # Panel b)
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]

    ax2[1].text(0.15, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10, weight='bold')
    ax2[1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2[1].plot(range(1,13),MeanCycTM5IS.NBP.iloc[indexrange]*factor*10**(-12),ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 4)    
    ax2[1].plot(range(1,13),MeanCycTM5prior.NBP.iloc[indexrange]*factor*10**(-12),ls=':',linewidth = lw,marker= '',color= 'grey',zorder = 4)    
    # SHADING
    ax2[1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    ax2[1].fill_between(range(1,13),(MeanCycTM5IS.NBP.iloc[indexrange]-StdCycTM5IS.NBP.iloc[indexrange])*factor*10**(-12),(MeanCycTM5IS.NBP.iloc[indexrange]+StdCycTM5IS.NBP.iloc[indexrange])*factor*10**(-12),color= 'black', alpha=0.1,zorder = 2)
    ax2[1].fill_between(range(1,13),(MeanCycTM5prior.NBP.iloc[indexrange]-StdCycTM5prior.NBP.iloc[indexrange])*factor*10**(-12),(MeanCycTM5prior.NBP.iloc[indexrange]+StdCycTM5prior.NBP.iloc[indexrange])*factor*10**(-12),color= 'black', alpha=0.2,zorder = 2)
    
    # SETTINGS
    
    #Panel a)
    
    ax2[0].legend(fontsize=6,loc = 4,ncol = 3)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    ax2[0].set_ylim(-350,350)
    #ax2[0].set_ylim(-280,280)
    #ax2[0].set_xticklabels(ax2[1,0].get_xticks(),fontsize=7)
    #ax2[0].set_yticklabels(ax2[1,0].get_yticks(),fontsize=7)
    ax2[0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[0].set_axisbelow(True)

    #Panel b)  
          
    secaxB = ax2[0].twinx()     
    secaxB.set_ylim(-350,350)
    #secaxB.set_ylim(-280,280)
    secaxB.set_yticklabels([])
    secaxB.tick_params(axis = 'y',direction = 'in')
    ax2[1].set_xticks(ticks = [3,6,9,12])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    #ax2[1].set_ylim(-280,280)
    ax2[1].set_ylim(-350,350)
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 

    plt.savefig("/Results/Plots/Africa/Paper2024/Figure4__wMSC"+RegionName+str(Numm)+".png", dpi=250,bbox_inches='tight')#,format = 'pdf')



# Figure 5 MIP Ensemble
if False:
    #Panel B
    Tm5OCOs = TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_LNLGIS')[0][0]]
    Tm5OCOs.insert(1,column ='Month',value = Tm5OCOs.apply(lambda x: x.MonthDate.month,axis = 1))
    Tm5OCOmean = Tm5OCOs[(Tm5OCOs.MonthDate >= datetime.date(minyearmean,1,1))].groupby(['Month'])['CO2_flux_nbe_total'].mean().reset_index()
    Tm5OCOmeanStDev = Tm5OCOs[(Tm5OCOs.MonthDate >= datetime.date(minyearmean,1,1))].groupby(['Month'])['CO2_flux_nbe_total'].std(ddof = 0).reset_index()

    #all Models without LoFi = provided model ensemble mean
    AllwoLofinames = ['UT','WOMBAT','OU','CSU','NIES','CAMS','CMS-Flux','JHU','CT','COLA','TM5-4DVAR','Baker','Ames']
    AllwoLofi_CTnames = ['UT','WOMBAT','OU','CSU','NIES','CAMS','CMS-Flux','JHU','COLA','TM5-4DVAR','Baker','Ames']
    
    #without LoFI is the right one as LoFI also not included in provided mean ensemble data
    AllMIPwoLofi = getEnsembleFluxMIP(PMIPOCOModel, PMIPpriorModel, MIPModelNames, MIPModelNamesPrior,AllwoLofinames)
 
    CAMSTM5BAKER = getEnsembleFluxMIP(PMIPOCOModel, PMIPpriorModel, MIPModelNames, MIPModelNamesPrior,['TM5-4DVAR','Baker','CAMS'])

    AllMIPwoLofiIS = getEnsembleFluxMIP(PMIPISModel, PMIPpriorModel, MIPModelNames, MIPModelNamesPrior,AllwoLofinames)
    
    AllMIPwoLofiPrior = getEnsembleFluxMIP(PMIPpriorModel, PMIPpriorModel, MIPModelNamesPrior, MIPModelNamesPrior,AllwoLofinames)
    AllMIPwoLofiCTPrior = getEnsembleFluxMIP(PMIPpriorModel, PMIPpriorModel, MIPModelNamesPrior, MIPModelNamesPrior,AllwoLofi_CTnames)
    
    
    shift = 0 #Lofi missing in prior -> index shifts
    
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.98*cm),
                            gridspec_kw={'width_ratios': [20, 4],'height_ratios':[1]})
    
    #Shading
    #monthly means
    ax2[0].plot([allSatModel.MonthDate.values[0],allSatModel.MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    
    ax2[0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    ax2[0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/GOSAT+IS',zorder = 4)# + str(round(ampl*factor, 2)),zorder = 4)
    
    ax2[0].plot(allModel.MonthDate,allModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label='In-situ-only inversions',zorder = 3)
    
    ax2[0].plot(PMIP.MonthDate,PMIP.Landtot*(44/12)*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label='MIP/OCO-2+IS')#+ str(round(ampl*(44/12)*factor, 2)),zorder = 4)
    #ax2[0].plot(PMIPTM54DVar.MonthDate,PMIPTM54DVar.Landtot*(44/12)*factor,ls='--',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm MIP_{TM5-4DVar}$',zorder = 4)
    
    ax2[0].plot(CAMSTM5BAKER.MonthDate,CAMSTM5BAKER.meannbp*(44/12)*factor,ls='-',linewidth = lw,marker= '',color= 'black',label='MIP/OCO-2+I'+r'$\rm S_{selected}$')#+ str(round(ampl*(44/12)*factor, 2)),zorder = 4)
    
    ax2[0].fill_between(CAMSTM5BAKER.MonthDate,
                        (CAMSTM5BAKER.minnbp)*(44/12)*factor,
                        (CAMSTM5BAKER.maxnbp)*(44/12)*factor,
                        linewidth = lw,color= 'black',alpha = 0.3,zorder = 2)#, label = '100%')#,label='All MIP wo LoFI',zorder = 2)
    
    
    #range
    ax2[0].fill_between(AllMIPwoLofi.MonthDate,
                        (AllMIPwoLofi.minnbp)*(44/12)*factor,
                        (AllMIPwoLofi.maxnbp)*(44/12)*factor,
                        linewidth = lw,color= 'black',alpha = 0.1)#, label = '100%')#,label='All MIP wo LoFI',zorder = 2)
    
    ax2[0].plot(PMIPIS.MonthDate,PMIPIS.Landtot*(44/12)*factor,ls=':',linewidth = lw,marker= '',color= 'black',label='MIP/IS',zorder = 3)
    
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]

    ax2[1].plot([0.5,12.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1].plot(range(1,13),allModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm Inverse~model_{in-situ}$',zorder = 3)
    ax2[1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    
    ax2[1].plot(range(1,13),PMIPMeanSais.Landtot.iloc[indexrange]*(44/12)*factor,ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 3)
    #ax2[1].plot(range(1,13),PMIPTM54DVarMeanSais.Landtot.iloc[indexrange]*(44/12)*factor,ls='--',linewidth = lw,marker= '',color= 'darkred',zorder = 4)
    ax2[1].plot(range(1,13),PMIPISMeanSais.Landtot.iloc[indexrange]*(44/12)*factor,ls=':',linewidth = lw,marker= '',color= 'black',zorder = 3)
    ax2[1].plot(range(1,13),(CAMSTM5BAKER.groupby(['Month'])[['TM5-4DVAR','Baker','CAMS']].mean().mean(axis = 1)).iloc[indexrange]*(44/12)*factor,ls='-',linewidth = lw,marker= '',color= 'black',zorder = 3)
    # SHADING 
    
    #ax2[1].fill_between(range(1,13),(allModelMeanSais.mean_nbp.iloc[indexrange]-allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allModelMeanSais.mean_nbp.iloc[indexrange]+allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'blue', alpha=0.3,zorder = 2)
    ax2[1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    
    ax2[1].fill_between(range(1,13),
                (AllMIPwoLofi.groupby(['Month'])[AllwoLofinames].mean().quantile(0,axis = 1)).iloc[indexrange]*(44/12)*factor,
                (AllMIPwoLofi.groupby(['Month'])[AllwoLofinames].mean().quantile(1,axis = 1)).iloc[indexrange]*(44/12)*factor,
                linewidth = lw,alpha=0.1,color= 'black',label=r'$\rm All~MIP~wo~LoFI$',zorder = 2)
    
    ax2[1].fill_between(range(1,13),
                (CAMSTM5BAKER.groupby(['Month'])[['TM5-4DVAR','Baker','CAMS']].mean().quantile(0,axis = 1)).iloc[indexrange]*(44/12)*factor,
                (CAMSTM5BAKER.groupby(['Month'])[['TM5-4DVAR','Baker','CAMS']].mean().quantile(1,axis = 1)).iloc[indexrange]*(44/12)*factor,
                linewidth = lw,alpha=0.3,color= 'black',label=r'$\rm All~MIP~wo~LoFI$',zorder = 2)
    
      
    # SETTINGS
    ax2[0].text(0.03, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=9)#, weight='bold')
    ax2[0].legend(fontsize=6,loc = 4,ncol = 3)#,loc=3)#, loc = 9)
    ##ax2[0].set_ylim(-1.65,1.8)
    ax2[0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0].tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2[0].tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2[0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    #ax2[0].set_ylim(-280,280)
    ax2[0].set_xlim(datetime.date(2014,12,15),datetime.date(2019,3,15))
    #ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2022,12,31))
    ax2[0].set_axisbelow(True)

    ax2[1].text(0.15, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=9)#, weight='bold')      
    secaxD = ax2[1].twinx()   
    secaxD.set_ylim(ax2[0].get_ylim())
    ##secaxD.set_ylim(-1.65,1.8)
    secaxD.set_yticklabels([])
    secaxD.tick_params(axis = 'y',direction = 'in')
    ax2[1].set_xticks(ticks = [3,6,9,12])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1,1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    ax2[1].set_ylim(ax2[0].get_ylim())
    ##ax2[1].set_ylim(-1.65,1.8)
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6, zorder = 0)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig("/Results/Plots/Africa/Paper2024/Figure5_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')


# Figure 6 GOME SIF TRENDY
if False:
    
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
    lw = 1
    cm = 1/2.54
    fig, ax = plt.subplots(figsize = (8*cm, 6*cm))
    colorl = ['grey','green','blue','red']
    for num, mod in enumerate(['ORCHIDEE','ORCHIDEEv3','CABLE-POP','OCN']):
        ax.plot(range(1,13),NormalizeMeanCycle(MeanCycle,mod+'gpp').iloc[indexrange],lw=lw,ls=':',marker= '',color= colorl[num],label=mod +' gpp')
    _ , AmplFac = NormalizeMeanCycle(getMSC(PMonth_means_GOMESIFv2grid,GOMEvar,year_min,year_max),GOMEvar,True)
    ax.errorbar(range(1,13),NormalizeMeanCycle(getMSC(PMonth_means_GOMESIFv2grid,GOMEvar,year_min,year_max),GOMEvar).iloc[indexrange],getMSC(PMonth_means_GOMESIFv2grid,'StDev',year_min,year_max)['StDev'].iloc[indexrange]/AmplFac,lw=lw,ls='',marker= '.',color= 'black',capsize=2,label=r'GOME-2 SIF')
    #add shading
    ax.fill_between(range(1,13),NormalizeMeanCycle(getMSC(PMonth_means_GOMESIFv2grid,GOMEvar,year_min,year_max),GOMEvar).iloc[indexrange]-getMSC(PMonth_means_GOMESIFv2grid,'StDev',year_min,year_max)['StDev'].iloc[indexrange]/AmplFac,
                                    NormalizeMeanCycle(getMSC(PMonth_means_GOMESIFv2grid,GOMEvar,year_min,year_max),GOMEvar).iloc[indexrange]+getMSC(PMonth_means_GOMESIFv2grid,'StDev',year_min,year_max)['StDev'].iloc[indexrange]/AmplFac,
                                    color= 'black',alpha=0.05)    
    
    ax.set_ylabel('Normalized',fontsize = 7)
    ax.set_xlabel('Month',fontsize = 7)
    ax.legend(loc=4,fontsize=6)
    ax.tick_params(axis = 'both',labelsize = 7, zorder = 0)
    ax.set_xticks(ticks = [3,6,9,12])
    ax.set_xticklabels(['9','12','3','6'])
    #ax.set_ylim(-3.5,4.5)
    #ax.set_yticks(ticks = [-1,0,1])
    #ax.set_yticklabels(['-1','0','1'])
    ax.axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    
    plt.savefig("/Results/Plots/Africa/Paper2024/GOME2SIFv2_"+GOMEvar+"_Norm_TrendyWocn_wStDevgrid_"+RegionName+str(Numm)+"_"+str(year_min)+"_"+str(year_max)+".png", dpi=300,bbox_inches='tight')


# Figure 7 Paper with Junly - June Year 
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(2, 2, 
                            figsize = (18.0*cm,12.2*cm),
                            gridspec_kw={'width_ratios': [4, 1],'height_ratios':[1,1]})

    #get std of MSC
    allSatGFED = pd.merge(allSatModel, PGFEDregion[['total_emission','Year','Month']], on = ['Year','Month'])
    allSatGFED.insert(0,column = 'NEE',value = allSatGFED.mean_nbp - allSatGFED.total_emission)
    StdCycallSatGFED = allSatGFED.groupby('Month')['NEE'].std(ddof=0).reset_index()
    AllMMTrendyNBP.insert(0,column ='meanneegoodModel',value = AllMMTrendyNBP.meanrhragoodModel-AllMMTrendyNBP.meangppgoodModel)
    StdCycgoodnee = AllMMTrendyNBP.groupby('Month')['meanneegoodModel'].std(ddof=0).reset_index()
    StdCycgoodnbp = AllMMTrendyNBP.groupby('Month')['meanNBPgoodModel'].std(ddof=0).reset_index()
    
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
    

     # Panel a)
    ax2[0,0].text(0.03, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0,0].transAxes,fontsize=10, weight='bold')
    # SHADING
    ## stDev for TRENDY
    #ax2[0,0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPbadModel-AllMMTrendyNBP.stdNBPbadModel)*factor,-(AllMMTrendyNBP.meanNBPbadModel+AllMMTrendyNBP.stdNBPbadModel)*factor,color= 'grey', alpha=0.3,zorder = 2)#,label='bad TRENDY')#,label = 'StDev')
    ax2[0,0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPgoodModel-AllMMTrendyNBP.stdNBPgoodModel)*factor,-(AllMMTrendyNBP.meanNBPgoodModel+AllMMTrendyNBP.stdNBPgoodModel)*factor,color= 'black', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    ## min max for satellite
    ax2[0,0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)     
    #monthly means
    ax2[0,0].plot([allSatModel.MonthDate.values[0],allSatModel.MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[0,0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBPgoodModel*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{selection}$',zorder = 4)
    #ax2[0,0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBPbadModel*factor,ls='-',linewidth = lw,marker= '',color= 'dimgrey',label=r'$\rm TRENDY_{others}$',zorder = 4)
    ax2[0,0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/GOSAT+IS',zorder = 4)
    

     #Panel b)
    ax2[0,1].text(0.15, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0,1].transAxes,fontsize=10, weight='bold')
    
    ax2[0,1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[0,1].plot(range(1,13),-MeanGood.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{all}$',zorder = 4)
    #ax2[0,1].plot(range(1,13),-MeanBad.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ax2[0,1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/GOSAT+IS',zorder = 4)
    # SHADING
    ax2[0,1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    #ax2[0,1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StDevGood.nbp.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StDevGood.nbp.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    #ax2[0,1].fill_between(range(1,13),-(MeanBad.nbp.iloc[indexrange]-StDevBad.nbp.iloc[indexrange])*factor,-(MeanBad.nbp.iloc[indexrange]+StDevBad.nbp.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[0,1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    #ax2[0,1].fill_between(range(1,13),-(MeanBad.nbp.iloc[indexrange]-StdCycbadnbp.meanNBPbadModel.iloc[indexrange])*factor,-(MeanBad.nbp.iloc[indexrange]+StdCycbadnbp.meanNBPbadModel.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    

    # Panel c)
    #monthly means
    
    # SHADING
    ax2[1,0].text(0.03, 0.9, r'$\rm \bf{C}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1,0].transAxes,fontsize=10, weight='bold')
    ## stDev for TRENDY
    ax2[1,0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPgoodModel-AllMMTrendyNBP.stdNBPgoodModel)*factor,-(AllMMTrendyNBP.meanNBPgoodModel+AllMMTrendyNBP.stdNBPgoodModel)*factor,color= 'black', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    ## min max for satellite
    ax2[1,0].fill_between(allSatGFED.MonthDate,(allSatModel.max_nbp-allSatGFED.total_emission)*factor,(allSatModel.min_nbp-allSatGFED.total_emission)*factor,color= 'red', alpha=0.3,zorder = 2)     
    
    ax2[1,0].plot([allSatModel.MonthDate.values[0],allSatModel.MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1,0].plot(AllMMTrendyNBP.MonthDate,AllMMTrendyNBP.meanneegoodModel*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{selection}~NEE$',zorder = 4)
    ax2[1,0].plot(allSatGFED.MonthDate,allSatGFED.NEE*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/GOSAT+IS - GFED',zorder = 4)
    
    #annual bar plot
    secAx = ax2[1,0].twinx()
    MeanAf = allSatModel[(allSatModel.JJYear >= 2009)&(allSatModel.JJYear <= 2017)].groupby(['JJYear'])['mean_nbp'].sum().reset_index()
    MeanAf = pd.merge(MeanAf, PGFEDregion[(PGFEDregion.JJYear >= 2009)&(PGFEDregion.JJYear <= 2017)].groupby(['JJYear'])['total_emission'].sum().reset_index(),on = ['JJYear'])
    MeanAfTrendyGood = AllMMTrendyNBP[(AllMMTrendyNBP.JJYear >= 2009)&(AllMMTrendyNBP.JJYear <= 2017)].groupby(['JJYear'])['meanrhragoodModel'].sum().reset_index()
    MeanAfTrendyGood = pd.merge(MeanAfTrendyGood, AllMMTrendyNBP[(AllMMTrendyNBP.JJYear >= 2009)&(AllMMTrendyNBP.JJYear <= 2017)].groupby(['JJYear'])['meangppgoodModel'].sum().reset_index(),on = ['JJYear'])
    MeanAf.insert(0, column = 'JJYearDate', value = MeanAf.apply(lambda x: datetime.date(int(x.JJYear),12,31), axis = 1))
    MeanAfTrendyGood.insert(0, column = 'JJYearDate', value = MeanAfTrendyGood.apply(lambda x: datetime.date(int(x.JJYear),12,31), axis = 1))
    #secAx.bar(MeanAf.JJYearDate, MeanAf.mean_nbp*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'red',alpha = 0.3)
    #secAx.bar(MeanAfTrendyGood.JJYearDate, (-1)*MeanAfTrendyGood.meanNBPgoodModel*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'black',alpha = 0.3)
    secAx.bar(MeanAf.JJYearDate, (MeanAf.mean_nbp-MeanAf.total_emission)*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'red',alpha = 0.3)
    secAx.bar(MeanAfTrendyGood.JJYearDate, (MeanAfTrendyGood.meanrhragoodModel-MeanAfTrendyGood.meangppgoodModel)*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'black',alpha = 0.3)
    
    #Panel d)
    ax2[1,1].text(0.15, 0.9, r'$\rm \bf{D}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1,1].transAxes,fontsize=10, weight='bold')
    
    ax2[1,1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1,1].plot(range(1,13),(MeanGood.rarh.iloc[indexrange]-MeanGood.gpp.iloc[indexrange])*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{good}$',zorder = 4)
    ax2[1,1].plot(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-GFEDMonthlyAnnualMean.total_emission.iloc[indexrange])*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label='TM5-4DVar/GOSAT+IS - GFED',zorder = 4)
    # SHADING
    ax2[1,1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]-StdCycallSatGFED.NEE.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]-GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]+StdCycallSatGFED.NEE.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    #ax2[1,1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StDevGood.nbp.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StDevGood.nbp.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    #ax2[1,1].fill_between(range(1,13),-(MeanBad.nbp.iloc[indexrange]-StDevBad.nbp.iloc[indexrange])*factor,-(MeanBad.nbp.iloc[indexrange]+StDevBad.nbp.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[1,1].fill_between(range(1,13),(MeanGood.rarh.iloc[indexrange]-MeanGood.gpp.iloc[indexrange]-StdCycgoodnee.meanneegoodModel.iloc[indexrange])*factor,(MeanGood.rarh.iloc[indexrange]-MeanGood.gpp.iloc[indexrange]+StdCycgoodnee.meanneegoodModel.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    
    #SETTINGS
    #Panel a)
    ax2[0,0].legend(fontsize=6,loc = 4,ncol = 3)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[0,0].grid(True,which = 'both', axis='x',zorder = 0)
    #ax2[0,0].set_ylim(-280,280)
    ax2[0,0].set_ylim(-350,350)
    ax2[0,0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0,0].set_xticklabels([])
    #ax2[0,0].set_yticklabels(ax2[0,0].get_yticks(),fontsize=7)
    ax2[0,0].tick_params(axis = 'x',length = 0.01, zorder = 0)
    ax2[0,0].tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2[0,0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[0,0].set_axisbelow(True)

    
    #Panel b)  
    secaxB = ax2[0,1].twinx()     
    #secaxB.set_ylim(-280,280)
    secaxB.set_ylim(-350,350)
    secaxB.set_yticklabels([])
    secaxB.tick_params(axis = 'y',direction = 'in')          
    ax2[0,1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    ax2[0,1].set_ylim(-350,350)
    #ax2[0,1].set_ylim(-280,280)
    ax2[0,1].set_xticklabels([])
    ax2[0,1].set_yticklabels([])
    ax2[0,1].tick_params(axis = 'both',length = 0.01, zorder = 0)
    ax2[0,1].set_xlim(0.5,12.5)
    ax2[0,1].set_axisbelow(True)
    
    
    #Panel c)
    ax2[1,0].legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[1,0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[1,0].set_xlabel(r'Date',fontsize=7)
    ax2[1,0].set_ylim(-350,350)
    ax2[1,0].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[1,0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[1,0].set_axisbelow(True)
    secAx.set_ylim(-500,500)
    ax2[1,0].yaxis.set_label_position("right")
    ax2[1,0].yaxis.tick_right()
    secAx.yaxis.set_label_position("left")
    secAx.yaxis.tick_left()
    secAx.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    #secAx.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    secAx.tick_params(axis = 'y',labelsize = 6, zorder = 0)
        
    secaxC = ax2[1,0].twinx() 
    secaxC.set_ylim(-350,350)
    secaxC.set_yticklabels([])
    secaxC.tick_params(axis = 'y',direction = 'in')  

    ax2[1,0].set_yticklabels([])
    ax2[1,0].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1,0].tick_params(axis = 'x',labelsize = 6)
    
    #Panel d) 
    ax2[1,1].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[1,1].yaxis.set_label_position("right")
    ax2[1,1].yaxis.tick_right()        
         
    ax2[1,1].set_xticks(ticks = [3,6,9,12])
    ax2[1,1].set_xticklabels(['9','12','3','6'])
    ax2[1,1].tick_params(axis = 'both',labelsize = 6)
    
    ax2[1,1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1,1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1,1].set_xlabel(r'Month',fontsize=7)
    #ax2[1,1].set_ylim(-280,280)
    ax2[1,1].set_ylim(-350,350)
    ax2[1,1].set_xlim(0.5,12.5)
    ax2[1,1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig("/Results/Plots/Africa/Paper2024/Fig7_completeABCD_JJAnnualNEE"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')

# Figure 8, needs to be run three times with different regions
if False:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.1*cm),#27
                            gridspec_kw={'width_ratios': [1, 2]})
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #secAx = ax2.twinx()
    
    norm = ''#'Norm'
    timerange = [7,8,10,12,2,4,6]
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
        
    factor = 12/44*1000

    clustertype = MeanGood.copy() #MeanGood, MeanBad, MeanLateResp, MeanEarlyGPP
    clustertypeerror = StDevGood.copy() #StDevGood, StDevBad, StDevLateResp, StDevEarlyGPP
    AllMMTrendyNBP.insert(0,column ='meanneegoodModel',value = AllMMTrendyNBP.meanrhragoodModel-AllMMTrendyNBP.meangppgoodModel)
    StdCycgoodnee = AllMMTrendyNBP.groupby('Month')['meanneegoodModel'].std(ddof=0).reset_index()

    #Panel A
    ax2[0].text(0.10, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10, weight='bold')
    ax2[0].plot([1,12],[0,0],color = 'grey',lw=lw,ls='-', marker='')
    ax2[0].plot(range(1,13),-clustertype['nbp'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'black',label='NBP')
    ax2[0].plot(range(1,13),-(clustertype['gpp'+norm].iloc[indexrange]-clustertype['rarh'+norm].iloc[indexrange])*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'black',label='NEE')
    #ax2[0].plot(range(1,13),clustertype['gpp'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'green',label='GPP')
    #ax2[0].plot(range(1,13),clustertype['ra'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'grey',label='ra')
    ax2[0].plot(range(1,13),clustertype['rh'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'slateblue',label='rh')
    #ax2[0].plot(range(1,13),(clustertype['rarh'+norm]).iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'brown',label='rh+ra')
    ax2[0].plot(range(1,13),(clustertype['npp'+norm]).iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'green',label='GPP-RA')
    
    #ax2[0].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*10**(-3)*factor,ls=':',marker= '',color= 'darkred',label=r'$\rm Inverse~Model_{GOSAT}$')

    ax2[0].set_xticklabels(timerange)
        
    
    var = 'nbp'+norm
    ax2[0].fill_between(range(1,13),
                            (-clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (-clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'grey',alpha = 0.4)
    
    
    
    var = 'rh'+norm
    ax2[0].fill_between(range(1,13),
                            (clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'slateblue',alpha = 0.4)
    
    var = 'npp'+norm
    ax2[0].fill_between(range(1,13),
                            (clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'green',alpha = 0.4)
   
    # Panel B
    ax2[1].text(0.06, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10, weight='bold')
    
    MeanAfTrendyRAGood = AllMMTrendyNBP.groupby(['JJYear'])['meanragoodModel'].sum().reset_index()
    MeanAfTrendyGPPGOod = AllMMTrendyNBP.groupby(['JJYear'])['meangppgoodModel'].sum().reset_index()
    MeanAfTrendyNBPGood = AllMMTrendyNBP.groupby(['JJYear'])['meanNBPgoodModel'].sum().reset_index()
    MeanAfTrendyRHGood = AllMMTrendyNBP.groupby(['JJYear'])['meanrhgoodModel'].sum().reset_index()
    MeanAfTrendyGPPRAGood = MeanAfTrendyGPPGOod.copy(deep = True)
    MeanAfTrendyGPPRAGood = MeanAfTrendyGPPRAGood.merge(MeanAfTrendyRAGood,how = 'left', on = 'JJYear')
    MeanAfTrendyGPPRAGood.insert(0,column = 'meangpp-ragoodModel', value = MeanAfTrendyGPPRAGood.meangppgoodModel - MeanAfTrendyGPPRAGood.meanragoodModel)
    
    MeanAfTrendyNEPGood = MeanAfTrendyGPPRAGood.copy(deep = True)
    MeanAfTrendyNEPGood = MeanAfTrendyNEPGood.merge(MeanAfTrendyRHGood,how = 'left', on = 'JJYear')
    MeanAfTrendyNEPGood.insert(0,column = 'meannepgoodModel', value = MeanAfTrendyNEPGood['meangpp-ragoodModel'] - MeanAfTrendyNEPGood.meanrhgoodModel)
    
    #data to be ploted
    plotdata = [MeanAfTrendyNEPGood,MeanAfTrendyGPPRAGood,MeanAfTrendyRHGood]#,MeanAfTrendyNBPGood]#MeanAfTrendyRAGood, MeanAfTrendyGPPGOod
    plotvar = ['meannepgoodModel','meangpp-ragoodModel','meanrhgoodModel']#,'meanNBPgoodModel']#'meanragoodModel','meangppgoodModel'
    labels = ['NEE','GPP-RA','RH']#,'NBP']#'RA','GPP'
    colors = ['black','darkgreen','slateblue']#,'grey']
    faclist = [-12/44,-12/44,12/44]#,-12/44]#12/44,-12/44
    
    width = 0.8/len(plotdata)
    for y in range(2009,2018):
        for i,pd in enumerate(plotdata):
            offset = width * i
            pd = pd[(pd.JJYear >= 2009)& (pd.JJYear <= 2017)]
            if y == 2009: #plot label
                #anomaly
                ax2[1].bar(y + offset,  (pd[(pd.JJYear == y)][plotvar[i]]-pd[plotvar[i]].mean())*faclist[i], width, color = colors[i],label = labels[i])
                #absolute
                #ax2[1].bar(y + offset,  (pd[(pd.JJYear == y)][plotvar[i]])*faclist[i], width, color = colors[i],label = labels[i])
                
            else:
                ax2[1].bar(y + offset,  (pd[(pd.JJYear == y)][plotvar[i]]-pd[plotvar[i]].mean())*faclist[i], width, color = colors[i])
                #ax2[1].bar(y + offset,  (pd[(pd.JJYear == y)][plotvar[i]])*faclist[i], width, color = colors[i])

        #plot NBP hatched on top
        if y == 2009: #plot label
            pd = MeanAfTrendyNBPGood[(MeanAfTrendyNBPGood.JJYear >= 2009)& (MeanAfTrendyNBPGood.JJYear <= 2017)]
            ax2[1].bar(y ,  (pd[(pd.JJYear == y)]['meanNBPgoodModel']-pd['meanNBPgoodModel'].mean())*faclist[0], width*0.9, facecolor = 'none',edgecolor= 'grey',hatch = '////',label = 'NBP')
        else:
            pd = MeanAfTrendyNBPGood[(MeanAfTrendyNBPGood.JJYear >= 2009)& (MeanAfTrendyNBPGood.JJYear <= 2017)]
            ax2[1].bar(y ,  (pd[(pd.JJYear == y)]['meanNBPgoodModel']-pd['meanNBPgoodModel'].mean())*faclist[0], width*0.9, facecolor = 'none',edgecolor= 'grey',hatch = '////')
            #ax2[1].bar(y + offset,  (pd[(pd.JJYear == y)][plotvar[i]])*faclist[i], width, color = colors[i])

    #SETTINGS
    #ax2[0].set_ylim([-120,300])
    #ax2[0].set_ylim([-230,600])
    ax2[0].set_xlabel(r'Month',fontsize=7)
    #ax2[0].set_ylabel(r'$\rm flux~[\mu gCO_2/(m^2 s)]$',fontsize=15)
    ax2[0].set_ylabel(r'$\rm CO_2~flux~[TgC/month]$',fontsize=7)
    #ax2[0].set_title('Comparison Cluster models')
    ax2[0].set_zorder(1)  # default zorder is 0 for ax1 and ax2[0]
    ax2[0].patch.set_visible(False)  # prevents ax1 from hiding ax2[0]
    ###ax2[0].yaxis.set_label_position("right")
    ax2[0].tick_params(axis = 'x',labelsize = 6)
    ax2[0].tick_params(axis = 'y',labelsize = 6)
    
    
    ax2[1].set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    #ax2[1].tick_params(axis = 'y', zorder = 0,labelsize = 6)
    ax2[1].legend(fontsize=6,ncol = 2,loc = 8)#,loc=3)#, loc = 9)
    ax2[1].set_xticks(np.array(range(2009,2018))+width)
    ax2[1].set_xticklabels(['2010','2011','2012','2013','2014','2015','2016','2017','2018'])
    ax2[1].tick_params(axis = 'x',labelsize = 6)
    ax2[1].tick_params(axis = 'y',labelsize = 6)
    #ax2[1].set_ylim(-180,150)
    #ax2[1].set_ylim(-200,230)
    
    ax2[1].set_xlabel(r'Year',fontsize=7)
    ax2[1].yaxis.set_label_position("right")
    ax2[1].yaxis.tick_right() 
    
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig("/Results/Plots/Africa/Paper2024/Figure8_"+RegionName+str(Numm)+'_'+str(year_min)+'_'+str(year_max)+".png", dpi=400,bbox_inches='tight')
    

#####SUPPLEMENT####
#S1 and S2 in Plot_TimeSeries.py
    
#S4 Fire
if False:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    factor =  12/44
    
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6.98*cm))
    
    ax2.plot(PMonth_means_GFAS.MonthDate,PMonth_means_GFAS.CO2fireE*10**(-12)*factor,color = 'red',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'GFAS')
    ax2.plot(PGFEDregion.MonthDate,PGFEDregion.total_emission*factor,color = 'orange',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'GFED')
    ax2.plot(PMonth_means_FINN.MonthDate,PMonth_means_FINN.CO2fireE*10**(-12)*factor,color = 'purple',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'FINN')
    ax2.set_ylabel(r'$\rm Fire~CO_2~emission~(TgC/month)$',fontsize=7)
    ax2.tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2.tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.legend(fontsize=6,loc = 1)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_axisbelow(True)

    plt.savefig("/Results/Plots/GFED_GFAS_FINN_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".png", dpi=400, bbox_inches = "tight")
    

#S5 MIP OCO-2 cosamples created with Plot_TimeSeries in Scripts_PhD.py

# S5 Annual NBP
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6.1*cm))

    MeanAf = allSatModel[(allSatModel.Year >= 2010)].groupby(['Year'])['mean_nbp'].sum().reset_index()
    MeanAf = pd.merge(MeanAf, PGFEDregion[(PGFEDregion.Year >= 2010)].groupby(['Year'])['total_emission'].sum().reset_index(),on = ['Year'])
    MeanAfTrendyGood = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= 2010)].groupby(['Year'])['meanrhragoodModel'].sum().reset_index()
    MeanAfTrendyGood = pd.merge(MeanAfTrendyGood, AllMMTrendyNBP[(AllMMTrendyNBP.Year >= 2010)].groupby(['Year'])['meangppgoodModel'].sum().reset_index(),on = ['Year'])
    MeanAf.insert(0, column = 'YearDate', value = MeanAf.apply(lambda x: datetime.date(int(x.Year),6,30), axis = 1))
    MeanAfTrendyGood.insert(0, column = 'YearDate', value = MeanAfTrendyGood.apply(lambda x: datetime.date(int(x.Year),6,30), axis = 1))
    #secAx.bar(MeanAf.YearDate, MeanAf.mean_nbp*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'red',alpha = 0.3)
    #secAx.bar(MeanAfTrendyGood.YearDate, (-1)*MeanAfTrendyGood.meanNBPgoodModel*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'black',alpha = 0.3)
    ax2.bar(MeanAf.YearDate, (MeanAf.mean_nbp)*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'red',alpha = 0.3,label = r'$\rm Inverse~model_{+GOSAT}$' )
    ax2.bar(MeanAfTrendyGood.YearDate, (MeanAfTrendyGood.meanrhragoodModel-MeanAfTrendyGood.meangppgoodModel+MeanAf.total_emission)*factor, width = datetime.date(int(2010),12,31)-datetime.date(int(2010),1,1), color = 'black',alpha = 0.3, label = r'$\rm TRENDY_{selection}~NEE$ + GFED')
    
    # SETTINGS
    #Panel c)
    ax2.legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)

    ax2.tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2009,9,15),datetime.date(2019,3,15))
    ax2.set_axisbelow(True)
    
    ax2.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    #secAx.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    ax2.tick_params(axis = 'y',labelsize = 6, zorder = 0)

    plt.savefig("/Results/Plots/Africa/Paper2024/FigS5X_AnnualNBP"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')


#S6 ERA5 Anomalies
if False:
    if Numm not in [756,772,783]:
        print('not implemented for other regions yet')
        sys.exit()
    ERA5 = pd.read_pickle("/ERA5/Africa/DataFrames/GDF1_ERA5Precip_Af"+str(Numm)+".pkl")
    ERA5t = pd.read_pickle("/ERA5/Africa/DataFrames/GDF1_t2m_Af"+str(Numm)+".pkl")
    
    #precipitation 
    MonthSumERA5 = getMonthMeans(ERA5,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'monthlytp')
    date = MonthSumERA5.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthSumERA5.insert(loc=1,column='MonthDate',value=date)
    #get annual values
    AnnualSumERA5 = MonthSumERA5.groupby(['Year'])['monthlytp'].sum().reset_index()

    MonthMeanERA5t = getMonthMeans(ERA5t,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'t2m')
    date = MonthMeanERA5t.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthMeanERA5t.insert(loc=1,column='MonthDate',value=date)
    MeanMonthERA5t = MonthMeanERA5t.groupby(['Month'])['t2m'].mean().reset_index()

    #get annual values
    AnnualMeanERA5t = MonthMeanERA5t.groupby(['Year'])['t2m'].mean().reset_index()
    
    #plot data as bars
    fig2, ax2 = plt.subplots(figsize = (7,2))
    #ax2.bar(AnnualSumERA5.Year,AnnualSumERA5.monthlytp*1000)#,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'Blue',label='Total area')
    ax2.bar(AnnualSumERA5.Year,AnnualSumERA5.monthlytp*1000-(AnnualSumERA5.monthlytp.mean()*1000),alpha = 0.3,zorder = 2,label = 'Precip.')#,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'Blue',label='Total area')
    if True:
        secAx = ax2.twinx()
        secAx.bar(AnnualMeanERA5t.Year,AnnualMeanERA5t.t2m-AnnualMeanERA5t.t2m.mean(),color = 'red',alpha = 0.3,zorder = 2,label = 'Temp.')#,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'Blue',label='Total area')
        secAx.legend(fontsize=12,loc=4)#,loc=3)#, loc = 9)
        secAx.set_ylabel('Annual Temperature \n Anomaly (K)',fontsize=12)
        #secAx.set_ylim(-1.1,1.1)
        #ax2.set_ylim(-180,180)
        secAx.set_ylim(-0.5,0.5)
        ax2.set_ylim(-70,70)
        ax2.plot([2009.5,2018.5],[0,0],color = 'grey',ls = '-',marker='',zorder = 0)
    ax2.legend(fontsize=12,loc=3)#,loc=3)#, loc = 9)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Year',fontsize=12)
    #ax2.set_ylabel('Annual \n precipitation \n (mm)',fontsize=12)
    ax2.set_ylabel('Annual precipitation \n Anomaly (mm)',fontsize=12)
    ax2.tick_params(axis = 'both',labelsize = 12, zorder = 0)
    ax2.set_axisbelow(True) #to put grid behind the plots, zorder in grid does not work as grid is part of the axis
    
    plt.savefig("/Results/Plots/Africa/Paper2024/FigureS6_"+RegionName+str(Numm)+"_"+str(year_min)+"-"+str(year_max)+".png", dpi=250,bbox_inches='tight')

#S7 ERA5 Monthly
if False:
    if Numm not in [756,772,783]:
        print('not implemented for other regions yet')
        sys.exit()
    ERA5 = pd.read_pickle("/ERA5/Africa/DataFrames/GDF1_ERA5Precip_Af"+str(Numm)+".pkl")
    ERA5t = pd.read_pickle("/ERA5/Africa/DataFrames/GDF1_t2m_Af"+str(Numm)+".pkl")
    
    #precipitation 
    MonthSumERA5 = getMonthMeans(ERA5,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'monthlytp')
    date = MonthSumERA5.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthSumERA5.insert(loc=1,column='MonthDate',value=date)

    MonthMeanERA5t = getMonthMeans(ERA5t,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'t2m')
    date = MonthMeanERA5t.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthMeanERA5t.insert(loc=1,column='MonthDate',value=date)
    
    #plot data as bars
    fig2, ax2 = plt.subplots(figsize = (10,5))
    ax2.bar(MonthSumERA5.MonthDate,MonthSumERA5.monthlytp*1000,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'Blue',label='Precipitation')
    if True:
        secAx = ax2.twinx()
        secAx.plot(MonthMeanERA5t.MonthDate,MonthMeanERA5t.t2m,color = 'red',label='Temp.')
        secAx.legend(fontsize=12,loc=2,ncol = 3)#,loc=3)#, loc = 9)
        secAx.set_ylabel('Mean Temperature',fontsize=12)
    
    ax2.legend(fontsize=12,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=12)
    ax2.set_ylabel('Mean monthly \n precipitation \n (mm)',fontsize=12)
    ax2.tick_params(axis = 'both',labelsize = 12, zorder = 0)
    ax2.set_axisbelow(True) #to put grid behind the plots, zorder in grid does not work as grid is part of the axis
    
    plt.savefig("/Results/Plots/Africa/Paper2024/FigureS7_"+RegionName+str(Numm)+".png", dpi=250,bbox_inches='tight')

# S8 MSC TRENDY selection with GPP and ra
if False:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (9.0*cm,6.1*cm))
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #secAx = ax2.twinx()
    
    norm = ''#'Norm'
    timerange = [7,8,10,12,2,4,6]
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
        
    factor = 12/44*1000

    clustertype = MeanGood.copy() #MeanGood, MeanBad, MeanLateResp, MeanEarlyGPP
    clustertypeerror = StDevGood.copy() #StDevGood, StDevBad, StDevLateResp, StDevEarlyGPP
    AllMMTrendyNBP.insert(0,column ='meanneegoodModel',value = AllMMTrendyNBP.meanrhragoodModel-AllMMTrendyNBP.meangppgoodModel)
    StdCycgoodnee = AllMMTrendyNBP.groupby('Month')['meanneegoodModel'].std(ddof=0).reset_index()

    #Panel A
    #ax2.text(0.10, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=10, weight='bold')
    ax2.plot([1,12],[0,0],color = 'grey',lw=lw,ls='-', marker='')
    ax2.plot(range(1,13),-clustertype['nbp'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'black',label='NBP')
    ax2.plot(range(1,13),-(clustertype['gpp'+norm].iloc[indexrange]-clustertype['rarh'+norm].iloc[indexrange])*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'black',label='NEP')
    ax2.plot(range(1,13),clustertype['gpp'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'green',label='GPP')
    ax2.plot(range(1,13),clustertype['ra'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls=':',marker= '',color= 'grey',label='ra')
    ax2.plot(range(1,13),clustertype['rh'+norm].iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'slateblue',label='rh')
    #ax2.plot(range(1,13),(clustertype['rarh'+norm]).iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'brown',label='rh+ra')
    ax2.plot(range(1,13),(clustertype['npp'+norm]).iloc[indexrange]*10**(-3)*factor,lw=lw,ls='-',marker= '',color= 'green',label='GPP-RA')
    
    #ax2.plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*10**(-3)*factor,ls=':',marker= '',color= 'darkred',label=r'$\rm Inverse~Model_{GOSAT}$')

    ax2.set_xticklabels(timerange)

    var = 'nbp'+norm
    ax2.fill_between(range(1,13),
                            (-clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (-clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'grey',alpha = 0.4)
    
    var = 'rh'+norm
    ax2.fill_between(range(1,13),
                            (clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'slateblue',alpha = 0.4)
    
    var = 'npp'+norm
    ax2.fill_between(range(1,13),
                            (clustertype[var].iloc[indexrange]-clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            (clustertype[var].iloc[indexrange]+clustertypeerror[var].iloc[indexrange])*10**(-3)*factor,
                            ls='-',color= 'green',alpha = 0.4)
   
    
    #SETTINGS
    ax2.legend(fontsize=6,loc='center left', bbox_to_anchor=(1, 0.5))
    #ax2.set_ylim([-120,300])
    #ax2.grid(True,which = 'both', axis='x')
    ax2.set_xlabel(r'Month',fontsize=7)
    #ax2.set_ylabel(r'$\rm flux~[\mu gCO_2/(m^2 s)]$',fontsize=15)
    ax2.set_ylabel(r'$\rm CO_2~flux~[TgC/month]$',fontsize=7)
    #ax2.set_title('Comparison Cluster models')
    ax2.set_zorder(1)  # default zorder is 0 for ax1 and ax2
    ax2.patch.set_visible(False)  # prevents ax1 from hiding ax2
    ###ax2.yaxis.set_label_position("right")
    ax2.tick_params(axis = 'x',labelsize = 6)
    ax2.tick_params(axis = 'y',labelsize = 6)
    
    plt.savefig("/Results/Plots/Africa/Paper2024/FigureS8_"+RegionName+str(Numm)+'_'+str(year_min)+'_'+str(year_max)+".png", dpi=400,bbox_inches='tight')
    
#S9 y axis split: need to be run with region 772 and 784    
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    mpl.rcParams['hatch.linewidth'] = 0.5
    fig2, ax2 = plt.subplots(2, 1, 
                            figsize = (18.0*cm,8*cm),
                            gridspec_kw={'width_ratios': [1],'height_ratios':[1,1]})
    
    clusterkind = 'good'

    ax2[1].plot([AllMMTrendyNBP.MonthDate.values[0],AllMMTrendyNBP.MonthDate.values[-1]],[0,0],color = 'grey',ls='-', marker='')
        
    y_gpp = AllMMTrendyNBP['meangpp-ra'+clusterkind+'Model']*factor
    y_rarh = (AllMMTrendyNBP['meanrh'+clusterkind+'Model'])*factor
    print(AllMMTrendyNBP.MonthDate)
    
    ax2[0].fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp,
                        ls='-',linewidth = lw,facecolor = 'none',edgecolor= 'green',hatch='//',
                        label = 'GPP-RA')
    ax2[0].fill_between(AllMMTrendyNBP.MonthDate,
                        y_rarh,
                        ls='-',linewidth = lw,facecolor = 'none',edgecolor= 'indigo',hatch='\\\\',
                        label = 'soil respiration')
    ax2[0].fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp, y_rarh,
                        ls='-',linewidth = lw,where =(y_gpp>y_rarh),color= 'green',
                        alpha = 0.4,interpolate=True)
    ax2[0].fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp,y_rarh,
                        ls='-',linewidth = lw,where =(y_gpp<y_rarh),color= 'indigo',
                        alpha = 0.6,interpolate=True)
    
    #y1 = -AllMMTrendyNBP['meanNBP'+clusterkind+'Model']*factor#*10**9*44/12
    y1 = -(AllMMTrendyNBP['meangpp-ra'+clusterkind+'Model']-AllMMTrendyNBP['meanrh'+clusterkind+'Model'])*factor
    #yStd = AllMMTrendyNBP['stdNBP'+clusterkind+'Model']*factor
    ax2[1].plot(AllMMTrendyNBP.MonthDate,y1,linewidth = lw,ls='-',marker= '',color= 'black',label='GPP - TER')
    
    ax2[1].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1>0),color= 'indigo',
                        alpha = 0.6,linewidth = lw,interpolate=True)
    ax2[1].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1<0),color= 'green'
                        ,alpha = 0.4,linewidth = lw,interpolate=True)
    
    #ax2[1].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='--',marker= '',color= 'darkred',label=r'$\rm Inverse~Model_{+GOSAT}$')
    ax2[1].plot(allSatModel.MonthDate,(allSatModel.mean_nbp-PGFEDregion.total_emission)*factor,ls='--',marker= '',color= 'darkred',label=r'$\rm Inverse~Model_{+GOSAT} - GFED$')
    
    ax2[0].legend(fontsize=6,loc = 1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[1].legend(fontsize=6,loc = 1,ncol = 3)#,loc=3)#, loc = 9)
    
    ax2[0].set_ylim(0,300)
    ax2[1].set_ylim(-150,150)
    ax2[0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[1].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0].set_xticklabels([])
    ax2[1].set_xlabel(r'Date',fontsize=7)
    #ax2.set_xticklabels(ax2.get_xticks(),fontsize=7)
    #ax2.set_yticklabels(ax2.get_yticks(),fontsize=7)
    ax2[0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[1].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0].tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[0].set_axisbelow(True)
    ax2[1].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[1].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig("/Results/Plots/Africa/Flux/TRENDYgood_GPP-RA_RH_NEE_TimeSeries_split_"+RegionName+str(Numm)+"withGOSAT-GFED.png", dpi=400,bbox_inches='tight')

#S10 FLUXNET Plots, see extra script
