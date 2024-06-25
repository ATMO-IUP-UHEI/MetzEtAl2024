#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

DS0 = xr.open_dataset("/LandCover/MODIS_MCD12Q1/MCD12C1.A2015001.061.2022166123617_Majority_Land_Cover_Type_1.nc")
DS = xr.Dataset(dict(landcover=(['lat','lon'],DS0.__xarray_dataarray_variable__.values),lat=(['lat'],np.array(range(1800,-1800,-1))/20-0.025),lon=(['lon'],np.array(range(-3600,3600))/20+0.025)))
DS = DS.set_coords('lat')
DS = DS.set_coords('lon')
DS = DS.sortby('lat')

DSA = DS.sel(lat=slice(-40,-10),lon=slice(5,55))

classes1 = { 0: 'water',
           1: 'evergreen needleleaf forest',
           2: 'evergreen broadleaf forest',
           3:'deciduous needleleaf forest',
           4:'deciduous broadleaf forest',
           5:'mixed forests',
           6:'closed shrubland',
           7:'open shrublands',
           8:'woody savannas',
           9:'savannas',
           10:'grasslands',
           11:'permanent wetlands',
           12:'croplands',
           13:'urban and built-up',
           14:'cropland/natural vegetation mosaic',
           15:'snow and ice',
           16:'barren or sparsely vegetated'}

classes2 = { 0: 'water',
           1: 'evergreen needleleaf forest',
           2: 'evergreen broadleaf forest',
           3:'deciduous needleleaf forest',
           4:'deciduous broadleaf forest',
           5:'mixed forests',
           6:'closed shrubland',
           7:'open shrublands',
           8:'woody savannas',
           9:'savannas',
           10:'grasslands',
           11:'',
           12:'croplands',
           13:'urban and built-up',
           14:'',
           15:'snow and ice',
           16:''}

classes3 = { 0: 'water',
           1: 'grasses/cereal', 
           2: 'shrubs',
           3:'broadleaf crops', 
           4:'savannas',
           5: 'evergreen needleleaf forest',
           6: 'evergreen broadleaf forest',
           7:'deciduous needleleaf forest',
           8:'deciduous broadleaf forest',
           9:'unvegetared',
           10:'urban'}

colorlist1 = ['#0032C8','#58481F','#009900','#70663E',
             '#4E751F','green','goldenrod','darkkhaki',
             'darkorange','orange','gold','lightblue',
             'peru','darkred','peachpuff','white','darkgrey']

colorlist3 = ['#0032C8','gold','goldenrod','peru',
              'orange','#4E751F','#009900','#70663E',
             '#58481F','grey','darkred']
text = ''
for i in range(0,11):
    text = text + classes3[10-i] +' \n'

fig, ax = plt.subplots(figsize = (6,4))    
DSA.landcover.plot(levels=np.array(range(0,17))+0.5,colors=colorlist1)
#plt.text(70, -42, text,fontsize=11)
plt.savefig('LandCoverMap_type1_Lat-10.png',dpi=300, bbox_inches='tight')