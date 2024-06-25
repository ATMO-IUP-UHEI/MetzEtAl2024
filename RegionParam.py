#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# part of the Code for Metz et al. 2024
# Author: Eva-Marie Metz
Function to get Region Name and Coordinates
"""

def getRegion(Numm):
    '''
    # Function to get Region Name and Coordinates based on identification number
    # arguments:
    #           Numm: identification number of region
    # returns:
    #           RegionName: name of region as string
    #           Long_min, Long_max, Lat_min, Lat_max: coordinates of bounding box
    '''

    #Regions not based on TRANSCOM Regions
    region_names = ['world',               
                ]


    #Region =   [0  ]
    Long_minl = [-180]
    Long_maxl = [180 ]
    Lat_minl =  [-90  ]
    Lat_maxl =  [90 ]
    



    #Regions based on TRANSCOM Regions
    Transcom = ['SAT','SA','AU','SATr','SATr_SE','I_P','SAT','SAT','SAT','SAT', #0-9
                'SAT','SAT','SA','SA','SA','SA','SA','AU','AU','AU', #10-19
                'AU','AU','AU','AU','AU','AU','AU','AU','AU','AU', #20-29
                'SAT','SAT','SAT','SAT','SA','SA','SA','Af','SAT_SATr','TA',#30-39
                'EU','NAS','As','NAN','EB','AU','AU','AU','AU','AU', #40-49
                'AU','AU','AU','AU','Af','Af','Af','Af','Af','Af',#50-59
                'Af','Af','Af','Af','SA','SA', 'SA', 'SAT','Af','Af', #60-69 
                'Af','Af','Af','NAB','NAT','SATr','SAT','NA','SA','EB', #70-79
                'ET','TA','AU','EU','Af'#80-89
                ]


    coordt = [[-56,0,-84,-32],
          [-36,0,7,57],
          [-48,-12,114,180], 
          [-12,15,-84,-32], 
          [-12,0,-84,-32], #4
          [-12,6,90,160],
          [-56,0,-80,-60], #Longitudinal investigations
          [-56,0,-75,-55],
          [-56,0,-70,-50],
          [-56,0,-65,-45],
          [-56,0,-60,-40], #10
          [-56,0,-55,-35],
          [-36,0,7,30],
          [-36,0,15,35],
          [-36,0,20,40],
          [-36,0,25,45], #15
          [-36,0,30,57],
          [-48,-12,114,135],
          [-48,-12,120,140],
          [-48,-12,125,145],
          [-48,-12,130,150], #20
          [-48,-12,135,160],
          [-48,-12,160,180],
          [-48,-12,114,125],
          [-48,-12,125,135],
          [-48,-12,135,145], #25
          [-48,-12,145,155], #------------------------------------
          [-24,-12,114,180], #Latitudinal investigations
          [-36,-24,114,180],
          [-48,-36,114,180],
          [-15,0,-84,-32],   #30
          [-25,-15,-84,-32],
          [-40,-25,-84,-32],
          [-60,-40,-84,-32],
          [-10,0,7,57],
          [-20,-10,7,57],    #35
          [-36,-20,7,57],
          [-36,40,-20,50],      
          [-56,15,-84,-32],    
          [-12,30,75,160],      
          [30,85,-30,180],  #40 
          [0,49,-180,-10],      
          [0,60,30,150],      
          [49,90,-180,-10],      
          [40,90,60,180] ,      
          [-30,-22,125,145],  #45    
          [-22,-10,114,180], 
          [-29,-22,114,180],
          [-50,-29,114,180],  
          [-50,-10,112,180], #49 
          [-50,-10,112,180], #50
          [-50,-10,112,180], #51 
          [-50,-10,112,180], #52 
          [-50,-10,112,180], #53 
          [10,35,-19,57], #54 
          [-10,10,-19,57], #55 
          [-35,-10,7,57], #56 
          [10,38,-19,57], #57 
          [-10,1,7,57], #58 
          [-17,-10,7,57], #59 
          [-25,-17,7,57], #60 
          [-35,-25,7,57], #61 
          [-35,-10,7,32], #62 
          [-35,-10,32,57], #63 
          [-35,1,7,32], #64 
          [-35,1,32,57], #65 
          [-36,2,7,57],#66 
          [-60,0,-90,-25],#67 
          [-35,-10,7,42], #68 
          [-20,-10,7,57], #69 
          [-35,-20,7,57], #70 
          [-35,-17,7,27], #71 
          [-35,-17,7,42], #71 
          [45,85,-180,-50], #73 
          [10,60,-140,-60], #74 
          [-20,30, -100, -40], #75 
          [-60,0,-80,-30], #76 
          [-1,40,-20,51], #77 
          [-35,2,8,57], #78 
          [30,90,-180,180], #79 
          [0,60,20,150], #80 
          [-20,30,75,180], #81 
          [-50,-10,110,180], #82 
          [33,75,-12,65], #83
          [-17,-10,7,42], #84
          ]

    if Numm <700:
        Long_min = Long_minl[Numm]
        Long_max = Long_maxl[Numm]
        Lat_min = Lat_minl[Numm]
        Lat_max = Lat_maxl[Numm]
    
        RegionName = region_names[Numm]

    elif Numm < 900:
        Long_min = coordt[Numm - 700][2]
        Long_max = coordt[Numm - 700][3]
        Lat_min = coordt[Numm - 700][0]
        Lat_max = coordt[Numm - 700][1]
    
        RegionName = Transcom[Numm - 700]
    else:
        Long_min = coordt[Numm - 900][2]
        Long_max = coordt[Numm - 900][3]
        Lat_min = coordt[Numm - 900][0]
        Lat_max = coordt[Numm - 900][1]
    
        RegionName = Transcom[Numm - 900]
    
    return RegionName, Long_min, Long_max, Lat_min, Lat_max

