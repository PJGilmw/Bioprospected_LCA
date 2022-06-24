# -*- coding: utf-8 -*-
"""
Created on Sat May 22 01:11:59 2021

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk
"""


'''
Demonstrates the behavior of the module estimating the solar power received by a given PBR geometry.

Execute the block on the influence of azimuth to reproduce the figure in SI I.

# API source : https://ec.europa.eu/jrc/en/PVGIS/docs/noninteractive

'''



import requests
import matplotlib.pyplot as plt
import pandas as pd
import decimal
import random
import numpy as np
import os

# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)

import Retrieving_solar_and_climatic_data_1 as solardata
import Functions_for_physical_and_biological_calculations_1 as functions









#Influence of the length of the PBR

# Assuming the PBR unit is a rectangle with width=length/3

df_resultslength= pd.DataFrame(columns=['Upper','Lower','Average'])

length_values=range(1,50)

for length in length_values:
    
    datacollection=solardata.Qreceived_bym2PBR_month(43.695, 1.922, 3, 90, 1.5, 0.03, 0.01, 0.2, length, length/3)

    averageday_persquaremeter=sum(datacollection[0]['Average'])
    
    upperday_persquaremeter=sum(datacollection[0]['Upper'])
    
    lowerday_persquaremeter=sum(datacollection[0]['Lower'])

    df_resultslength.loc[length]=[upperday_persquaremeter,lowerday_persquaremeter,averageday_persquaremeter]

plt.plot(length_values,df_resultslength)

plt.ylabel('Wh collected per day')

plt.xlabel('Length of the PBR (width=length/3)')

plt.savefig('../Plot/Sensitivity_length_Qm2',dpi=600)





###Influence of the Tube diameter



df_resultsdiameter=pd.DataFrame(columns=['Upper','Lower','Average'])

diameter_values=range(1,11)

diameter_values=[a/100 for a in diameter_values]


for diameter in diameter_values:
    
    
    datacollection=solardata.Qreceived_bym2PBR_month(43.695, 1.922, 3, 90, 1.5, diameter, 0.01, 0.2, 25, 25/3)

    averageday_persquaremeter=sum(datacollection[0]['Average'])
    
    upperday_persquaremeter=sum(datacollection[0]['Upper'])
    
    lowerday_persquaremeter=sum(datacollection[0]['Lower'])

    df_resultsdiameter.loc[diameter]=[upperday_persquaremeter,lowerday_persquaremeter,averageday_persquaremeter]

plt.plot(diameter_values,df_resultsdiameter)
plt.ylabel('Wh collected per day')
plt.xlabel('Tube diameter in m')
plt.savefig('../Plot/Sensitivity_Tube_diameter',dpi=600)




###Influence of the gap between tubes


df_resultsgap=pd.DataFrame(columns=['Upper','Lower','Average'])

gap_values=range(1,30)

gap_values=[a/100 for a in gap_values]

for gap in gap_values:
    
    datacollection=solardata.Qreceived_bym2PBR_month(43.695, 1.922, 3, 90, 1.5, 0.03,gap, 0.2, 25, 25/3)

    averageday_persquaremeter=sum(datacollection[0]['Average'])
    
    upperday_persquaremeter=sum(datacollection[0]['Upper'])
    
    lowerday_persquaremeter=sum(datacollection[0]['Lower'])

    df_resultsgap.loc[gap]=[upperday_persquaremeter,lowerday_persquaremeter,averageday_persquaremeter]

plt.plot(gap_values,df_resultsgap)
plt.ylabel('Wh collected per day')
plt.xlabel('Gap between rows in m')
plt.savefig('../Plot/Sensitivity_horizontal_distance',dpi=600)



###Influence of the horizontal distance between stacks


df_resultshori_dist=pd.DataFrame(columns=['Upper','Lower','Average'])

hori_dist_values=range(30,300,10)

hori_dist_values=[a/100 for a in hori_dist_values]

for hori_dist in hori_dist_values:
    
    datacollection=solardata.Qreceived_bym2PBR_month(43.695, 1.922, 3, 90, 1.5, 0.60,0.01, hori_dist, 25, 25/3)

    averageday_persquaremeter=sum(datacollection[0]['Average'])
    
    upperday_persquaremeter=sum(datacollection[0]['Upper'])
    
    lowerday_persquaremeter=sum(datacollection[0]['Lower'])

    df_resultshori_dist.loc[hori_dist]=[upperday_persquaremeter,lowerday_persquaremeter,averageday_persquaremeter]

plt.plot(hori_dist_values,df_resultshori_dist)
plt.ylabel('Wh collected per day')
plt.xlabel('Horizontal distance between stacks')
plt.savefig('../Plot/Sensitivity_horizontal_distance',dpi=600)





### Azimuth frontal influence

# Qreceived_bym2PBR_month(lat, 
#                             long,
#                             month,
#                             azimuthfrontal,
#                             height,
#                             tubediameter,
#                             gapbetweentubes,
#                             horizontaldistance,
#                             length_of_PBRunit,
#                             width_of_PBR_unit)

df_resultsazimuth=pd.DataFrame(columns=['Upper','Lower','Average'])
azi_values=[0,45,90,135,180,-45,-90,-135]

for azi in azi_values:
    
    datacollection=solardata.Qreceived_bym2PBR_month(43.695, 1.922, 3, azi, 1.5, 0.1,0.01, 0.2, 25, 25/3)

    averageday_persquaremeter=sum(datacollection[0]['Average'])
    upperday_persquaremeter=sum(datacollection[0]['Upper'])
    lowerday_persquaremeter=sum(datacollection[0]['Lower'])

    df_resultsazimuth.loc[azi]=[upperday_persquaremeter,lowerday_persquaremeter,averageday_persquaremeter]

# plt.plot(azi_values, df_resultsazimuth)

df_resultsazimuth.plot()
plt.ylabel('Wh collected per day')
plt.xlabel('Azimuth frontal')
plt.legend(fontsize=8,loc='best')

plt.savefig('../Plot/Sensitivity_Azimuth',dpi=600)
