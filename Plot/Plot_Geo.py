# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:49:03 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk
"""

'''
Generates the maps of Europe with average impacts scores used Figure3.

To ensure uniqueness of the files, each geodatframe file is given a unique ID 
corresponding to the time of production during the simulation. 
This ensures that each simulation has its own file if several simulations are performed.

For this reason, you have to fill in the name of the file for you  simulation
'''

# CHANGE THE NAME OF THE FILES.

namefile_gdfpoints="..\Outputs-Multi\gdfpoints_outputs_5_19_058642.pkl"


namefile_geodataframe="..\Outputs-Multi\geodataframe_output_5_19_058642.pkl"


''' Execute the whole script to produce the maps.''' 




import pickle
import os


currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)



import pandas as pd
import decimal
from random import *
from math import*
import csv
import copy
import numpy as np
import random


import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
#import pysal as ps
from shapely.geometry import Polygon, mapping, Point


import pickle





# Define function which imports pickle objects
def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    
        

"""Import geodataframes"""



# Multidimensional  Sampling


gdf_points_multi=importpickle("..\Outputs-Multi\gdfpoints_outputs_5_19_058642.pkl")


gdf_points_multi=importpickle(namefile_gdfpoints)


geodataframe_output_multi=importpickle("..\Outputs-Multi\geodataframe_output_5_19_058642.pkl")
geodataframe_output_multi=importpickle("..\Outputs-Multi\geodataframe_output_5_19_058642.pkl")



"""Plot"""








#GWP


#Multi

f, ax = plt.subplots(1,figsize=(8, 8),dpi=500)
geodataframe_output_multi.plot(axes=ax,alpha=1,column='Mean_GWP100',legend=True)
gdf_points_multi.plot(axes=ax,alpha=1,markersize=15,color="red")  

plt.xlabel('longitude (°)',size=17)
plt.ylabel('latitude (°)',size=17)
plt.legend()
#plt.tick_params(labelsize=20)
cb_ax = f.axes[1] 
cb_ax.tick_params(labelsize=15)
plt.show()

f.savefig('Geo_GWP_Fig3.jpeg',dpi=600)

plt.show()



#TETP




#Multi

f, ax = plt.subplots(1,figsize=(8, 8),dpi=500)
geodataframe_output_multi.plot(axes=ax,alpha=1,column='Mean_TETPinf',legend=True)
gdf_points_multi.plot(axes=ax,alpha=1,markersize=15,color="red")  
plt.xlabel('longitude (°)',size=17)
plt.ylabel('latitude (°)',size=17)
plt.legend()
#plt.tick_params(labelsize=20)
cb_ax = f.axes[1] 
cb_ax.tick_params(labelsize=15)
plt.show()

f.savefig('GEO_TETP_Fig3.jpeg',dpi=600)





#WDP




#Multi

f, ax = plt.subplots(1,figsize=(8, 8),dpi=500)
geodataframe_output_multi.plot(axes=ax,alpha=1,column='Mean_WDP',legend=True)
gdf_points_multi.plot(axes=ax,alpha=1,markersize=15,color="red")  
plt.xlabel('longitude (°)',size=17)
plt.ylabel('latitude (°)',size=17)
plt.legend()
#plt.tick_params(labelsize=20)
cb_ax = f.axes[1] 
cb_ax.tick_params(labelsize=15)
plt.show()

f.savefig('GEO_WDP_Fig3.jpeg',dpi=600)




# FEP




#Multi

f, ax = plt.subplots(1,figsize=(8, 8),dpi=500)
geodataframe_output_multi.plot(axes=ax,alpha=1,column='Mean_FEP',legend=True)
geodataframe_output_multi.plot(axes=ax,alpha=1,markersize=15,color="red")  

plt.xlabel('longitude (°)',size=17)
plt.ylabel('latitude (°)',size=17)
plt.legend()
#plt.tick_params(labelsize=20)
cb_ax = f.axes[1] 
cb_ax.tick_params(labelsize=15)
plt.show()

f.savefig('GEO_FEP_Fig3.jpeg',dpi=600)










