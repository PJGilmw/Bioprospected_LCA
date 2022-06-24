# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 23:43:11 2022

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script to execute to conduct the Multi-dimensional sampling and associated calculations.

'''

""" Choose Simulation parameters"""
""" The current values are the ones used in the article.

################
#####WARNING### :
    The code was built to be run efficiently in parallel on a server/computer with numerous cores.
This script was run on a remote server with 64 cores, 256 GB RAM and it took around 7 days to produce the results with the chosen samples' sizes.
A laptop with fewer cores will take much more time.

In addition, the results will occupy around 25 GB on the disk.

To test the code, one can use lower sample sizes. 

An internet connection is necessary as the code downloads climatic data from an API.

"""

Size_sample = 256 # Size of the Sobol sample for the strain-compound pairs 
size_montecarlo= 150 # Size of the Monte Carlo Sample for the uncertain techno-operational parameters
cultivationperiod=[4,5,6,7,8,9] # 1=January, 2= February etc.
fraction_max_yield=0.3 # Fraction of maximum yiled achieved by the cultivated strains
senstype="SOBOL"
latsec= 2 # Latitude gap (°) between 2 random locations (yield 28 locations)

""" Execute the script."""





import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint

import os
import sys
# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)

import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy
import numpy as np
import random


import bw2data
import bw2io
from bw2data.parameters import *
import brightway2 as bw


from SALib.test_functions import Ishigami
import math
from SALib.sample import saltelli
from SALib.sample import fast_sampler
from SALib.analyze import sobol
from SALib.analyze import fast
import SALib

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d

import shapefile as shp
import geopandas as gpd
#import pysal as ps
from shapely.geometry import Polygon, mapping, Point


import Cultivation_simul_2nd as cultsimul
import Functions_for_physical_and_biological_calculations_2nd as functions
import Main_simulations_functions_Multi_2nd as mainfunc
import Map_2nd as geo
import Retrieving_solar_and_climatic_data_2nd as solardata

import ray




"""Model Parameters"""





# Primary Parameters dictionnaries

# Description of the parameters are given in the appendix.
# Here values can be changed for parameters with unique values (no distributions)
# Values with distributions will be overwritten.

#Biological parameters


Biodict = {'rhoalgae': 1070, # kg.m-3
           'lipid_af_dw': 0.3,   # .
           'ash_dw': 0.05,  # .
           'MJ_kglip': 36.3,  # MJ.kg-1
           'MJ_kgcarb': 17.3,  # MJ.kg-1
           'MJ_kgprot': 23.9,  # MJ.kg-1
           'PAR': 0.45,  # .
           'losspigmentantenna': 0.21,  # .
           'quantumyield': 8,  # molphotons.moloxygen-1
           'lossvoltagejump': 1-1267/1384,  # .
           'losstoATPNADPH': 1-590/1267,  # .
           'losstohexose': 1-469/590,  # .
           'lossrespiration': 0.20,  # .
           'bioact_fraction_molec': 0.1,  # .
           'random_no3': 1,   # .
           'Topt': 25,  # °C
           'T_plateau': 10,  # °C
           'dcell': 5*10**-6,  # m
           'incorporation_rate': 0.10,  # .
           'nutrient_utilisation': 0.9, # .
           'co2_utilisation': 0.85,  # .
           'phospholipid_fraction': 0.20,  # .

            'random_market_subst': 0.1,  # .
            
            # A value of 1 for Thermoregulation at night
            
            'random_night_monitoring': 0,
            
                         
            'random_bioclass': 0, 
            
            'fertiliser_substitutability_AD': 0.8,
                                 
            'carbon_degradibilidy_AD':0.7,
           }

# Physical parameters

Physicdict = {'hconv': 6.6189,  # W.m−2 .K−1
              'Cp': 4.186,  # kJ.(kg.K)-1
              'rhomedium': 1000,  # kg.m-3
              'rhowater': 1000,  # kg.m-3
              'Cw': 2256,     # kJ.kg-1
              'CH4_LHV':50}  # MJ.kg-1


#Geographic parameters

#  Granada 37.189, -3.572
#  Aalborg 57.109, 10.193


Locationdict = {'lat': 37.189, #°
                'long': -3.572, #°
                'depth_well': 10,  # m
                'azimuthfrontal': 90} # °


#Techno-operational parameters

Tech_opdict = {'height': 1.5,  # m
                   'tubediameter': 0.03,  # m
                   'gapbetweentubes': 0.01,  # m
                   'horizontaldistance': 0.2,  # m
                   'length_of_PBRunit': 30,    # m
                   'width_of_PBR_unit': 30,   # m
                   'biomassconcentration': 1.4,  # kg.m-3
                   'flowrate': 0.38,     # m.s-1
                   'centrifugation_efficiency': 0.98,  # .
                   'pumpefficiency': 0.90,  # .
                   'slurry_concentration': 0.15,  # .
                   'water_after_drying': 0.05,  # gwater.g dbio-1
                   'recyclingrateaftercentrifuge': 0.3, # gwater.g dbio-1
                   'rhosuspension': 1015,  # kg.m-3 
                   'roughness': 0.0000015,  # m
                   'cleaningvolumeVSfacilityvolume': 4,  # .
                   'concentration_hypo': 2*10**-3,  # kg.m-3
                   'concentration_hydro': 30,  # kg.m-3
                   'boilerefficiency': 0.75,  # .
                   'glass_life_expectancy': 50, # years
                   

                   'extraction': 'yes',
                   
                   'heat_pump': 'yes',
                   
                   'COP': 3,
                   
                   'heat_input_per_kg_biomass_AD':0.68, # kWh.kg dw -1 
                                 
                   "elec_input_per_kg_biomass_AD":0.108 # kWh.kg dw -1


                   }




# Parameters distribution directories 

# Exact same parameters dictionnaries the normal ones but instead of one value,
# each parameter is assigned a list containing  :
   # [Distribution,min,max,mode,sd]

# Distribution :
#   - 'unique' if no distribution. The value indicated indicated in
#     the normal dictionnary will be considered.
#   - 'unif' for uniform, uses min and max
#   - 'triang, uses mim max and mode with mode as a fracion of max-min

#

# Biological parameters



Biodict_distributions = {'rhoalgae': ['unif', [0, 1020, 1250, 0, 0]], 
                           
                         'lipid_af_dw': ['triang', [0, 0.1, 0.7, 0.2, 0]],    # see dist_lip_williams
                         
                         'ash_dw': ['unif', [0, 0.01, 0.10, 0.10, 0]],
                         
                         'MJ_kglip': ['unique', [36.3, 0, 0, 0, 0]],
                       
                         'MJ_kgcarb': ['unique', [17.3, 0, 0, 0, 0]],
                        
                         'MJ_kgprot': ['unique', [23.9, 0, 0, 0, 0]],
                         
                         'PAR': ['unique', [0.45, 0, 0, 0, 0]],
                        
                         'losspigmentantenna': ['unif', [0, 0.16, 0.24, 0, 0]],
                        
                         'quantumyield': ['unif', [0, 8, 11, 0, 0]],
                         
                         'lossvoltagejump': ['unique', [1-1267/1384, 0, 0, 0, 0]],
                       
                         'losstoATPNADPH': ['unique', [1-590/1267, 0, 0, 0, 0]],
                 
                         'losstohexose': ['unique', [1-469/590, 0, 0, 0, 0]],
                      
                         'lossrespiration': ['unif', [0, 0.16, 0.24, 0, 0]],
                    
                         'bioact_fraction_molec': ['unif', [0, 0.01, 1, 0, 0]],
                         
                         'random_no3': ['unif', [0, 0, 1, 0, 0]],
                         
                         'Topt': ['unif', [0, 15, 35, 25, 0]],
                         
                         'T_plateau': ['unif', [0, 5, 10, 0, 0]],
                      
                         'dcell': ['unif', [0, 0.000001, 0.00005, 0.000005, 0]],
                
                         'incorporation_rate': ['unique', [0, 0.01, 0.15, 0, 0]],

                         'nutrient_utilisation': ['unif', [0, 0.75, 0.90, 0, 0]],
             
                         'co2_utilisation': ['triang', [0, 0.5, 0.95, 0.9, 0]],
            
                         'phospholipid_fraction': ['unif', [0, 0.1, 0.6, 0, 0]],

                         'random_market_subst': ['unif', [0, 0, 1, 0, 0]],
                         
                         'random_night_monitoring': ['unif', [0, 0, 1, 0, 0]],
                         
                         'random_bioclass': ['unif', [0, 0, 1, 0, 0]], 

                         'fertiliser_substitutability_AD': ['unif', [0.8, 0.6, 0.9, 0, 0]],
                                 
                         'carbon_degradibilidy_AD': ['unif', [0.7, 0.5, 0.6, 0, 0]],
                         }


# Physical parameters

Physicdict_distributions = {
    'hconv': ['unif', [0, 5, 10, 0, 0]],
    
    'Cp': ['unique', [4.186, 0, 0, 0, 0]],
     
    'rhomedium': ['unique', [1000, 0, 0, 0, 0]],
    
    'rhowater': ['unique', [1000, 0, 0, 0, 0]],
    
    'Cw': ['unique', [2256, 0, 0, 0, 0]], #kJ.kg-1
    
    'CH4_LHV':['unique', [50, 0, 0, 0, 0]]}  # MJ.kg-1



Locationdict_distributions = {'lat': ['unique', [43.695, 0, 0, 0, 0]],
                              
                              'long': ['unique', [1.922, 0, 0, 0, 0]],
                                                            
                              'depth_well': ['unif', [0, 5, 25, 0, 0]],
                              
                              'azimuthfrontal': ['unique', [90, 0, 0, 0, 0]]}




Tech_opdict_distributions = {'height': ['unique', [1.5, 0, 0, 0, 0]],   
                              
                                 'tubediameter': ['unif', [0, 0.03, 0.1, 0, 0]],
                                 
                                 'gapbetweentubes': ['unif', [0, 0.01, 0.1, 0, 0]],
                                 
                                 'horizontaldistance': ['unif', [0, 0.2, 0.8, 0, 0]],
                                 
                                 'length_of_PBRunit': ['unique', [20, 10, 35, 0, 0]], 
                                 
                                 'width_of_PBR_unit': ['unique', [20, 3, 10, 0, 0]], 
                                 
                                 'biomassconcentration': ['unif', [0, 1, 7, 0, 0]], 
                                 
                                 'flowrate': ['unif', [0, 0.2, 1, 0, 0]],
                                 
                                 'centrifugation_efficiency': ['unique', [0.98, 0, 0, 0, 0]], 
                                 
                                 'pumpefficiency': ['unique', [0.9, 0, 0, 0, 0]], 
                                                                  
                                 'slurry_concentration': ['unique', [0, 0.10, 0.20, 0, 0]], 
                                 
                                 'water_after_drying': ['unique', [0, 0.02, 0.07, 0, 0]], 
                                 
                                 'recyclingrateaftercentrifuge': ['unique', [0, 0.2, 0.4, 0, 0]],
                                 
                                 'rhosuspension': ['unif', [0, 1000, 1200, 0, 0]], 
                                 
                                 'roughness': ['unique', [0.0000015, 0, 0, 0, 0]], 
                                 
                                 'cleaningvolumeVSfacilityvolume': ['unique', [0, 2, 6, 0, 0]], 
                                 
                                 'concentration_hypo': ['unique', [2*10**-3, 0, 0, 0, 0]], 
                                 
                                 'concentration_hydro': ['unique', [30, 0, 0, 0, 0]], 
                                 
                                 'glass_life_expectancy': ['unique', [50, 0, 0, 0, 0]], 

                                 'boilerefficiency': ['unique', [0.75, 0, 0, 0, 0]],  
                                                                  
                                 'extraction': ['binary', ['yes', 0, 0, 0, 0]],
                                 
                                 'heat_pump': ['binary', ['yes', 0, 0, 0, 0]],
                   
                                 'COP': ['unique', [3, 2, 6, 0, 0]],
                                 

                                 
                                 'heat_input_per_kg_biomass_AD': ['unique', [0.68, 0, 0, 0, 0]], # kWh.kg dw -1
                                 
                                 'elec_input_per_kg_biomass_AD': ['unique', [0.108, 0, 0, 0, 0]] # kWh.kg dw -1
                                 
                                 
                                 
                                 }











# Managing Brightway projects and databases

# Loading the right project

bw.projects.set_current("Microalgae_Sim")



# Loading Ecoinvent
Ecoinvent = bw.Database('ecoinvent 3.6 conseq')


# Loading foreground database

MICAH = bw.Database('Microalgae_foreground_2')

# Loading biosphere

biosph = bw.Database('biosphere3')



# Loading necessary activites from the database to use them as inputs to the functions

for act in MICAH:
    # The original wastewater treatment activity
    if act['name'] == 'treatment of wastewater, average, capacity 1E9l/year PBR':
        wastewater = MICAH.get(act['code'])
        
for act in MICAH:
    # The original Molecule production activity
    if act['name'] == 'Molecule production PBR':
        # print('oik')
        Molprod = MICAH.get(act['code'])   
        
        
# Names of IA methods

methods_selected = ([
     ('ReCiPe Midpoint (H) V1.13', 'terrestrial ecotoxicity', 'TETPinf'),
     ('ReCiPe Midpoint (H) V1.13', 'climate change', 'GWP100'),
     ('ReCiPe Midpoint (H) V1.13', 'freshwater eutrophication', 'FEP'),
     ('ReCiPe Midpoint (H) V1.13', 'water depletion', 'WDP')])






# Accessory functions

def createFolder(directory):
    '''Creates a folder/directory with the path given as input'''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def export_pickle(var, name_var):
    '''Saves a pickle in the working driectory and
    saves the object in input across python sessions'''

    with open(name_var+'.pkl', 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

def export_pickle_2(var, name_var, namefolder_in_root):
    '''Saves a pickle in the working directory and
    saves the object in input across python sessions'''

    path_object = "../"+namefolder_in_root+"/"+name_var+".pkl"
    with open(path_object, 'wb') as pickle_file:
        pickle.dump(var, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)



def importpickle(path):
    with open(path, 'rb') as pickle_load:
        obj = pickle.load(pickle_load)
    return obj    



# Function to modify the electricty impacts according to the national mix

# Necessary
list_processes_electricity = []
for exc in list(Molprod.exchanges()):
#    print(exc)
        if exc['type']!='production':
            
            exchange1 = MICAH.get(exc['input'][1])  
            
            exchange1_name=exchange1['name']
            print(exchange1)
            for exc in list(exchange1.exchanges()):
                
               # "MJ heat biogas AD PBR" has an input of the modified anaerobic digestion actvitiy 
               #from the foreground database, it is not electricity and it can't be found with Ecoinvent.get() 
                
                if exc['type']=='technosphere' and exchange1_name!="MJ heat biogas AD PBR":
                    
                    act_background = Ecoinvent.get(exc['input'][1])  
                    
                    name_background = act_background['name']
                    
                    #print(name_background)
                    if 'electricity' in name_background:
                        print('ok')
                        list_processes_electricity.append(exchange1_name)
          
 

def modif_impact_elec_national_mix(dict_mono_technosphere_lcas, list_processes_electricity, listimpactforamix):
    """Function which modifies the impact of the foreground activites according to the impact of a mix """

    for input_ in dict_mono_technosphere_lcas:
        
        #print(input_)
        if input_ in list_processes_electricity:
            
            dict_mono_technosphere_lcas[input_] = listimpactforamix
    
    return dict_mono_technosphere_lcas    









''' Initialization LCA'''


list_cfs= [bw.Method((meth)).load() for meth in methods_selected]


# Calculating the Impacts for 1 unit of each of the technosphere inputs 
# to the molecule production

list_foreground_technosphere_inputs_FU = []

list_foreground_technosphere_inputs_names = []

for exc in list(Molprod.exchanges()):

        if exc['type']!='production':
            
            exchange1 = MICAH.get(exc['input'][1])  
            
            name_exchange = exchange1['name']  
            
            list_foreground_technosphere_inputs_names.append(name_exchange)
        
            list_foreground_technosphere_inputs_FU.append({exchange1 : 1})
    

my_calculation_setup = {'inv': list_foreground_technosphere_inputs_FU, 'ia': methods_selected}

bw.calculation_setups['mono_technosphere_inputs'] = my_calculation_setup



# Calculating the impacts for all chosen methods

mlca_techno = bw.MultiLCA('mono_technosphere_inputs')  

res = mlca_techno.results

dict_mono_technosphere_lcas = { name : results for (name, results) in zip(list_foreground_technosphere_inputs_names,res) }




# Wastewater treatment activity is broken down to technoshpere inputs and biosphere outputs
# Technosphere inputs per cubic meter of wastewater constant for any microalga
# Biosphere outputs depend on the composition of the waste water and are then 
# calculated after the cultivation simulation



#Technosphere of wastewater treatment


list_wastewater_technosphere_inputs_FU = []

list_wastewater_technosphere_inputs_names = []

for exc in list(wastewater.exchanges()): # original activity

        if exc['type']=='technosphere':
            
            exchange1 = Ecoinvent.get(exc['input'][1])  # Full name
            
            name_exchange = exchange1['name']  # Name
            
            list_wastewater_technosphere_inputs_names.append(name_exchange)
        
            list_wastewater_technosphere_inputs_FU.append({exchange1 : exc['amount']})
    


my_calculation_setup = {'inv': list_wastewater_technosphere_inputs_FU, 'ia': methods_selected}

bw.calculation_setups['mono_technosphere_wastewater_inputs'] = my_calculation_setup

mlca_waste_water_techno = bw.MultiLCA('mono_technosphere_wastewater_inputs')  


# Impact for the technosphere inputs  for 1 cubic meter of  to waste water 

res_wastewater_technospere = mlca_waste_water_techno.results


dict_mono_technosphere_wastewater_lcas = { name : results for (name, results) in zip(list_wastewater_technosphere_inputs_names,res_wastewater_technospere) }


# Sum the impacts for all technosphere inputs to the treatment of 1 cubic meter
                                                                                            
list_sum_impacts_technosphere_waste_water = []

for meth_index in range(len(methods_selected)):
    
    sum_impact = sum([dict_mono_technosphere_wastewater_lcas[flow][meth_index] for flow in dict_mono_technosphere_wastewater_lcas ])
    
    list_sum_impacts_technosphere_waste_water.append(sum_impact)
    
    
# Update the dictionnary with the impacts associated to the technosphere inputs to 1 cubic meter of wastewater treatment.  

dict_mono_technosphere_lcas['Wastewater treatment PBR'] = list_sum_impacts_technosphere_waste_water         






# Wastewater treatment activity is broken down to technoshpere inputs and biosphere outputs
# Technosphere inputs per cubic meter of wastewater constant for any microalga
# Biosphere outputs depend on the composition of the waste water and are then 
# calculated after the cultivation simulation



#Technosphere of wastewater treatment


list_wastewater_technosphere_inputs_FU = []

list_wastewater_technosphere_inputs_names = []

for exc in list(wastewater.exchanges()): # original activity

        if exc['type']=='technosphere':
            
            exchange1 = Ecoinvent.get(exc['input'][1])  # Full name
            
            name_exchange = exchange1['name']  # Name
            
            list_wastewater_technosphere_inputs_names.append(name_exchange)
        
            list_wastewater_technosphere_inputs_FU.append({exchange1 : exc['amount']})
    


my_calculation_setup = {'inv': list_wastewater_technosphere_inputs_FU, 'ia': methods_selected}

bw.calculation_setups['mono_technosphere_wastewater_inputs'] = my_calculation_setup

mlca_waste_water_techno = bw.MultiLCA('mono_technosphere_wastewater_inputs')  


# Impact for the technosphere inputs  for 1 cubic meter of  to waste water 

res_wastewater_technospere = mlca_waste_water_techno.results


dict_mono_technosphere_wastewater_lcas = { name : results for (name, results) in zip(list_wastewater_technosphere_inputs_names,res_wastewater_technospere) }


# Sum the impacts for all technosphere inputs to the treatment of 1 cubic meter
                                                                                            
list_sum_impacts_technosphere_waste_water = []

for meth_index in range(len(methods_selected)):
    
    sum_impact = sum([dict_mono_technosphere_wastewater_lcas[flow][meth_index] for flow in dict_mono_technosphere_wastewater_lcas ])
    
    list_sum_impacts_technosphere_waste_water.append(sum_impact)
    
    
# Update the dictionnary with the impacts associated to the technosphere inputs to 1 cubic meter of wastewater treatment.  

dict_mono_technosphere_lcas['Wastewater treatment PBR'] = list_sum_impacts_technosphere_waste_water         


 

# dict_mono_technosphere_lcas now contains the impact for 1 unit of each process of the LCI 
# except for the direct biosphere flows from wasterwater treatment which will be added for each LCI calculation.
                                               




# Initializing the LCI dictionnary

# Each flow has the same name of the corresponding activity 
# in the original database

LCIdict = {'market for ammonium sulfate, as N PBR': 0,  
           'market for calcium nitrate PBR': 0,
           'P source production PBR': 0,
           'Hydrogen peroxyde PBR': 0,
           'K source production PBR': 0,
           'Land PBR': 0,
           'Glass PBR': 0,
           'Hypochlorite PBR': 0,
           'Microalgae CO2 PBR': 0,
           'Heating kWh PBR': 0,
           'Cooling kWh PBR': 0,
           'Electricity centrifuge kWh PBR': 0,
           'Electricity mixing kWh PBR': 0,
           'Electricity pumping kWh PBR': 0,
           'Electricity drying kWh PBR': 0,
           'Electricity cell disruption kWh PBR': 0,
           'Electricity cell disruption kWh PBR': 0,
           'Electricity aeration kWh PBR':0,
           'Feed energy PBR': 0,
           'Feed protein PBR': 0,
           'LT Fishmeal PBR': 0,
           'Rapeseed oil PBR': 0,
           'Wheat PBR': 0,
           'Wheat gluten PBR': 0,
           'Fish oil PBR': 0,
           'Soyabean meal PBR': 0,
           'Poultry meal PBR': 0,
           'Hemoglobin meal PBR': 0,
           'Electricity cell disruption kWh PBR': 0,
           'Co solvent Extraction PBR': 0,
           'Extraction electricity kWh PBR': 0,
           'Mg source production PBR': 0,
           'Wastewater treatment PBR': 0,
           'Water Cleaning PBR': 0,
           'Water(Cultivation) PBR': 0,
           'CO2 direct emissions PBR': 0,
           "N fertiliser substitution AD PBR":0,
           "P fertiliser substitution AD PBR":0,
           "CO2 substitution biogas burning AD PBR":0,
           "K fertiliser substitution AD PBR":0,
           "Mg fertiliser substitution AD PBR":0,
           "MJ heat biogas AD PBR":0,
           "Electricity kWh AD PBR":0}


# Assigning each process to a a broader category for contribution analysis

# Names of categories

categories_contribution = ['Thermoregulation',
                           'Direct land occupation',
                           'Post harvest processing',
                           'Harvesting',
                           'Culture mixing',
                           'Water Pumping',
                           'Substitution by co-products',
                           'Use of glass',
                           'Nutrients consumption',
                           'CO2 consumption',
                           'Water consumption and treatment',
                           'Cleaning']


# List of processes to assign to categories (same order)

processes_in_categories = [['Cooling kWh PBR', 'Heating kWh PBR'],
                           ['Land PBR'],
                           ['Co solvent Extraction PBR', 'Extraction electricity kWh PBR',
                               'Electricity cell disruption kWh PBR', 'Electricity drying kWh PBR'],
                           ['Electricity centrifuge kWh PBR'],
                           ['Electricity mixing kWh PBR'],
                           ['Electricity pumping kWh PBR'],
                           ['Feed energy PBR', 'Feed protein PBR', 'Fish oil PBR', 'Hemoglobin meal PBR', 'LT Fishmeal PBR',
                            'Poultry meal PBR', 'Rapeseed oil PBR', 'Soyabean meal PBR', 'Wheat gluten PBR', 'Wheat PBR',
                            "N fertiliser substitution AD PBR","P fertiliser substitution AD PBR","CO2 substitution biogas burning AD PBR",
                            "K fertiliser substitution AD PBR","MJ heat biogas AD PBR","Electricity kWh AD PBR","Mg fertiliser substitution AD PBR"],
                           ['Glass PBR'],
                           ['market for ammonium sulfate, as N PBR', 'market for calcium nitrate PBR', 'P source production PBR',
                             'K source production PBR', ],
                           ['Microalgae CO2 PBR','CO2 direct emissions PBR'],
                           ['Water(Cultivation) PBR','Wastewater treatment PBR'],
                           ['Hypochlorite PBR', 'Hydrogen peroxyde PBR','Water Cleaning PBR']]



# Uploading the necessary csv files

# Fish feed composition

fishfeed_table = pd.read_csv("../Data/Feed_composition.csv", sep=";",
                             header=0, encoding='unicode_escape', engine='python')

# Cleaning

fishfeed_table = fishfeed_table[0:8][[
    'Ingredient', 'kg.kg feed-1', 'Lipid', 'Protein', 'Carb', 'Ash', 'Water']]

# Elemental composition of macronutrients

elemental_contents = pd.read_csv("../Data/elemental_contents.csv",
                                 sep=";",
                                 header=0,
                                 encoding='unicode_escape',
                                 engine='python')

# Cleaning

elemental_contents = elemental_contents.iloc[:, 1:]





""" Simulation Functions"""



def pre_download_climatic_data(latsec,methods_selected,biosph,MICAH,Ecoinvent,cultivationperiod ):
    """Function which generates the grid of locations"""
    
    # Generating the grid of points and create the associated geodataframes containing the impacts of the electricity mixes.    
    geodataframe_metrop, gdf_points = geo.geodataframe_initialize(latsec,methods_selected,biosph,MICAH,Ecoinvent )

    

    # We make sure to download all the climatic data for the grid:


    for country_index in range(geodataframe_metrop.shape[0]):
        
        list_points = geodataframe_metrop['random_points'].iloc[country_index]
    
        for point_index in range(len(list_points)):
        
        # corresponding rank (row) of thhis point in the geodataframe containing points
        
            
        
            point = list_points[point_index]
            #print(point)
            

            for month in cultivationperiod:
                
                # This will download the climatic data
        
                solardata.Qreceived_bym2PBR_month(round(point[1],3),
                                                        round(point[0],3),
                                                        month,
                                                        90,
                                                        1.5,
                                                        0.3,
                                                        0.05,
                                                        0.3,
                                                        30,
                                                        30)
    
    # Export the grid of locations
    
    
    x = datetime.datetime.now()
    
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    name_geo_countries ='geodataframe_input'+"_"+month+"_"+day+"_"+microsec
    
    name_geo_points ='gdfpoints_input'+"_"+month+"_"+day+"_"+microsec


    export_pickle_2(geodataframe_metrop, name_geo_countries, "Outputs-Multi")
    
    export_pickle_2(gdf_points, name_geo_points, "Outputs-Multi")


    
    return geodataframe_metrop, gdf_points               
    




def resetdic(dic):
    """Function which resets a distribution dictionnary by setting all distributions to "unique" (no distribution) """
    dic_reset = copy.deepcopy(dic)
    for a in dic_reset:
        
        if dic_reset[a][0]!='unique' and dic_reset[a][0]!='binary':
            
            dic_reset[a][0]="unique"

    return dic_reset




def geo_simulations_with_geodfasinput_sampling_2(Tech_opdict,  # To modify at the beginning
                    Biodict,  # To modify at the beginning
                    Locationdict,  # To modify at the beginning
                    Physicdict,  # To modify at the beginning
                    Tech_opdict_distributions,  # To modify at the beginning
                    Biodict_distributions,  # To modify at the beginning
                    Locationdict_distributions,  # To modify at the beginning
                    Physicdict_distributions,  # To modify at the beginning
                    LCIdict,  # Do not modify
                    Size_sample,
                    size_montecarlo,
                    cultivationperiod,
                    fraction_max_yield,
                    elemental_contents,  # To modify in the csv
                    fishfeed_table,  # To modify in the csv
                    methods_selected,  # To modify at the beginning
                    list_cfs,
                    categories_contribution,  # To modify at the beginning
                    processes_in_categories,
                    senstype,
                    latsec,
                    biosph,
                    MICAH,
                    Ecoinvent,
                    list_processes_electricity,
                    geodataframe_metrop, 
                    gdf_points):
    """Function which performs the Multi-dimensional sampling, collects and processes the results."""
    
    
    # Create distribution dictionnaries with all distributions set to "unique"
    
    Tech_opdict_distributions_unique = resetdic(Tech_opdict_distributions)
    Locationdict_distributions_unique = resetdic(Locationdict_distributions)
    Physicdict_distributions_unique = resetdic(Physicdict_distributions)
            
    # Generate the strain-compounds pairs 
    
    output_sample_algae = mainfunc.sampling_func(Tech_opdict_distributions_unique,
                  Biodict_distributions,
                  Locationdict_distributions_unique, 
                  Physicdict_distributions_unique,
                  Size_sample, 
                  senstype)
    
    
    sample_algae = output_sample_algae[0]
    
    names_param_algae = output_sample_algae[1]
    
    problem_algae = output_sample_algae[-1]
    
    

    
    time1=time()
    
    
    # Start Ray.
    ray.shutdown()
    
    ray.init()
    
    # Keep track of position in geodataframes
    point_rank_in_points_gdf = -1
    
    number_countries = geodataframe_metrop.shape[0]
    
    number_points = gdf_points.shape[0]
    
    number_algae = len(sample_algae)
    
    #Initialization for collection of results
    # For sensitivity analysis
    
    # Empty array to collect sensitivity indices
    results_for_algae_tot_weighted = np.array([[0]*len(methods_selected)]*number_algae,dtype="float32")
    
    # Not used in the article
    results_for_algae_tot_non_weighted = np.array([[0]*len(methods_selected)]*number_algae,dtype="float32")

    
    table_prod = np.array([[0]*number_points]*number_algae)
    
    
    results_sobol_fast_weighted=[]
    
    # Not used in the article
    results_sobol_fast_non_weighted=[]
    
    for country_index in range(number_countries):    
        # For 1 country

        

        # Modify impacts of 1kWh of electricty to match the national mix
        
        list_imp_mix_nat = geodataframe_metrop['Imp_elec'].iloc[country_index]
        
        dict_mono_technosphere_lcas_nat = modif_impact_elec_national_mix(dict_mono_technosphere_lcas, list_processes_electricity, list_imp_mix_nat)
        
        # Update location dict with lat, long
        
        list_points = geodataframe_metrop['random_points'].iloc[country_index]
        
        name_CNTR = geodataframe_metrop['CNTR_CODE'].iloc[country_index]
    
        index_point = -1
            
        len_list_points = len(list_points)
        
        
        
        # Initialize minimum values for different IC in the country
        list_min_country=[float('inf')]*len(methods_selected)
        list_max_country=[-float('inf')]*len(methods_selected)
       

        
        
        
        # Initialize the list of productions achieved by each location to weigh the mix shares
        list_prod_points = [0]*len_list_points
        
        for point_index in range(len_list_points):
            # For 1 location in the country

            


            # corresponding rank (row) of this point in the geodataframe containing points
            
            point_rank_in_points_gdf += 1
            
            point = list_points[point_index]
            
            
            index_point += 1
    
            # Change the locationdict accordingly
            
            Locationdict = {'lat': round(point[1],3), #°
                    'long': round(point[0],3), #°
                    'Twell': 10,  # °C  # Will be modified during stochastic sampling
                    'depth_well': 10,  # m  # Will be modified during stochastic sampling
                    'azimuthfrontal': 90} # ° # Will be modified during stochastic sampling
            
            

            
            
            # Initialize minimum values for different IC for this point

            list_min_point = [float('inf')]*len(methods_selected)
            list_max_point = [-float('inf')]*len(methods_selected)
            
            

            a=time()     
                
            
            # Perform Multi-dimensional sampling for this location point, 
            # Calculate LCAs and  collect results.  
            # Calculates results for all Strain-compounds on this point
            
            res_point = mainfunc.final_function_simulations_sampling_2_SPECIAL(dict_mono_technosphere_lcas_nat,
                                             Tech_opdict,  # To modify at the beginning
                                             Biodict,  # To modify at the beginning
                                             Locationdict,  # To modify at the beginning
                                             Physicdict,  # To modify at the beginning
                                             Tech_opdict_distributions,  # To modify at the beginning
                                             Biodict_distributions,  # To modify at the beginning
                                             Locationdict_distributions,  # To modify at the beginning
                                             Physicdict_distributions,  # To modify at the beginning
                                             LCIdict,  # Do not modify
                                             size_montecarlo,  # Sample size
                                             cultivationperiod,  # Paper value
                                             fraction_max_yield,  # Paper value
                                             elemental_contents,  # To modify in the csv
                                             fishfeed_table,  # To modify in the csv
                                             methods_selected,  # To modify at the beginning
                                             list_cfs,
                                             categories_contribution,  # To modify at the beginning
                                             processes_in_categories,  # To modify at the beginning
                                             dict_mono_technosphere_lcas,
                                             sample_algae,
                                             names_param_algae)  
              

            
            
            
            b = time()-a  # For indication
                
                

            
            # Reminder: Output of final_function_simulations :
                

            # return (sample_mc, # 0
            #         results_table_df,  #1
            #         list_desc_stat_algae,  #2
            #         list_total_desc_stat_contri_df)  #3
            

               
                
            # Results of all simulations with LCIA results
            
            results_table_df_point =  res_point[1]
            
            
             # Statistical Description of results
            
            list_desc_stat_algae_for_this_point = res_point[2]
            

            # Statistical Description of contributions
                
            list_total_desc_stat_contri_df = res_point[3]           

            # Save memory
            del res_point
            
            # Add columns coord and location
    
            
            results_table_df_point['lat']=point[1]
            
            results_table_df_point['long']=point[0]
            
            results_table_df_point['CNTR']=geodataframe_metrop['CNTR_CODE'].iloc[country_index]
            
            
            #print(type(results_table_df_point["GWP100"].iloc[2]))
            
            # Add columns with impact for 1 kWh of electricity 
            
            for meth_index in range(len(methods_selected)):
                
                meth = methods_selected[ meth_index]
                
                name_col = meth[-1] + "_1kWh"
                
                results_table_df_point[name_col] = list_imp_mix_nat[meth_index]
                
                    
            
           
            suffix_file = name_CNTR+'_'+str(index_point)
            
            name_file_results = 'results_table_df_point' + suffix_file
    
                    
            

            results_table_df_point.to_csv('../Outputs-Multi/'+name_file_results+'.csv', sep=';', encoding='utf-8')
    
        
    
                

                

            
            
    
            # Browse the dataframes containing the statistical description for each Strain-compound for this point
            for index_alga in range(number_algae):   
                
                
                # Statistical descriptions of contributions:
                
                    # Add columns coord and location
        
                desc_contributions_alga_for_this_point = list_total_desc_stat_contri_df[index_alga]
                
                desc_contributions_alga_for_this_point['lat']=point[1]
                
                desc_contributions_alga_for_this_point['long']=point[0]
                
                desc_contributions_alga_for_this_point['CNTR']=geodataframe_metrop['CNTR_CODE'].iloc[country_index]
                
                suffix_file = name_CNTR+'_'+str(index_point)+'_'+str(index_alga)

                
                name_file_desc_contri = 'desc_contributions_point' + suffix_file
                
                desc_contributions_alga_for_this_point.to_excel('../Outputs-Multi/'+name_file_desc_contri+'.xlsx', 
                             sheet_name='Sheet1',
                             na_rep='', 
                             float_format=None,
                             columns=None,
                             header=True, 
                             index=True,
                             index_label=None,
                             startrow=0,
                             startcol=0,
                             engine=None,
                             merge_cells=True,
                             encoding=None,
                             inf_rep='inf',
                             verbose=True,
                             freeze_panes=None,
                             storage_options=None)
            
                
                
                
                # Statistical description of LCI and LCA results
                
                desc_LCI_for_alga_for_this_point = list_desc_stat_algae_for_this_point[index_alga]
                
                prod_mean = desc_LCI_for_alga_for_this_point["Total production harvested kg dw.m-2"].loc['mean']
                    
                # Save the value for the weighing
                
                table_prod[index_alga,point_rank_in_points_gdf]=prod_mean
                
                
                # Update min and max and update the list of averages 
                
                for meth_index in range(len(methods_selected)):
                    
                    
                    name_meth = methods_selected[meth_index][-1]
    
                    
                    # Update min and max for the country
                    
                    if list_min_country[meth_index]>desc_LCI_for_alga_for_this_point[name_meth].loc['min']:
                       
                        list_min_country[meth_index] = desc_LCI_for_alga_for_this_point[name_meth].loc['min']
                    
                    if list_max_country[meth_index]<desc_LCI_for_alga_for_this_point[name_meth].loc['max']:
                        
                        list_max_country[meth_index] = desc_LCI_for_alga_for_this_point[name_meth].loc['max']
        
                    
                    # Update min and max for the point

                    if list_min_point[meth_index]>desc_LCI_for_alga_for_this_point[name_meth].loc['min']:
                       
                        list_min_point[meth_index] = desc_LCI_for_alga_for_this_point[name_meth].loc['min']
                    
                    if list_max_point[meth_index]<desc_LCI_for_alga_for_this_point[name_meth].loc['max']:
                        
                        list_max_point[meth_index] = desc_LCI_for_alga_for_this_point[name_meth].loc['max']
        

            
        
                    

        
        
        
                    #Update list averages for strain-compounds

                   
                    results_for_algae_tot_weighted[index_alga,meth_index] += desc_LCI_for_alga_for_this_point[name_meth].loc['50%']*prod_mean
                    

                    
                    results_for_algae_tot_non_weighted[index_alga,meth_index] += desc_LCI_for_alga_for_this_point[name_meth].loc['50%']
                    
                    
                    # Update geodataframes
                    
                    name_col_std='Std'+'_'+name_meth
                    
                    name_col_mean='Mean'+'_'+name_meth
                    
                    name_col_min='Min'+'_'+name_meth
                    
                    name_col_max='Max'+'_'+name_meth
                    
        
                    # Geodataframe countries 
                    
                    geodataframe_metrop[name_col_mean].iloc[country_index]=(geodataframe_metrop[name_col_mean].iloc[country_index] 
                                                         + desc_LCI_for_alga_for_this_point[name_meth].loc['mean']/(len_list_points*number_algae))
                    
                    #Geodataframe points 
                    
                    gdf_points[name_col_mean].iloc[point_rank_in_points_gdf]=(gdf_points[name_col_mean].iloc[point_rank_in_points_gdf]
                                                                              +desc_LCI_for_alga_for_this_point[name_meth].loc['mean']/number_algae)
                    
                    
                    

            
                    
            # All STrain-compounds were drawn : update columns for the point     
          
            
          
            # Means have already been alculated during the sampling
      
            # Update min and max for the point
            
            for meth_index in range(len(methods_selected)):
                
                
                name_meth = methods_selected[meth_index][-1]

                
                name_col_std='Std'+'_'+name_meth
                
                name_col_mean='Mean'+'_'+name_meth
                
                name_col_min='Min'+'_'+name_meth
                
                name_col_max='Max'+'_'+name_meth
            

                #Geodataframe points

                #Min
                
                gdf_points[name_col_min].iloc[point_rank_in_points_gdf]=gdf_points[name_col_min].iloc[point_rank_in_points_gdf]+list_min_point[meth_index]
                
                
                # Max
        
                #Geodataframe points
        
                gdf_points[name_col_max].iloc[point_rank_in_points_gdf]=gdf_points[name_col_max].iloc[point_rank_in_points_gdf]+list_max_point[meth_index]
        
                
        
        
                
           
                   
        
                
        # End of loop for points in country    
        
        
        # Update columns for country
        
        for meth_index in range(len(methods_selected)):
                
                # Update min and max for the country
                
                
                name_meth = methods_selected[meth_index][-1]
                
                name_col_std='Std'+'_'+name_meth
                
                name_col_min='Min'+'_'+name_meth
                name_col_max='Max'+'_'+name_meth            
                
                # Geodataframe countries


                # + allows to conserve consistency in the object type (float32)

                geodataframe_metrop.iloc[country_index, geodataframe_metrop.columns.get_loc(name_col_min)] = geodataframe_metrop.iloc[country_index, geodataframe_metrop.columns.get_loc(name_col_min)]+ list_min_country[meth_index]
                geodataframe_metrop.iloc[country_index, geodataframe_metrop.columns.get_loc(name_col_max)] = geodataframe_metrop.iloc[country_index, geodataframe_metrop.columns.get_loc(name_col_min)]+list_max_country[meth_index]
  
    
  
    
  
    
    
    ray.shutdown()
    
    # Weigh impacts with productivities or not

    for alga_index in range(number_algae):
        
        # Weigh impact with productivities 
        # Not used in the article
        results_for_algae_tot_weighted[alga_index,:] = results_for_algae_tot_weighted[alga_index,:]/sum(table_prod[alga_index,:])
                

        results_for_algae_tot_non_weighted[alga_index,:] = results_for_algae_tot_non_weighted[alga_index,:]/number_points
                   


    # Sensitivity analysis
    
    for meth_index in range(len(methods_selected)):  # For each column of LCIA results = for each method (Impact category)
        
        name_meth = methods_selected[meth_index][-1]
        
        output_weighted = results_for_algae_tot_weighted[:, meth_index] # Call the result column for this impact category
        
        output_non_weighted = results_for_algae_tot_non_weighted[:, meth_index] # Call the result column for this impact category


        # # Rearranging into a monodimensional array
         
        # array_output = pd.DataFrame(output).to_numpy()
    
        # flat_array = array_output.flat
    
        # output_list_clean = []
    
        #  for a in flat_array:
        #      output_list_clean.append(a)
    
        #  output_clean = np.array(output_list_clean)
    
        # Performing the sensitivity analysis
   
        if senstype == 'SOBOL':
            
            sobol_fast_res_weighted = SALib.analyze.sobol.analyze(
                problem_algae, output_weighted, calc_second_order=False)            
            
            sobol_fast_res_non_weighted = SALib.analyze.sobol.analyze(
                problem_algae, output_non_weighted, calc_second_order=False)
            
        elif senstype == 'FAST':
            
            sobol_fast_res_weighted = SALib.analyze.fast.analyze(
                problem_algae, output_weighted)            
            
            sobol_fast_res_non_weighted = SALib.analyze.fast.analyze(
                problem_algae, output_non_weighted)            
   
    
    
    
        results_sobol_fast_weighted.append([name_meth, sobol_fast_res_weighted])
        
        results_sobol_fast_non_weighted.append([name_meth, sobol_fast_res_non_weighted])

    
    
    
    
    # COnverting dicitonnaries into dataframes of  Sensitivity indices
    
    
    name_columns = [IC[0]+'__'+ type_ind  for IC in results_sobol_fast_weighted for type_ind in list(results_sobol_fast_weighted[0][1])]                

    count=-1
    if senstype == 'SOBOL':
        
        name_columns=[IC[0]+'__'+ type_ind  for IC in results_sobol_fast_weighted for type_ind in list(results_sobol_fast_weighted[0][1])]                

        sensi_multi_weighted = pd.DataFrame(np.zeros((len(results_sobol_fast_weighted[0][1]['ST']),
                                             len(results_sobol_fast_weighted)*4)),
                     columns=name_columns)

        sensi_multi_non_weighted = pd.DataFrame(np.zeros((len(results_sobol_fast_non_weighted[0][1]['ST']),
                                             len(results_sobol_fast_non_weighted)*4)),
                     columns=name_columns)        
        # put values in datframes (1 column = 1 Impact category, n rows=n paramters)

        
        for IC_index in range(len(results_sobol_fast_weighted)) :
            
            count += 1
            
            IC_weighted = results_sobol_fast_weighted[IC_index]
            IC_non_weighted = results_sobol_fast_non_weighted[IC_index]
        
            sensi_multi_weighted.iloc[:,count*4]=IC_weighted[1]['S1']
            sensi_multi_weighted.iloc[:,count*4+1]=IC_weighted[1]['S1_conf']
            sensi_multi_weighted.iloc[:,count*4+2]=IC_weighted[1]['ST']
            sensi_multi_weighted.iloc[:,count*4+3]=IC_weighted[1]['ST_conf']   
            
            sensi_multi_non_weighted.iloc[:,count*4]=IC_non_weighted[1]['S1']
            sensi_multi_non_weighted.iloc[:,count*4+1]=IC_non_weighted[1]['S1_conf']
            sensi_multi_non_weighted.iloc[:,count*4+2]=IC_non_weighted[1]['ST']
            sensi_multi_non_weighted.iloc[:,count*4+3]=IC_non_weighted[1]['ST_conf'] 
        
    elif senstype == 'FAST': # No confidence interval with this Salib version
        
        name_columns=[IC[0]+'__'+ type_ind  for IC in results_sobol_fast_weighted for type_ind in list(results_sobol_fast_weighted[0][1])[:-1]]                

        

        
        
        sensi_multi_weighted = pd.DataFrame(np.zeros((len(results_sobol_fast_weighted[0][1]['ST']),
                                             len(results_sobol_fast_weighted)*2)),
                     columns=name_columns)

        sensi_multi_non_weighted = pd.DataFrame(np.zeros((len(results_sobol_fast_non_weighted[0][1]['ST']),
                                             len(results_sobol_fast_non_weighted)*2)),
                     columns=name_columns) 
        
        # put values in datframes (1 column = 1 Impact category, n rows=n paramters)

        
        for IC_index in range(len(results_sobol_fast_weighted)) :
            
            count += 1
            
            IC_weighted = results_sobol_fast_weighted[IC_index]
            IC_non_weighted = results_sobol_fast_non_weighted[IC_index]
        
     
            sensi_multi_weighted.iloc[:,count*2]=IC_weighted[1]['S1']
            sensi_multi_weighted.iloc[:,count*2+1]=IC_weighted[1]['ST']
            
            sensi_multi_non_weighted.iloc[:,count*2]=IC_non_weighted[1]['S1']
            sensi_multi_non_weighted.iloc[:,count*2+1]=IC_non_weighted[1]['ST']    
    
  
    
    # put values in datframes (1 column = 1 Impact category, n rows=n paramters)
    
    sensi_multi_weighted['Parameter'] = names_param_algae
        
    sensi_multi_non_weighted['Parameter'] = names_param_algae

    
    sensi_multi_weighted_melt=pd.melt(sensi_multi_weighted,id_vars=['Parameter'],value_vars=name_columns,var_name="Type")
    
    sensi_multi_non_weighted_melt=pd.melt(sensi_multi_non_weighted,id_vars=['Parameter'],value_vars=name_columns,var_name="Type")

     
    
    
    
    sensi_multi_non_weighted_melt[["Impact_Category","Type indice"]] = sensi_multi_non_weighted_melt["Type"].str.split("__",1,expand=True)
    
    sensi_multi_weighted_melt[["Impact_Category","Type indice"]] = sensi_multi_weighted_melt["Type"].str.split("__",1,expand=True)
    
    
    # Export sensitivity indices
    
    x = datetime.datetime.now()
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))

    name_file_sensi_multi_weighted = 'sensi_multi_sampling_2_weighted' + "_"+month+"_"+day+"_"+microsec
    name_file_sensi_multi_non_weighted = 'sensi_multi_sampling_2_non_weighted' + "_"+month+"_"+day+"_"+microsec
    
    sensi_multi_weighted_melt.to_excel('../Outputs-Multi/'+name_file_sensi_multi_weighted+'.xlsx', 
                  sheet_name='Sheet1',
                  na_rep='', 
                  float_format=None,
                  columns=None,
                  header=True, 
                  index=True,
                  index_label=None,
                  startrow=0,
                  startcol=0,
                  engine=None,
                  merge_cells=True,
                  encoding=None,
                  inf_rep='inf',
                  verbose=True,
                  freeze_panes=None,
                  storage_options=None)
    
    
    sensi_multi_non_weighted_melt.to_excel('../Outputs-Multi/'+name_file_sensi_multi_non_weighted+'.xlsx', 
                  sheet_name='Sheet1',
                  na_rep='', 
                  float_format=None,
                  columns=None,
                  header=True, 
                  index=True,
                  index_label=None,
                  startrow=0,
                  startcol=0,
                  engine=None,
                  merge_cells=True,
                  encoding=None,
                  inf_rep='inf',
                  verbose=True,
                  freeze_panes=None,
                  storage_options=None)    
    
    
    
    # Export geodataframes 
    
    month=str(x.month)
    day=str(x.day)
    microsec=str(x.strftime("%f"))
                 
    name_geo_countries ='geodataframe_output'+"_"+month+"_"+day+"_"+microsec
    
    name_geo_points ='gdfpoints_outputs'+"_"+month+"_"+day+"_"+microsec


    export_pickle_2(geodataframe_metrop, name_geo_countries, "Outputs-Multi")
    
    export_pickle_2(gdf_points, name_geo_points, "Outputs-Multi")
    
    
    

    
    # Export sample Strain-compounds
    
    dataframe_sample = pd.DataFrame(sample_algae,columns=names_param_algae)
    
    dataframe_sample['alga'] = dataframe_sample.index
    
    
    name_sample_algae ="algae_sample" +"_"+month+"_"+day+"_"+microsec
    
    dataframe_sample.to_excel('../Outputs-Multi/'+name_sample_algae+'.xlsx', 
                             sheet_name='Sheet1',
                             na_rep='', 
                             float_format=None,
                             columns=None,
                             header=True, 
                             index=True,
                             index_label=None,
                             startrow=0,
                             startcol=0,
                             engine=None,
                             merge_cells=True,
                             encoding=None,
                             inf_rep='inf',
                             verbose=True,
                             freeze_panes=None,
                             storage_options=None) 
    
    # Export Simulation information 
    
    info_simu = pd.DataFrame(np.array([[dataframe_sample.shape[0]]]))
    
    name_info_simu = "info_simu" +"_"+month+"_"+day+"_"+microsec

    info_simu.to_excel('../Outputs-Multi/'+name_info_simu+'.xlsx', 
                         sheet_name='Sheet1',
                         na_rep='', 
                         float_format=None,
                         columns=None,
                         header=True, 
                         index=True,
                         index_label=None,
                         startrow=0,
                         startcol=0,
                         engine=None,
                         merge_cells=True,
                         encoding=None,
                         inf_rep='inf',
                         verbose=True,
                         freeze_panes=None,
                         storage_options=None) 
    
    time2=time()
    
    timetot=time2-time1
    
    print("timetot",timetot)


    return geodataframe_metrop, gdf_points,sample_algae
    



    
    
    
"""CALL THE FUNCTIONS"""








#####
""" This generates the same grid of locations as the one use di nthe article.
 Indeed, the random seed has been fixed to 1 in the Script Map-2nd  in line 186."""
#####
geodataframe_metrop,gdf_points = pre_download_climatic_data(latsec,methods_selected,biosph,MICAH,Ecoinvent,cultivationperiod )



""" Launch the simulations"""






    
# """ Or  Download a grid of locations."""

# def importpickle(path):
#     "Function which imports pickle objects"
#     with open(path, 'rb') as pickle_load:
#         obj = pickle.load(pickle_load)
#     return obj    
        

# """Import the locations grid used in the article"""


# geodataframe_metrop = importpickle("..\Locations grid\gdfpoints_input_original.pkl")

# gdf_points = importpickle("..\Locations grid\geodataframe_input_original.pkl")

    
    
    
  


print('Now calculating')







geodataframe_metrop_2, gdf_points__2,sample_algae_2 = geo_simulations_with_geodfasinput_sampling_2(Tech_opdict,  # To modify at the beginning
                    Biodict,  # To modify at the beginning
                    Locationdict,  # To modify at the beginning
                    Physicdict,  # To modify at the beginning
                    Tech_opdict_distributions,  # To modify at the beginning
                    Biodict_distributions,  # To modify at the beginning
                    Locationdict_distributions,  # To modify at the beginning
                    Physicdict_distributions,  # To modify at the beginning
                    LCIdict,  # Do not modify
                    Size_sample,
                    size_montecarlo,
                    cultivationperiod,
                    fraction_max_yield,
                    elemental_contents,  # To modify in the csv
                    fishfeed_table,  # To modify in the csv
                    methods_selected,  # To modify at the beginning
                    list_cfs,
                    categories_contribution,  # To modify at the beginning
                    processes_in_categories,
                    senstype,
                    latsec,
                    biosph,
                    MICAH,
                    Ecoinvent,
                    list_processes_electricity,
                    geodataframe_metrop, 
                    gdf_points)

