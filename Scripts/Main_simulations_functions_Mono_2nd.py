# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:45:34 2021

@author: Pierre Jouannais, Department of Planning, DCEA, Aalborg University
pijo@plan.aau.dk

"""
'''
Script containing the main functions calling the other scripts to calculate the LCAs in Multi-dimensional sampling.


'''

import datetime
from time import *
import requests
import pickle
import cProfile
from scipy.integrate import odeint
import numpy as np
import os
import pandas as pd
import decimal
from random import *
import pstats
from itertools import *
from math import*
import csv
import copy

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

import Cultivation_simul_2nd as cultsimul
import Functions_for_physical_and_biological_calculations_2nd as functions





import ray
import time



# Set working directory to file location 
# (works only when executing the whole file and not only sections (Run Current cell))

currentfolder=os.path.dirname(os.path.realpath(__file__))
os.chdir(currentfolder)










def waste_water_impact_biosphere(my_conc_N,
                      my_conc_P,
                      my_conc_C,
                      my_conc_Mg,
                      my_conc_K,
                      my_conc_S,
                      methods,
                      list_cfs):
    

    '''Function which :
        -Establishes the mass balances which match the input
        concentrations for a given microalgal waste water.
        
        -Calculates the impact of these emissions for each impact category ( method)

        Inputs :
            
            #my_conc_N, my_conc_P, my_conc_C, my_conc_Mg, my_conc_K, my_conc_S:
                concentrations in g.L-1 of different elements in the actual
                microalgal wastewater entering the treatment
                
            #methods: The list of impact assessment methods.    

        Outputs :
            
            #list_sum_impacts_biosphere_waste_water : list containing the total 
            impacts due to biosphere emissions for the treatment of 1 cubic meter of wastewater.
            1 element per Impact category

        '''
    
    # Defining materials from original wastewater activity
    
    
    # All original values of the output flows classified by element
    # ['ID of the exchange, amount of the exchange']
    
    N_original_flows = [['ae70ca6c-807a-482b-9ddc-e449b4893fe3', 0.00049],
                        ['0017271e-7df5-40bc-833a-36110c1fe5d5', 0.000644],
                        ['6dc1b46f-ee89-4495-95c4-b8a637bcd6cb', 0.0001146],
                        ['d068f3e2-b033-417b-a359-ca4f25da9731', 0.00067568],
                        ['b61057a3-a0bc-4158-882e-b819c4797419', 2.393e-05],
                        ['13331e67-6006-48c4-bdb4-340c12010036', 0.011027],
                        ['9990b51b-7023-4700-bca0-1a32ef921f74', 0.00059438],
                        ['7ce56135-2ca5-4fba-ad52-d62a34bfeb35', 0.048295]]
    
    P_original_flows = [['490b267b-f429-4d9a-ac79-224e37fb4d58', 7.2652e-05],
                        ['1727b41d-377e-43cd-bc01-9eaba946eccb', 0.0027476],
                        ['329fc7d8-4011-4327-84e4-34ff76f0e42d', 2.705e-05],
                        ['198ce8e3-f05a-4bec-9f7f-325347453326', 6.2034e-07]]
    
    C_original_flows = [['f65558fb-61a1-4e48-b4f2-60d62f14b085', 0.0072992],
                        ['8734eb08-50cf-4f5a-8d1a-db76d38efe3c', 4.8293e-05],
                        ['725c7923-0ed8-43e5-b485-fad7e34bef08', 4.8293e-05],
                        ['62859da4-f3c5-417b-a575-8b00d8d658b1', 0.012346],
                        ['73ed05cc-9727-4abf-9516-4b5c0fe54a16', 0.17202],
                        ['960c0f37-f34c-4fc1-b77c-22d8b35fd8d5', 0.0075377],
                        ['9afa0173-ecbd-4f2c-9c5c-b3128a032812', 0.0001613],
                        ['baf58fc9-573c-419c-8c16-831ac03203b9', 0.00050213]]
    
    S_original_flows = [['bfc0bf1c-e5e2-4702-a502-08c892031837', 0.0011039],
                        ['d4049741-cef2-4edd-a3af-9728b9e3a568', 0.0010988],
                        ['8c52f40c-69b7-4538-8923-b371523c71f5', 0.000884],
                        ['37d35fd0-7f07-4b9b-92eb-de3c27050172', 0.14465]]
    
    Mg_original_flows = [['ce9fd912-233a-4807-a33e-0323b1e4a7a2', 0.00014782],
                         ['ebfe261d-ab0d-4ade-8743-183c8c6bdcc6', 2.205e-07],
                         ['e8475907-2081-4fd5-9526-bfcef88380db', 0.00039974],
                         ['7bdab722-11d0-4c42-a099-6f9ed510a44a', 0.0051478]]
    
    K_original_flows = [['1653bf60-f682-4088-b02d-6dc44eae2786', 0.0003989]]
    
    Al_original_flows = [['2baa4381-b781-4f5e-90de-508b0fa3fd1f', 0.0010518],
                         ['97e498ec-f323-4ec6-bcc0-d8a4c853bae3', 6.228e-05],
                         ['01056d4b-f9b0-4dfc-b8d9-8407c8376efb', 0.00031181],
                         ['6f0b8b7c-3888-4174-b7e3-916d42d678ee', 6.5822e-07]]
    
    Na_original_flows = [['1fc409bc-b8e7-48b2-92d5-2ced4aa7bae2', 0.002186]]
    
    Ca_original_flows = [['ac066c02-b403-407b-a1f0-b29ad0f8188f', 0.045852],
                         ['ae28c923-a4a3-4f00-b862-1ae6e748efb9', 0.0012412],
                         ['a912f450-5233-489b-a2e9-8c029fab480f', 2.3777e-06],
                         ['f16fa1da-e426-4820-bf9d-71595c22283b', 0.0035605]]
    
    Fe_original_flows = [['9b6d6f07-ebc6-447d-a3c0-f2017d77d852', 0.0017779],
                         ['7c335b9c-a403-47a8-bb6d-2e7d3c3a230e', 0.0036009],
                         ['db364689-e1a3-4629-8835-e6c59d6daf09', 0.009475],
                         ['32cd0492-c0cb-4898-a2b1-675eedc5b688', 1.2671e-07]]
    Cl_original_flows = [['5e050fab-1837-4c42-b597-ed2f376f768f', 0.040484]]
    
    
    
    # The inputs of different elements in the orginal wastewater treament activity.
    # Added to the waste water as treatment.
    
    added_treatment_N = 2.61388*10**-5
    added_treatment_P = 0
    added_treatment_C = 0
    added_treatment_S = 0.003339284
    added_treatment_Al = 0.000497558
    added_treatment_Fe = 0.0098005
    added_treatment_Na = 9.44323*10**-5
    added_treatment_Ca = 4.17657*10**-7
    added_treatment_Cl = 0.010468958
    added_treatment_Mg = 0
    added_treatment_K = 0
    
    # The total outputs (and thus inputs) of elements in the original activity.
    # Including the added treatments.
    
    totalNoutput = 0.020770454172634983
    totalPoutput = 0.0009270353321544698
    totalCoutput = 0.0746397575218131
    totalSoutput = 0.05012543333333333
    totalAloutput = 0.0014265482200000001
    totalFeoutput = 0.01485392671
    totalNaoutput = 0.002186
    totalCaoutput = 0.050656077699999996
    totalCloutput = 0.040484
    totalMgoutput = 0.0056955805
    totalKoutput = 0.0003989
    
    
    # The actual inputs of elements contained in the waste
    # water of the original activity.
    
    # If the value is negative it means that the mass balance in the original
    # activity was not respected and we assume the element was not in the
    # incoming waste water.
    
    absolute_input_N = max(totalNoutput-added_treatment_N, 0)
    absolute_input_C = max(totalCoutput-added_treatment_C, 0)
    absolute_input_P = max(totalPoutput-added_treatment_P, 0)
    absolute_input_S = max(totalSoutput-added_treatment_S, 0)
    absolute_input_Al = max(totalAloutput-added_treatment_Al, 0)
    absolute_input_Fe = max(totalFeoutput-added_treatment_Fe, 0)
    absolute_input_Na = max(totalNaoutput-added_treatment_Na, 0)
    absolute_input_Ca = max(totalCaoutput-added_treatment_Ca, 0)
    absolute_input_Cl = max(totalCloutput-added_treatment_Cl, 0)
    absolute_input_K = max(totalKoutput-added_treatment_K, 0)
    absolute_input_Mg = max(totalMgoutput-added_treatment_Mg, 0)
    
    
    total_flows=(N_original_flows
                 + P_original_flows
                 + C_original_flows
                 + S_original_flows
                 + Mg_original_flows
                 + K_original_flows
                 + Al_original_flows
                 + Na_original_flows
                 + Ca_original_flows
                 + Fe_original_flows
                 + Cl_original_flows)
    
    

    
    # Initialize the dicitonnary that will contain the impacts associated to each substance in the wastewater
    
    dictionnary_original_flows= {flow[0] : [0]*len(methods) for flow in total_flows} 
    
    # Collect the characterization factors for each impact category
    
    #list_cfs= [bw.Method((meth)).load() for meth in methods]



    meth_index = -1
    
    for meth in methods:
        
        meth_index += 1
        
        
        cfs_dictionnary = { subst[0][1] : subst[1] for subst in list_cfs[meth_index]}
        
        # For all the biosphere flows in the original waste water activity 
        for flow in total_flows: 
            
            if flow in N_original_flows: # If this flow contains nitrogen
                
                # We assume that  the added treatment is
                # equally shared among the flows of a same element.
                
                original_added_treatment = (flow[1] * added_treatment_N/totalNoutput) 
                
                if flow[0] in cfs_dictionnary: # if there is a cf for this flow in this method
        
                    # The impact is cf * new value of the flow in the microalgal wastewater
                    
                    # The new value of the flow is : 
                        
                    # (Original flow - part that comes from the treatment)
                    # * (new concentration waste water/original concentration waste water)
                    # + Share of the treatment ending up in this flow.
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]] * ((flow[1] - original_added_treatment)
                                                 * my_conc_N/absolute_input_N
                                                 + original_added_treatment)
                    
                    # Update the total impact associated to this flow for the right method

                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter
                  
            elif flow in P_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_P/totalPoutput) 
                
                if flow[0] in cfs_dictionnary:

                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_P/absolute_input_P
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter
           
            elif flow in C_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_C/totalCoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_C/absolute_input_C
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter


    
            elif flow in S_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_S/totalSoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_S/absolute_input_S
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter


            elif flow in Mg_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Mg/totalMgoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Mg/absolute_input_Mg
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in K_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_K/totalKoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_K/absolute_input_K
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in Al_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Al/totalAloutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Al/absolute_input_Al
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

                
            elif flow in Na_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Na/totalNaoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Na/absolute_input_Na
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

            elif flow in Ca_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Ca/totalCaoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Ca/absolute_input_Ca
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

            elif flow in Fe_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Fe/totalFeoutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Fe/absolute_input_Fe
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

    
            elif flow in Cl_original_flows:
                
                original_added_treatment = (flow[1] * added_treatment_Cl/totalCloutput) 
                
                if flow[0] in cfs_dictionnary:
                    
                    impact_percubic_meter = cfs_dictionnary[flow[0]]*((flow[1]-original_added_treatment)
                                                 * my_conc_Cl/absolute_input_Cl
                                                 + original_added_treatment)
                    
                    dictionnary_original_flows[flow[0]][meth_index] = impact_percubic_meter

        
        # Summing the impacts of all flows for 1 cubic meter
        
        list_sum_impacts_biosphere_waste_water =[]
        
        for meth_index in range(len(methods)):
            
            sum_impact = sum([dictionnary_original_flows[flow][meth_index] for flow in dictionnary_original_flows ])


            list_sum_impacts_biosphere_waste_water.append(sum_impact)
            
            
            
            
            
        #wastewater_copy.save()
    
        return list_sum_impacts_biosphere_waste_water


     






'''Main calculating functions '''

def LCI_one_strain_uniquevalues(Biodict,
                                Physicdict,
                                Tech_opdict,
                                Locationdict,
                                LCIdict,
                                months_suitable_for_cultivation,
                                fraction_maxyield,
                                fishfeed_table,
                                elemental_contents):
    
    '''Calculate the LCI for one set of parameters given in input by simulating 
    the cultivation and scaling the values to the FU.
    
    Inputs:
        #Biodict : Dictionnary with biological parameters
        #Physicdict : Dictionnary with physical parameters
        #Tech_opdict : Dictionnary with techno-operational parameters
        #LCIdict : Initialized LCI dictionnary
        #months_suitable_for_cultivation : Months for cultivation ; 
        list of month numbers : [a,b,c]
        
        #fraction_maxyield : Fraction of the maximum yield achieved ; .
        #fishfeed_table : DataFrame with fish feed composition
        #elemental_contents : Table with elemental compositons of macronutrients
        

    Outputs:
        # LCIdict_updated : Dictionnary containing the calculated LCI
        
        Other variables resulting from the simulation  for further investigation 
        and review of the code (Non necessary):
            
        # surfaceyield : Areal yield ; kg.m-2 
        # volumetricyield : Volmetric yield ; kg.m-3 
        # optimization_performance : Fish feed Substitution alrogithm performance
        # needed_dbio_check : Non necessary
        # substitution_check : Non necessary
        # total_production_kg_dw : Total production ; kg dw
        # total_production_harvested_kg_dw : Actual harvested production harvested ; kg dw
        # total_production_loss_kg_dw : Production not harvested  ; kg dw
        # conc_waste_water_nutrient_N, : N Concentration in wastewater ; kg.m-3
        # Same for all elements
        #conc_waste_water_biomass : Biomass Concentration  in wastewater ; kg.m-3
        #bioact_molec_dbio : Molecule Concentration  in the dried biomass ; kg. kg dbio -1
        # min_centrifugation_rate_m3_h : Obsolete
        # max_centrifugation_rate_m3_h : Obsolete
        # totalwatercentrifuged : Total volume of water centrifuged ; m3
        # tubelength : Total tube length over 1m2 ; m
        # facilityvolume : Cultivation volume over 1m2 ; m3
        # exchange_area : Exchange area tube/air over 1m2 ; m
        # totalcooling_thermal : Thermal Cooling needed per m2
        (not via Heat Exchanger) ; kWh
        
    '''

    
        
    
    LCIdict_updated = LCIdict.copy()
    
    # Collecting all values from the dictionnaries and creating local variables

    # Tech_opdict 
    
    height = Tech_opdict['height']
    tubediameter = Tech_opdict['tubediameter']
    gapbetweentubes = Tech_opdict['gapbetweentubes']
    horizontaldistance = Tech_opdict['horizontaldistance']
    length_of_PBRunit = Tech_opdict['length_of_PBRunit']
    width_of_PBR_unit = Tech_opdict['width_of_PBR_unit']
    biomassconcentration = Tech_opdict['biomassconcentration']
    flowrate = Tech_opdict['flowrate']
    centrifugation_efficiency = Tech_opdict['centrifugation_efficiency']
    pumpefficiency = Tech_opdict['pumpefficiency']
    slurry_concentration = Tech_opdict['slurry_concentration']
    water_after_drying = Tech_opdict['water_after_drying']
    recyclingrateaftercentrifuge = Tech_opdict['recyclingrateaftercentrifuge']
    roughness = Tech_opdict['roughness']
    rhosuspension = Tech_opdict['rhosuspension']
    cleaningvolumeVSfacilityvolume = Tech_opdict['cleaningvolumeVSfacilityvolume']
    concentration_hypo = Tech_opdict['concentration_hypo']
    concentration_hydro = Tech_opdict['concentration_hydro']
    boilerefficiency = Tech_opdict['boilerefficiency']
    glass_life_expectancy = Tech_opdict['glass_life_expectancy']
    prob_night_monitoring = Tech_opdict['prob_night_monitoring']
    extraction = Tech_opdict['extraction']
    random_market_subst = Tech_opdict['random_market_subst']
    heat_pump = Tech_opdict['heat_pump']
    COP = Tech_opdict['COP']
    
    fertiliser_substitutability_AD = Tech_opdict['fertiliser_substitutability_AD']
    
    carbon_degradibilidy_AD = Tech_opdict['carbon_degradibilidy_AD']
    
    

    heat_input_per_kg_biomass_AD= Tech_opdict['heat_input_per_kg_biomass_AD']
    
    elec_input_per_kg_biomass_AD= Tech_opdict['elec_input_per_kg_biomass_AD']
    
    
    

    # print("prob_night_monitoring",prob_night_monitoring)

    # print("tubediameter",tubediameter)

    # Physicdict

    Cp = Physicdict['Cp']
    hconv = Physicdict['hconv']
    rhomedium = Physicdict['rhomedium']  
    rhowater = Physicdict['rhowater']  
    Cw = Physicdict['Cw']  
    CH4_LHV = Physicdict['CH4_LHV']  

    # Biodict

    lipid_af_dw = Biodict['lipid_af_dw']
    rhoalgae = Biodict['rhoalgae']
    MJ_kglip = Biodict['MJ_kglip']
    MJ_kgcarb = Biodict['MJ_kgcarb']
    MJ_kgprot = Biodict['MJ_kgprot']
    PAR = Biodict['PAR']
    losspigmentantenna = Biodict['losspigmentantenna']
    quantumyield = Biodict['quantumyield']
    lossvoltagejump = Biodict['lossvoltagejump']
    losstoATPNADPH = Biodict['losstoATPNADPH']
    losstohexose = Biodict['losstohexose']
    lossrespiration = Biodict['lossrespiration']
    bioact_fraction_molec = Biodict['bioact_fraction_molec']
    prob_no3 = Biodict['prob_no3']
    Topt = Biodict['Topt']
    T_plateau = Biodict['T_plateau']
    
    # print("T_plateau",T_plateau)

    # Conversion to Tmax, Tmin for simpler calculation
    Tmax = Topt+T_plateau/2
    Tmin = Topt-T_plateau/2

    dcell = Biodict['dcell']
    incorporation_rate = Biodict['incorporation_rate']
    ash_dw = Biodict['ash_dw']
    nutrient_utilisation = Biodict['nutrient_utilisation']
    co2_utilisation = Biodict['co2_utilisation']
    phospholipid_fraction = Biodict['phospholipid_fraction']

    # Locationdict

    lat = Locationdict['lat']
    long = Locationdict['long']
    depth_well = Locationdict['depth_well']
    azimuthfrontal = Locationdict['azimuthfrontal']
    
    # print("lat",lat)

    # We assume that the plant will be purposely located next to a water source with Twell<TMin  :
    
    Twell=np.random.uniform(5, Tmin,size=None)   
    
    
        

    # Qualitative parameters are determined based on the probabilities
    # and a new entry is created in the LCI dict

    # Nsource

    if random() < prob_no3: # Then the source is Nitrate
        Nsource = 'no3'
    else:
        Nsource = 'nh3'

    LCIdict_updated['Nsource'] = Nsource


   

    dice =random()
    if dice<(1/3):
        #print('okC1')
        biochemicalclass='lip'
        
    elif (1/3)<dice<(2/3)  :
        #print('okC2')
        biochemicalclass='prot'
        
    else:
        #print('okC2')
        biochemicalclass='carb'


    # The molecule is a carbohydrate

    #biochemicalclass = 'carb'

    LCIdict_updated['Bio_class'] = biochemicalclass

    #print(random_market_subst)
    # Thermoregulation at night
    if random() < prob_night_monitoring:
        
        night_monitoring = 'yes'

    else:
        night_monitoring = 'no'

    LCIdict_updated['night_monitoring'] = night_monitoring

    # Market for substitution
   # random_market_subst=0.8

    if random_market_subst < (1/3):
        
        market_for_substitution = 'animal feed'

    elif (1/3)<random_market_subst< (2/3):
        
        market_for_substitution = 'fish feed'  
        
    elif random_market_subst > (2/3):
        
        market_for_substitution = 'anaerobic digestion' 

    #print(market_for_substitution)

    LCIdict_updated['market_for_substitution'] = market_for_substitution


    # Collecting PBR geometry

    geom = functions.PBR_geometry(height,
                                  tubediameter,
                                  gapbetweentubes,
                                  horizontaldistance,
                                  length_of_PBRunit,
                                  width_of_PBR_unit)
    tubelength = geom[1]
    facilityvolume = geom[0]
    exchange_area = geom[-1]

    # LCI values which do not depend on the cultivation simulation

    # Calculating biomass composition at different levels

    biomass_composition = functions.biomasscompo(
        lipid_af_dw,
        ash_dw,
        water_after_drying,
        phospholipid_fraction,
        elemental_contents)


    # ash-free biomass composition (af_dw)
    prot_af_dw = biomass_composition[0]
    carb_af_dw = biomass_composition[1]

    # Including ash (dw)
    lip_dw = biomass_composition[2]
    prot_dw = biomass_composition[3]
    carb_dw = biomass_composition[4]

    # After harvesting and drying  (dbio)
    lip_dbio = biomass_composition[5]
    prot_dbio = biomass_composition[6]
    carb_dbio = biomass_composition[7]
    ash_dbio = biomass_composition[8]

    # Elementary composition

    C_af_dw = biomass_composition[9]
    N_af_dw = biomass_composition[10]
    P_af_dw = biomass_composition[11]
    K_af_dw = biomass_composition[12]
    Mg_af_dw = biomass_composition[13]
    S_af_dw = biomass_composition[14]

    C_dw = biomass_composition[15]
    N_dw = biomass_composition[16]
    P_dw = biomass_composition[17]
    K_dw = biomass_composition[18]
    Mg_dw = biomass_composition[19]
    S_dw = biomass_composition[20]

    # Calculating the absolute bioactive molecule content in the dried biomass
    if biochemicalclass == 'lip':
        bioact_molec_dbio = bioact_fraction_molec * lip_dbio

    elif biochemicalclass == 'carb':
        bioact_molec_dbio = bioact_fraction_molec * carb_dbio

    elif biochemicalclass == 'prot':
        bioact_molec_dbio = bioact_fraction_molec * prot_dbio

    # Nutrients

    # Nitrogen consumption 

    # considering ash content

    N_demand = N_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    # Recycling part of the nutrients with supernatant
    N_input = (N_demand / nutrient_utilisation + N_demand *
               recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)

    N_waste = N_input - N_demand

    #Only the correct source of N is updated

    if Nsource == 'nh3':
        # Already as N in Ecoinvent
        LCIdict_updated['market for ammonium sulfate, as N PBR'] = N_input
        
        LCIdict_updated['market for calcium nitrate PBR'] = 0

    if Nsource == 'no3':
        LCIdict_updated['market for ammonium sulfate, as N PBR'] = 0
        
        # Conversion from N to Calcium nitrate 
        LCIdict_updated['market for calcium nitrate PBR'] = N_input/0.15


    # Phosphorus consumption 

    P_demand = P_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    P2O5_demand = P_demand/0.4366     # Conversion P to P2O5

    P2O5_input = (P2O5_demand / nutrient_utilisation + P2O5_demand *
                  recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recylcing

    P_waste = P2O5_input*0.4366 - P_demand

    LCIdict_updated['P source production PBR'] = P2O5_input

    # C
    
    C_demand = C_dw * ((1 / bioact_molec_dbio) * (1 / (1 - water_after_drying)))
    
    CO2_demand = C_demand * (44/12)  # Conversion C to CO2

    CO2_input = CO2_demand / co2_utilisation
    
    CO2_direct_emission = CO2_input - CO2_demand

    LCIdict_updated['Microalgae CO2 PBR'] = CO2_input
    LCIdict_updated['CO2 direct emissions PBR'] = CO2_direct_emission

    # K

    K_demand = K_dw * ((1/bioact_molec_dbio) * (1/(1 - water_after_drying)))

    K2O5_demand = K_demand*1.2 # Conversion to K2O5

    K2O5_input = (K2O5_demand / nutrient_utilisation + K2O5_demand *
                  recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recycling

    K2O5_waste = K2O5_input - K2O5_demand
    K_waste = K2O5_input/1.2 - K_demand

    LCIdict_updated['K source production PBR'] = K2O5_input

    # Mg
    
    Mg_demand = Mg_dw * ((1 / bioact_molec_dbio)*(1/(1 - water_after_drying)))
    MgSO4_demand = Mg_demand * (120.4/24.3)

    MgSO4_input = (MgSO4_demand / nutrient_utilisation + MgSO4_demand *
                   recyclingrateaftercentrifuge) / (1 + recyclingrateaftercentrifuge)  # Recycling

    Mg_input = MgSO4_input * (24.3/120.4)

    Mg_waste = Mg_input-Mg_demand

    LCIdict_updated['Mg source production PBR'] = MgSO4_input

    # S
    
    S_demand = S_dw * ((1/bioact_molec_dbio) * (1/(1-water_after_drying)))

    S_input = MgSO4_input*(32/120.4)
    
    if Nsource == 'nh3':  # Then ammonium sulfate also brings sulfate
    
        # N input --> (NH4)2SO4 input --> S input
        
        S_input += (N_input/0.21) * 0.24
        
    S_waste = S_input-S_demand




    # Cultivation simulation

    # Intializing variables
    totalcooling_thermal = 0
    totalcooling = 0
    totalheating = 0
    totalproduction = 0
    totalproduction_loss = 0
    totalproduction_harvested = 0
    totalwaterpumpedfromthefacility = 0
    totalwaterpumpedfromthewell = 0
    total_elec_centrifuge = 0
    total_elec_mixing = 0
    totalwatercentrifuged = 0

    meantemp_at_harvest_time_cultivationperiod = 0
    min_centrifugation_rate_list = []
    max_centrifugation_rate_list = []

    # Simulating an average day for each month of the cultivation period
    for month in months_suitable_for_cultivation:

        # Calling the cultivation simulation function
        simulation_averageday = cultsimul.cultivation_simulation_timestep10(hconv,
                                                                            Twell,
                                                                            depth_well,
                                                                            lat,
                                                                            long,
                                                                            azimuthfrontal,  
                                                                            month,  
                                                                            Cp,  
                                                                            height,
                                                                            tubediameter,
                                                                            gapbetweentubes,
                                                                            horizontaldistance,
                                                                            length_of_PBRunit,
                                                                            width_of_PBR_unit,  
                                                                            rhoalgae, 
                                                                            rhomedium,
                                                                            rhosuspension, 
                                                                            dcell, 
                                                                            Tmax,
                                                                            Tmin,
                                                                            Biodict, 
                                                                            ash_dw,
                                                                            Nsource,  
                                                                            fraction_maxyield,  
                                                                            biomassconcentration,
                                                                            flowrate,  
                                                                            centrifugation_efficiency,
                                                                            pumpefficiency,  
                                                                            slurry_concentration,
                                                                            water_after_drying,
                                                                            recyclingrateaftercentrifuge,
                                                                            night_monitoring,
                                                                            elemental_contents)

        # Collecting results and multiplying by 
        # average number of days in a month : 30.4

        monthly_heating_energy = simulation_averageday[1]*30.4 #kWh
      
        monthly_waterpumped_from_the_facility = simulation_averageday[2]*30.4  # L

        monthly_waterpumped_from_the_well = simulation_averageday[3]*30.4  # L

        monthly_production = simulation_averageday[4]*30.4  # g dw

        monthly_production_harvested = simulation_averageday[5]*30.4  # g dw

        monthly_production_loss = simulation_averageday[6]*30.4  # g dw

        monthly_volumetric_yield = simulation_averageday[7]*30.4  # g dw.m-3

        monthly_energy_tocentrifuge = simulation_averageday[8]*30.4  # kWh
        
        collectedtemperatureevolution = simulation_averageday[9]  # °C

        monthly_cooling_energy_thermal = simulation_averageday[17]*30.4  # kWh

        monthly_cooling_energy = simulation_averageday[17]*30.4  # kWh

        # Collection min and max centrifugation rate (Obsolete)

        # list centrifugation rate L.s-1
        list_centrifugation_rate_wholeunit = simulation_averageday[20]

        # list centrifugation rate L.s-1
        water_centrifuged = simulation_averageday[21]*30.4

        
        list_centrifugation_rate_wholeunit_not_0 = [
            i for i in list_centrifugation_rate_wholeunit if i != 0]
        
        min_centrifugation_rate_list.append(
            min(list_centrifugation_rate_wholeunit_not_0))

        max_centrifugation_rate_list.append(
            max(list_centrifugation_rate_wholeunit_not_0))
        
        # Temperature : Some processes have temperature as an input for their 
        # electricity consumption 
        #


        # For Mixing : Mixing is needed at any time of the day and night.
        meantemp_daytotal = (sum(collectedtemperatureevolution) / 
            len(collectedtemperatureevolution))

        monthly_Electricity_mixing_day = functions.mixing_perday(
            rhosuspension,
            tubediameter,
            pumpefficiency,
            flowrate,
            roughness,
            meantemp_daytotal,
            biomassconcentration,
            tubelength)[0] * 30.4  # MJ.m-2.month
        
        # For Drying : Drying requires to heat the slurry to 100 C and the
        # energy will depend on the initial temperature : temperature of the
        # culture at harvest time (9PM)
        
        temp_at_harvest_time = collectedtemperatureevolution[3780]




        # Summing over months in the loop
        totalproduction += monthly_production  # g dw
        
        totalproduction_harvested += monthly_production_harvested  # g dw
        
        totalproduction_loss += monthly_production_loss  # g dw

        totalcooling_thermal += monthly_cooling_energy_thermal # kWh

        totalcooling += monthly_cooling_energy  # kWh
        
        totalheating += monthly_heating_energy  # kWh
        
        total_elec_centrifuge += monthly_energy_tocentrifuge  # kWh
        
        total_elec_mixing += monthly_Electricity_mixing_day/3.6  # conversion to kWh

        totalwaterpumpedfromthefacility += monthly_waterpumped_from_the_facility  # L
        
        totalwaterpumpedfromthewell += monthly_waterpumped_from_the_well  # L

        totalwatercentrifuged += water_centrifuged  # L

        # Collecting the mean temperature over the cultivation period

        # For drying
        meantemp_at_harvest_time_cultivationperiod += temp_at_harvest_time/len(months_suitable_for_cultivation) # °C

    # End of the loop
    # Collecting min and max centrifugation rate during cultivation period (Obsolete)

    min_centrifugation_rate_m3_h = min(
        min_centrifugation_rate_list)*3.6  # m3.h-1
    
    max_centrifugation_rate_m3_h = max(
        max_centrifugation_rate_list)*3.6  # m3.h-1

    # Total production conversion to kg

    total_production_kg_dw = totalproduction/1000  # kg dw

    total_production_harvested_kg_dw = totalproduction_harvested/1000

    total_production_loss_kg_dw = totalproduction_loss/1000

    # Adding the energy for the initial heating of the well water 

    # Water of the well is heaten to Tmin
    if Twell < Tmin:
        initalheating = facilityvolume*Cp*(Tmin-Twell)/3.6  # kWh

    else: # If the water is already warm enough
        initalheating = 0  # KwH

    #Updating LCI with calculated values
    
    if heat_pump=="yes": # use COP
    
        LCIdict_updated['Heating kWh PBR'] = ((totalheating + initalheating)/total_production_harvested_kg_dw)*(
            1/bioact_molec_dbio)*(1/(1-water_after_drying))/COP 
        
        LCIdict_updated['Cooling kWh PBR'] = (
            totalcooling_thermal/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))/(COP-1)
    
    else:  #electric heater and cooling heat exchanger
        LCIdict_updated['Heating kWh PBR'] = ((totalheating + initalheating)/total_production_harvested_kg_dw)*(
            1/bioact_molec_dbio)*(1/(1-water_after_drying))
        
        LCIdict_updated['Cooling kWh PBR'] = (
            totalcooling/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))
            
    
    
    LCIdict_updated['Electricity centrifuge kWh PBR'] = (
        total_elec_centrifuge/total_production_harvested_kg_dw) * (1/bioact_molec_dbio)*(1/(1-water_after_drying))

    LCIdict_updated['Electricity mixing kWh PBR'] = (
        total_elec_mixing/total_production_harvested_kg_dw) * (1/bioact_molec_dbio)*(1/(1-water_after_drying))


    # In this version, no electricity consumption is assumed for aeration
    LCIdict_updated['Electricity aeration kWh PBR'] = 0 

    # Pumping water from the well and facility
    
    #Calling the function  for depth = well depth
    energy_perm3_fromthewell = functions.pumping_per_m3(
        rhowater, depth_well, pumpefficiency)
    
    # Pumping from the facility
    energy_perm3_fromthefacility = functions.pumping_per_m3(
        rhowater, 1, pumpefficiency)

    # Water pumped from the well
    initialpumping = facilityvolume*energy_perm3_fromthewell

    pumpingforcleaning = (cleaningvolumeVSfacilityvolume
                          * facilityvolume
                          * energy_perm3_fromthewell)

    # Energy consumption for pumping water during the cultivation
    pumping_during_cultiv = (totalwaterpumpedfromthefacility
                             * energy_perm3_fromthefacility
                             + totalwaterpumpedfromthewell
                             * energy_perm3_fromthewell)/1000  # MJ Conversion L to m3

    totalenergypumping = initialpumping+pumpingforcleaning + pumping_during_cultiv

    LCIdict_updated['Electricity pumping kWh PBR'] = ((
        (totalenergypumping/3.6)
        / total_production_harvested_kg_dw)
        * (1/bioact_molec_dbio)
        * (1/(1 - water_after_drying)))  # kWh

    # Glass consumption

    # Assuming a constant wall thickness of 2 mm.
    glass_perm2 = exchange_area * 0.002 # m3 of glass

    glass_volume_perkgmolecule = ((glass_perm2/total_production_harvested_kg_dw) 
                                  * (1/bioact_molec_dbio)
                                  * (1/(1 - water_after_drying))
                                  * 1/(glass_life_expectancy))  # m3

    glass_mass_perkgmolec = glass_volume_perkgmolecule * 2700  # kg # 2700 kg.m-3

    LCIdict_updated['Glass PBR'] = (glass_mass_perkgmolec
                                    *1/(glass_life_expectancy))

    # Drying
    
    water_to_vaporize_perkilo_dbio = ((1/slurry_concentration)
        * (1 - slurry_concentration)
        * (1 - water_after_drying))  # L. kg-1 dbio

    # /1000 for conversion from kJ to MJ and /3.6 from MJ to kWh
    Electricity_drying_perkg = (water_to_vaporize_perkilo_dbio
                                * (Cw + Cp*(100-meantemp_at_harvest_time_cultivationperiod))
                                / (boilerefficiency*1000))/3.6  # kWh.kg dbio-1

    LCIdict_updated['Electricity drying kWh PBR'] = (Electricity_drying_perkg *
        (1/bioact_molec_dbio))  # Scaled up to 1 kg of molecule kWh

    # Water consumption

    initialfilling = facilityvolume  # m3
     
    
    refillingduringcultivation = totalwaterpumpedfromthewell/1000 # m3 

    totalwater = (initialfilling 
                  + refillingduringcultivation 
                  + cleaningvolumeVSfacilityvolume*initialfilling)  # m3
    
    # Water used for cultivation and not for cleaning
    totalwater_cultivation = refillingduringcultivation + initialfilling  # m3

    totalwater_cultivation_perkgmolecule = ((totalwater_cultivation/total_production_harvested_kg_dw) 
                                            * (1/bioact_molec_dbio)
                                            * (1/(1-water_after_drying)))  # m3

    LCIdict_updated['Water(Cultivation) PBR'] = totalwater_cultivation_perkgmolecule # m3

    #Cleaning
    totalwater_cleaning_perkgmolecule = ((cleaningvolumeVSfacilityvolume*initialfilling/total_production_harvested_kg_dw) 
                                            * (1/bioact_molec_dbio)
                                            * (1/(1-water_after_drying)))

    LCIdict_updated['Water Cleaning PBR'] = totalwater_cleaning_perkgmolecule  # m3

    # Wastewater

    # All water used for cultivatiion - what has been vaporized during drying   (scaled per kg molecule)
    
    # / 1000 to convert  water_to_vaporize_perkilo_dbio from L to m3
    totalwater_towaste_perkgmolecule = (totalwater_cultivation_perkgmolecule
                                        - water_to_vaporize_perkilo_dbio
                                        * (1/bioact_molec_dbio) / 1000)  # m3

    #  Negative sign as waste treatment activity (specific to brightway)
    LCIdict_updated['Wastewater treatment PBR'] = - totalwater_towaste_perkgmolecule
    
    # Not scaled to molecule for easier wastewater concentration calculation
    totalwater_towaste = (totalwater_cultivation 
                          - water_to_vaporize_perkilo_dbio*total_production_harvested_kg_dw 
                          * (1-water_after_drying) / 1000) # m3

    # Average Concentration waste water in biomass
    conc_waste_water_biomass = total_production_loss_kg_dw/totalwater_towaste  # kg.m-3

    # kg.m-3 or g.L-1  Waste nutrient per kg molecule produced/total wastewater per kg molecule produced
    # Includes the elements in the biomass
    conc_waste_water_nutrient_N = ((N_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * N_dw)) # kg.m-3 or g.L-1
    

    conc_waste_water_nutrient_P = ((P_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * P_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_K = ((K_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * K_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_Mg =((Mg_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * Mg_dw)) # kg.m-3 or g.L-1
    
    conc_waste_water_nutrient_S = ((S_waste/totalwater_towaste_perkgmolecule 
                                    + conc_waste_water_biomass * S_dw)) # kg.m-3 or g.L-1
    

    # Carbon only in biomass, CO2 is degazed
    conc_waste_water_C = C_dw * conc_waste_water_biomass  # kg.m-3


    # Land

    LCIdict_updated['Land PBR'] = (
        1/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # m2


    # Cleaning substances

    # Half of the water with 1 substance, half with the other one
    
    #Hypochlorite
    totalhypo = ((cleaningvolumeVSfacilityvolume*facilityvolume)/2) * concentration_hypo  # kg  
    
    #Hydrogen peroxide
    totalhydro = ((cleaningvolumeVSfacilityvolume*facilityvolume)/2) * concentration_hydro  # kg   

    LCIdict_updated['Hydrogen peroxyde PBR'] = (
        totalhydro/total_production_harvested_kg_dw)*(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # kg
    
    LCIdict_updated['Hypochlorite PBR'] = (
        totalhypo/total_production_harvested_kg_dw) *(1/bioact_molec_dbio)*(1/(1-water_after_drying))  # kg


    #Extraction and substitution

    if extraction == 'yes':
        # 1 kWh to disrupt 1 kg of microalgal biomass (kg dbio)
        LCIdict_updated['Electricity cell disruption kWh PBR'] = 1 * (1/bioact_molec_dbio)  # kWh.kg-1

        # for anaerobic digestion, we keep the wet biomass 
        bioact_molec_dw = bioact_molec_dbio/(1-water_after_drying)

        # If extraction, then the 
        # remaining biomass composition is changed according to the biochemical class of the extracted molecule

        if biochemicalclass == 'lip':
            
            lip_dbio_after_extract = (lip_dbio - bioact_molec_dbio)/(1-bioact_molec_dbio)
            
            lip_dw_after_extract = (lip_dw - bioact_molec_dw)/(1-bioact_molec_dw)
          
        
            carb_dbio_after_extract = carb_dbio / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = carb_dw /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = prot_dbio / (1-bioact_molec_dbio)
            
            prot_dw_after_extract = prot_dw /(1-bioact_molec_dw)


            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)


            water_dbio_after_extract = water_after_drying /(1-bioact_molec_dbio)
            


        elif biochemicalclass == 'carb':
            
            lip_dbio_after_extract = lip_dbio / (1-bioact_molec_dbio)
            
            lip_dw_after_extract = lip_dw/(1-bioact_molec_dw)

            
            carb_dbio_after_extract = (carb_dbio-bioact_molec_dbio) / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = (carb_dw-bioact_molec_dw) /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = prot_dbio / (1-bioact_molec_dbio)
           
            prot_dw_after_extract = prot_dw /(1-bioact_molec_dw)

            
            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)
            
            
            water_dbio_after_extract = water_after_drying / (1-bioact_molec_dbio)

        elif biochemicalclass == 'prot':
            
            lip_dbio_after_extract = lip_dbio / (1-bioact_molec_dbio)

            lip_dw_after_extract = lip_dw/(1-bioact_molec_dw)
            
            
            carb_dbio_after_extract = carb_dbio / (1-bioact_molec_dbio)
            
            carb_dw_after_extract = carb_dw /(1-bioact_molec_dw)

            
            prot_dbio_after_extract = (prot_dbio-bioact_molec_dbio) / (1-bioact_molec_dbio)
            
            prot_dw_after_extract = (prot_dw-bioact_molec_dw) /(1-bioact_molec_dw)


            ash_dbio_after_extract = ash_dbio / (1-bioact_molec_dbio)
            
            ash_dw_after_extract = ash_dw /(1-bioact_molec_dw)
           
            water_dbio_after_extract = water_after_drying / (1-bioact_molec_dbio)


        # After extraction, the substitution will occur with the new composition of the biomass
        
        # Call the function which calculates the masses of subsituted fish feed ingredient
        substitution = functions.optimization_for_fishfeed_substitution(fishfeed_table,
                                                                        lip_dbio_after_extract,
                                                                        prot_dbio_after_extract,
                                                                        carb_dbio_after_extract,
                                                                        water_dbio_after_extract, 
                                                                        ash_dbio_after_extract,
                                                                        incorporation_rate,
                                                                        MJ_kgcarb,
                                                                        MJ_kgprot,
                                                                        MJ_kglip)

    # if no extraction (molecule given to fish directly), 
    #  the biomass composition stays the same (Obsolete)
    else:  
        LCIdict_updated['Electricity cell disruption kWh PBR'] = 0  # kWH.kg-1
        LCIdict_updated['Extraction electricity kWh PBR'] = 0
        LCIdict_updated['Co solvent Extraction PBR'] = 0

        substitution = functions.optimization_for_fishfeed_substitution(fishfeed_table,
                                                                        lip_dbio,
                                                                        prot_dbio,
                                                                        carb_dbio,
                                                                        water_after_drying,
                                                                        ash_dbio,
                                                                        incorporation_rate,
                                                                        MJ_kgcarb,
                                                                        MJ_kgprot,
                                                                        MJ_kglip)

    # Choose the market that the dependent coproducts enter

    if market_for_substitution == 'animal feed':  

        # Model substitution 1  Animal Feed

        # kg #the same subsitution occurs for every kilo
        feedprot_m1 =  substitution[0] * (1/bioact_molec_dbio - 1) # kg
        feedenergy_m1 =  substitution[1] * (1/bioact_molec_dbio - 1)  # MJ

        LCIdict_updated['Feed energy PBR'] = -feedenergy_m1
        LCIdict_updated['Feed protein PBR'] = -feedprot_m1

        optimization_performance = 'No optimization'
        substitution_check = 'No optimization'

    # Model substitution 2 Fish Feed
    
    elif market_for_substitution == 'fish feed':

        LCIdict_updated['Feed energy PBR'] = 0  # Do not use Model 1
        LCIdict_updated['Feed protein PBR'] = 0

        if extraction == 'yes':  # Then the substituion only takes place with the remaining biomass
        
            # vect_substitution is a list containing the masses of fish feed ingredient substituted by the remaining biomass
            # (1/bioact_molec_dbio-1) = remaining biomass after extraction of the FU : 1kg of molecule
            
            # substitution[5] is the list of masses of fish feed ingredient 
            # replaced by the given biomass composition. in kg. kg dbio-1 (sum =1 kg)

            vect_substitution = substitution[5]*(1/bioact_molec_dbio-1)  # kg

            # evaluation of the optimized recipe
            optimization_performance = substitution[-1]*(1/bioact_molec_dbio-1)

        else:  # then the molecule incorporated in the biomass takes part in the substutition

            # substitution[5] is the list of masses of fish feed ingredient
            # replaced by the given biomass composition. in kg. kg dbio-1 (sum =1 kg)
            vect_substitution = substitution[5]*(1/bioact_molec_dbio)  # kg
            
            # evaluation of the optimized recipe
            optimization_performance = substitution[-1]*(1/bioact_molec_dbio)

        # Adding the ingredients of the fish feed to the LCI for substitution

        substitution_check = 0 # Obsolete
        
        for a in range(0, len(fishfeed_table['Ingredient'])):
            
            # Ingredients are ranked in the same order in the vector and in the fish feed table
            LCIdict_updated[fishfeed_table['Ingredient'][a]] = vect_substitution[a]
            
            substitution_check = (substitution_check
                                  + LCIdict_updated[fishfeed_table['Ingredient'][a]]) #Obsolete


    # biogaz

    elif market_for_substitution == 'anaerobic digestion':

        LCIdict_updated['Feed energy PBR'] = 0  # Do not use Model 1
        LCIdict_updated['Feed protein PBR'] = 0
        LCIdict_updated["Electricity drying kWh PBR"] = 0 # No drying, directly on wet biomass

        remaining_biomass_after_extraction = 1/bioact_molec_dw-1 # kg

        optimization_performance = 'No optimization'
        substitution_check = 'No optimization'


        [N_substituted,
         P_substituted,
         Mg_substituted,
         K_substituted,
         MJ_contained_in_produced_biogas]=functions.anaerobic(prot_dw_after_extract,
                                                    lip_dw_after_extract,
                                                    carb_dw_after_extract,
                                                    ash_dw_after_extract,
                                                    fertiliser_substitutability_AD,
                                                    carbon_degradibilidy_AD,
                                                    elemental_contents,
                                                    phospholipid_fraction,
                                                    CH4_LHV)
                                                    
                                                
                                                
        LCIdict_updated["N fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * N_substituted # as N      
        LCIdict_updated["P fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * P_substituted/0.4366      
        LCIdict_updated["Mg fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * Mg_substituted * (120.4/24.3)
        LCIdict_updated["K fertiliser substitution AD PBR"]   = -remaining_biomass_after_extraction * K_substituted * 1.2   
        
        
        # Inputs to AD

        
        LCIdict_updated["MJ heat biogas AD PBR"] = remaining_biomass_after_extraction * (heat_input_per_kg_biomass_AD - MJ_contained_in_produced_biogas) 
        
        
        LCIdict_updated["Electricity kWh AD PBR"] = elec_input_per_kg_biomass_AD * remaining_biomass_after_extraction   

        # Combustion of the methane used for AD leads to CO2 emission captured and injected in culture

        total_kg_methane_burnt_for_AD = remaining_biomass_after_extraction * heat_input_per_kg_biomass_AD/(CH4_LHV/3.6)  # kg CH4. kg molecule -1
        
       
        total_moles_methanes_burnt_for_AD = total_kg_methane_burnt_for_AD/(16*10**-3)   # moles CH4 burnt.kg molecule-1
        
        # mole ratio = volume ratio for perfect gases  (same molar volumes)
        # biogas from ecoinvent has 0,64 m3 CH4 per m3 biogas

        total_moles_CO2_in_biogas_for_AD = total_moles_methanes_burnt_for_AD*0.36/(0.64)  # moles CO2 burnt.kg molecule-1
        
        # all moles of CH4 become moles of CO2
        total_kg_CO2_resulting_from_burning_biogas_for_AD =( total_moles_methanes_burnt_for_AD + total_moles_CO2_in_biogas_for_AD)*44*10**-3 # 44 g.mol CO2 -1
        
        
        
        
        LCIdict_updated["CO2 substitution biogas burning AD PBR"] = - total_kg_CO2_resulting_from_burning_biogas_for_AD

 
    # Yields

    numberofcultivationdays = len(months_suitable_for_cultivation)*30.4 # days
    
    volumetricyield = (total_production_kg_dw 
                       / (facilityvolume*1000*numberofcultivationdays))  # kg.L-1.d-1
    
    surfaceyield = total_production_kg_dw/numberofcultivationdays  # kg.days-1


    # Check mass balance (Obsolete)

    needed_dbio_check = 1/(lipid_af_dw * (1-ash_dw) *
                           (1-water_after_drying) * bioact_molec_dbio)
    
    
    
    
    
    
    

    return [LCIdict_updated,
            surfaceyield,
            volumetricyield,
            optimization_performance,
            needed_dbio_check,
            substitution_check,
            total_production_kg_dw,
            total_production_harvested_kg_dw,
            total_production_loss_kg_dw,
            conc_waste_water_nutrient_N,
            conc_waste_water_nutrient_P,
            conc_waste_water_nutrient_K,
            conc_waste_water_nutrient_Mg,
            conc_waste_water_biomass,
            conc_waste_water_C,
            conc_waste_water_nutrient_S,
            bioact_molec_dbio,
            min_centrifugation_rate_m3_h,
            max_centrifugation_rate_m3_h,
            totalwatercentrifuged,
            tubelength,
            facilityvolume,
            exchange_area,
            totalcooling_thermal,
            Twell]  # []



def sampling_func(Tech_opdict_distributions,
                  Biodict_distributions,
                  Locationdict_distributions, 
                  Physicdict_distributions,
                  size, 
                  type_sens):
    '''Function which returns a random sample for the input space of 
    non-constant parameters according to the Sensitivity analysis alorgithm (Sobol or Fast)
    
    Inputs:
        
        # All parameters distributions dictionnaries : 
                  Tech_opdict_distributions,
                  Biodict_distributions,
                  Locationdict_distributions, 
                  Physicdict_distributions
                  
        #size : Size of the sample 
        ( Final number of combiantions =size*number uncertain parameters for FAST,
         size*(number uncertain parameters+2) for Sobol)
        
        #type_sens: "SOBOL" or "FAST"
    
    Outputs :
        
        -sample :  Generated sample. Array with 1 row = 1 combination of uncertain parameters
        -names_param :  List of names of the uncertain parameters
        -names_param_op :  List of names of the uncertain parameters from Tech_opdict_distributions
        -names_param_bio :  List of names of the uncertain parameters from Biodict_distributions
        -names_param_geo :  List of names of the uncertain parameters from Locationdict_distributions
        -names_param_phy :  List of names of the uncertain parameters from Physicdict_distributions
        -problem :  Problem format for the sensitivity analysis from Salib library (Saltelli)
    '''

    # Creation of the problem

    # Creating a copy of the distribtutions dictionnaires containing only the
    # parameters with distributions(not the constant ones)

    # Tech_opdict
    Tech_opdict_distributions_input = Tech_opdict_distributions.copy()
    
    to_be_deleted_op = []
    
    # Collecting the parameters that are constant
    
    for param in Tech_opdict_distributions_input:
        if Tech_opdict_distributions_input[param][0] == 'unique' or Tech_opdict_distributions_input[param][0] == 'binary':
            to_be_deleted_op.append(param)
    
    # Deleting the parameters that are constant
    
    for a in to_be_deleted_op:
        Tech_opdict_distributions_input.pop(a)

    # Biodict
    #print('5555555555')
    Biodict_distributions_input = Biodict_distributions.copy()

    to_be_deleted_bio = []

    for param in Biodict_distributions_input:
        if Biodict_distributions_input[param][0] == 'unique' or Biodict_distributions_input[param][0] == 'binary':
            to_be_deleted_bio.append(param)

    for a in to_be_deleted_bio:
        Biodict_distributions_input.pop(a)

    # Geography

    Locationdict_distributions_input = Locationdict_distributions.copy()

    to_be_deleted_geo = []

    for param in Locationdict_distributions_input:
        if Locationdict_distributions_input[param][0] == 'unique' or Locationdict_distributions_input[param][0] == 'binary':
            to_be_deleted_geo.append(param)

    for a in to_be_deleted_geo:
        Locationdict_distributions_input.pop(a)

    # Physics

    Physicdict_distributions_input = Physicdict_distributions.copy()

    to_be_deleted_phy = []

    for param in Physicdict_distributions_input:
        if Physicdict_distributions_input[param][0] == 'unique' or Physicdict_distributions_input[param][0] == 'binary':
            to_be_deleted_phy.append(param)

    for a in to_be_deleted_phy:
        Physicdict_distributions_input.pop(a)

    # Collecting names, bounds , dists to create Saltelli problem
    names_param = []
    bounds = []
    dists = []
    






    # 1) Operational
    
    names_param_op = []

    count_col = -1
    col_with_triang = [] # Triangular distributions will need post processing
    actual_lower_bounds = []

    for param in Tech_opdict_distributions_input:
        
        count_col += 1
        
        names_param_op.append(param)
        
        distrib = Tech_opdict_distributions_input[param][0]
        
        dists.append(distrib)
        #print(distrib)
        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Tech_opdict_distributions_input[param]
                           [1][1], Tech_opdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Tech_opdict_distributions_input[param]
                           [1][3], Tech_opdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0

            width = (Tech_opdict_distributions_input[param][1][2] 
                     - Tech_opdict_distributions_input[param][1][1])
            
            peak_lowerbound = (Tech_opdict_distributions_input[param][1][3] 
                               - Tech_opdict_distributions_input[param][1][1])
            
            bounds.append([width, peak_lowerbound/width])

            # Collect column with triang distribution to shift the values eventually (otherwise lower bound=0)
            col_with_triang.append(count_col)
            
            actual_lower_bounds.append(Tech_opdict_distributions_input[param][1][1])

    # 2) Biological

    names_param_bio = []

    for param in Biodict_distributions_input:

        count_col += 1
        
        names_param_bio.append(param)
        
        distrib = Biodict_distributions_input[param][0]
        
        dists.append(distrib)
        
        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Biodict_distributions_input[param][1]
                           [1], Biodict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Biodict_distributions_input[param][1]
                           [3], Biodict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0

            width = (Biodict_distributions_input[param][1][2] 
                     - Biodict_distributions_input[param][1][1])
            
            peak_lowerbound = (Biodict_distributions_input[param][1][3] 
                               - Biodict_distributions_input[param][1][1])
            
            bounds.append([width, peak_lowerbound/width])

            # Collect column with triangle distribution to shift the values eventually (otherwise lower bound=0)
            col_with_triang.append(count_col)
            
            actual_lower_bounds.append(Biodict_distributions_input[param][1][1])

    # 3) Geography
    
    names_param_geo = []

    for param in Locationdict_distributions_input:

        count_col += 1

        names_param_geo.append(param)

        distrib = Locationdict_distributions_input[param][0]

        dists.append(distrib)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Locationdict_distributions_input[param]
                           [1][1], Locationdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Locationdict_distributions_input[param]
                           [1][3], Locationdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0

            width = (Locationdict_distributions_input[param][1][2] 
                     - Locationdict_distributions_input[param][1][1])
            
            peak_lowerbound = (Locationdict_distributions_input[param][1][3]
                               - Locationdict_distributions_input[param][1][1])

            bounds.append([width, peak_lowerbound/width])

            # Collect column with triang distribution to shift the values eventually (otherwise lower bound=0)
            col_with_triang.append(count_col)
            
            actual_lower_bounds.append(Locationdict_distributions_input[param][1][1])

    # 3) Physics
    
    names_param_phy = []

    for param in Physicdict_distributions_input:
        
        distrib = Physicdict_distributions_input[param][0]
        
        dists.append(distrib)

        count_col += 1
        
        names_param_phy.append(param)

        if distrib == 'unif':  # then bounds = upper bound, lower bound

            bounds.append([Physicdict_distributions_input[param]
                           [1][1], Physicdict_distributions_input[param][1][2]])

        elif distrib == 'norm':  # then bounds = mean, sd

            bounds.append([Physicdict_distributions_input[param]
                           [1][3], Physicdict_distributions_input[param][1][4]])

        elif distrib == 'triang':  # then bounds = width, location of the peak in % of width, Assume lower bound = 0

            width = (Physicdict_distributions_input[param][1][2] 
                     - Physicdict_distributions_input[param][1][1])
            
            peak_lowerbound = (Physicdict_distributions_input[param][1][3] - Physicdict_distributions_input[param][1][1])
            
            bounds.append([width, round(peak_lowerbound/width,2)])

            # Collect column with triang distribution to shift the values eventually (otherwise lower bound=0)
            col_with_triang.append(count_col)

            actual_lower_bounds.append(Physicdict_distributions_input[param][1][1])

    names_param = names_param_op + names_param_bio + names_param_geo + names_param_phy
    
    
    
    #print to be deleted
   # for i in range(len(bounds)):
        #print(i)
        #print(dists[i],bounds[i])
    
    #print(dists)
    problem = {'num_vars': len(names_param),  # number of variables
               'names': names_param,  
               'bounds': bounds,
               'dists': dists}

    #print(problem)
    if type_sens == 'SOBOL':
        sample = SALib.sample.saltelli.sample(problem,
                                              size,
                                              calc_second_order=False)
    if type_sens == 'FAST':
        sample = SALib.sample.fast_sampler.sample(problem,
                                                  size)

    # Shift the values for the triangular distributions, otherwise lowerbound=0
    for index_col in range(len(col_with_triang)):

        sample[:, col_with_triang[index_col]] = sample[:,col_with_triang[index_col]]+actual_lower_bounds[index_col]

    return sample, names_param, names_param_op, names_param_bio, names_param_geo, names_param_phy, problem






@ray.remote
def calculateLCA_1param_parallel(constant_inputs,param_set) :   
        """Function which calculates one LCA and returns the LCIA results, given a set of input parameters.
        The function is built to be called in parrallel with ray."""      
    
        (Tech_opdict,
         Biodict,
         Locationdict,
         Physicdict,
         methods,
         list_cfs,
         dict_mono_technosphere_lcas,
         processes_in_categories,
         LCIdict,
         months_suitable_for_cultivation,
         fraction_maxyield, 
         fishfeed_table,
         elemental_contents,
         names_param_op,
         names_param_bio,
         names_param_geo,
         names_param_phy,
         names_param)=constant_inputs    
    
        # Copy that will be modified for calculation of the LCA for this set of parameters
        
        #print("processes_in_categories",processes_in_categories)
        
        
        new_dict_mono_technosphere_lcas = copy.deepcopy(dict(dict_mono_technosphere_lcas))

    
        # Update the dictionnaries whith the values of the sample

        # 1 ) Tech_opdict
        
        for param in Tech_opdict:  
            # We browse the parameters to look for the uncertain ones
            # which need to be updated

            for index in range(len(names_param_op)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:
                    
                    # Then it is an uncertain paramater and its value is 
                    # updated with the one from the generated sample
                    
                    Tech_opdict[param] = param_set[index]
                    
        # We will do the same thing for other dictionnaries but there is no need to browse all possible parameters, 
        # just the ones from other dictionnaries which are left.
        
        new_start = len(names_param_op)

        # 2) Biodict

        for param in Biodict:  
            # We browse the parameters to look for the variable ones
            # which need to be updated

            for index in range(new_start, new_start+len(names_param_bio)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Biodict[param] = param_set[index]

        new_start = new_start+len(names_param_bio)

        # 3) Locationdict

        for param in Locationdict:  
            # We browse the parameters to look for the variable ones
            # that need to be updated
            for index in range(new_start, new_start+len(names_param_geo)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Locationdict[param] = param_set[index]

        new_start = new_start+len(names_param_geo)
        
        # 4) Physicdict
        
        for param in Physicdict: 
            # We browse the parameters to look for the variable ones
            # that need to be updated

            for index in range(new_start, new_start+len(names_param_phy)):

                # Looking for the corresponding parameter in the saltelli set
                if names_param[index] == param:

                    Physicdict[param] = param_set[index]

        


        # Calculates LCI for this set with the updated dictionnaries
        
        LCI = LCI_one_strain_uniquevalues(Biodict,
                                          Physicdict,
                                          Tech_opdict,
                                          Locationdict,
                                          LCIdict,
                                          months_suitable_for_cultivation,
                                          fraction_maxyield, 
                                          fishfeed_table,
                                          elemental_contents)
        
        # Collecting the results of the function
        LCIdict_collected = LCI[0]


        surfaceyield = LCI[1]  # kg dw .d-1
        
        volumetricyield = LCI[2]  # kg dw.L-1.d-1
        
        optimization_performance = LCI[3]  # kg dw.d-1
        
        total_production_kg_dw = LCI[6]  # kg dw

        total_production_harvested_kg_dw = LCI[7]  # kg dw

        total_production_loss_kg_dw = LCI[8]  # kg dw

        conc_waste_water_nutrient_N = LCI[9]  # kg.m-3

        conc_waste_water_nutrient_P = LCI[10]  # kg.m-3
        
        conc_waste_water_nutrient_K = LCI[11]  # kg.m-3

        conc_waste_water_nutrient_Mg = LCI[12]  # kg.m-3

        conc_waste_water_biomass = LCI[13]  # kg.m-3

        conc_waste_water_C = LCI[14]  # kg.m-3

        conc_waste_water_nutrient_S = LCI[15]  # kg.m-3

        bioact_molec_dbio = LCI[16]  # .

        min_centrifugation_rate_m3_h = LCI[17]

        max_centrifugation_rate_m3_h = LCI[18]

        totalwater_centrifuged = LCI[19]

        totalwater_centrifuged = LCI[19]

        tubelength = LCI[20]

        facilityvolume = LCI[21]

        exchange_area = LCI[22]

        totalcooling_thermal = LCI[23]

        Twell = LCI[24]

        # List containing simualtions values which are not LCI or LCIA
        # (same order as their names in names_suppl_info)
        
        values_simu = [bioact_molec_dbio,
                       surfaceyield,
                       tubelength,
                       facilityvolume,
                       exchange_area,
                       totalcooling_thermal,
                       volumetricyield,
                       total_production_kg_dw,
                       total_production_harvested_kg_dw,
                       total_production_loss_kg_dw,
                       conc_waste_water_nutrient_N,
                       conc_waste_water_nutrient_P,
                       conc_waste_water_nutrient_K,
                       conc_waste_water_nutrient_Mg,
                       conc_waste_water_biomass,
                       conc_waste_water_C,
                       conc_waste_water_nutrient_S,
                       min_centrifugation_rate_m3_h,
                       max_centrifugation_rate_m3_h,
                       totalwater_centrifuged,
                       Twell]

        values_LCI = [LCIdict_collected[i] for i in LCIdict_collected]

        names_LCI = [a for a in LCIdict_collected]
        

        # Now calculating the LCIA for this row
        
        # Calling the fuction that calculates the impacts associated to the emissions 
        # of the wastewater treatment of 1 cubic meter of the wastewater

        list_sum_impacts_biosphere_waste_water= waste_water_impact_biosphere(
                      conc_waste_water_nutrient_N,
                      conc_waste_water_nutrient_P,
                      conc_waste_water_C,
                      conc_waste_water_nutrient_Mg,
                      conc_waste_water_nutrient_K,
                      conc_waste_water_nutrient_S,
                      methods,
                      list_cfs)
        
        # Adding the biosphere flows for 1 cubic meter of wastewater
        
        new_dict_mono_technosphere_lcas['Wastewater treatment PBR'] =[a+b for (a,b) in zip(new_dict_mono_technosphere_lcas['Wastewater treatment PBR'],list_sum_impacts_biosphere_waste_water) ]

        # Multipliying each impact per unit of input processes by the input amount in the calculated LCI
        
        for i in new_dict_mono_technosphere_lcas:
            
            new_dict_mono_technosphere_lcas[i]=[LCIdict_collected[i]*a for a in new_dict_mono_technosphere_lcas[i]]
        
        # Calculating total LCA by summing
        
        
        list_LCA_res =[]
        
        
        for meth_index in range(len(methods)):
            
            sum_impact = sum([new_dict_mono_technosphere_lcas[flow][meth_index] for flow in new_dict_mono_technosphere_lcas ])


            list_LCA_res.append(sum_impact)
            
            
            
        # Row to add to the results dataframe (Uncertain parameters values, LCI figures, Other values about the simulation, LCIA)
       

        row_to_add = list(param_set) + values_LCI + values_simu + list_LCA_res

        # Contributions

        rows_table_contributions =[[0]*len(processes_in_categories) for n in range(len(methods))]
        

        
        
        #print("new_dict_mono_technosphere_lcas",new_dict_mono_technosphere_lcas)
        
        # Contribution per process category
        for process in new_dict_mono_technosphere_lcas :
            
            # browsing the categories
            for index_content_categ in range(len(processes_in_categories)):
                
    
                # if this process belongs to category
                if process in processes_in_categories[index_content_categ]:
                    
                    # Then we add this value to the corresponding column in the  list_tables_contribution
                    for meth_index in range(len(methods)): #we do this for all methods 
                        
                        
                        rows_table_contributions[meth_index][index_content_categ] = (
                            rows_table_contributions[meth_index][index_content_categ] 
                            + new_dict_mono_technosphere_lcas[process][meth_index])    

        
        return (row_to_add, rows_table_contributions, names_LCI)     








def final_function_simulations(dict_mono_technosphere_lcas,
                               Tech_opdict,
                               Biodict,
                               Locationdict,
                               Physicdict,
                               Tech_opdict_distributions,
                               Biodict_distributions,
                               Locationdict_distributions,
                               Physicdict_distributions,
                               LCIdict,
                               size,
                               months_suitable_for_cultivation,
                               fraction_maxyield,
                               elemental_contents, 
                               fishfeed_table, 
                               methods,
                               list_cfs,
                               categories_contribution, 
                               processes_in_categories,
                               type_sens):
    '''Function which calls all other functions and generates the LCA results, 
    uncertainty, sensitivity and contribution analysis.
    
    Inputs:

        #All normal parameters dictionnaries:
            -Tech_opdict
            -Biodict
            -Locationdict
            -Physicdict
            
        # All parameters distributions dictionnaries : 
            -Tech_opdict_distributions,
            -Biodict_distributions,
            -Locationdict_distributions, 
            -Physicdict_distributions
            
        #LCIdict: the  LCI dictionnary
        #size : Size of the sample 
        ( Final number of combiantions =size*number uncertain parameters for FAST,
         size*(number uncertain parameters+2) for Sobol)
        
        #months_suitable_for_cultivation : Months for cultivation ; 
        list of month numbers : [a,b,c]
        
        #fraction_maxyield : Fraction of the maximum yield achieved ; .
        #elemental_contents : DataFrame with elemental compositons of macronutrients
        #fishfeed_table : DataFrame with fish feed composition
        #methods : List of Impact categories to apply for the LCA
        #categories_contribution : Names of process categories considered for the contribution analysis
        #processes_in_categories : List of processes to assign to categories (same order)
        #type_sens: "SOBOL" or "FAST"
    
    
    Outputs :

        #sample :  Randomly generated sample. Array with 1 row = 1 combination of uncertain parameters (1 iteration)
        #problem_sobol_FAST :  Sobol or Fast problem as generated by SAlib
        #results_table_df : Final dataframe with 1 row per simulation (iteration).
        Each row contains the values of the uncertain parameters, the LCI figures,
        key figures about the simulation (yields etc.) and the LCIA.
        
        #results_table :  Same but as an numpy array
        #results_sobol_fast :  List of lists of type : [name impact category, Dataframe with indexes for each index] 
        sensitivity indexes for uncertain parameters for one impact category 
        
        #list_tables_contribution_df_melted :  List of melted dataframes with 
        each dataframe containing contributions for each processes for one impact category 
        
        #list_tables_contribution_abs_df_melted  : List of melted dataframes with 
        each dataframe containing contributions for each processes for one impact category 
        Contribution calculated as share oh the sum of absolute values.
        
        #all_methods_contribution_abs_df_melted : Dataframes merging the previous ones (all impact categories)
            
        #list_tables_contribution_df :  Same as previous but not melted

        #list_opti_perfo : list of performances obtained by the optimization 
        algorithnm for each iteration. (Obsolete)
        
        #result_LCI_rescaled :  Dataframe with 1 row per simulation (iteration) 
        but LCI figures are scaled to 1kg of dried biomass instead of 1 kg of molecule
        
        #sensi_multi_melt : Dataframe with Total_order Sensitivity index for each parameter and each impact category
        #desc_stat_results : Dataframe with statitstical description of the simulation's outputs.
        #total_desc_stat_contri_df : Dataframe with statitstical description of the contribution analysis.
    '''

    # Columns that should not be processed numerically
    
    columns_not_float = ['Nsource', 'night_monitoring',
                         'Bio_class', 'market_for_substitution']
    
    # Columns that are not LCI nor LCIA but other simulation's outputs
    
    names_suppl_info = ['bioact_molec_dbio',
                    'Areal productivity kg.m-2.d-1',
                    'tube length m',
                    'PBR volume m3',
                    'exchange area m2',
                    'total cooling (thermal kWh)',
                    'Volumetric productivity kg.L-2.d-1',
                    'Total production kg dw.m-2',
                    'Total production harvested kg dw.m-2',
                    'Total_production_loss_kg dw.m-2',
                    'Conc wastewater g N.L-1',
                    'Conc wastewater g P.L-1',
                    'Conc wastewater g K.L-1',
                    'Conc wastewater g Mg.L-1',
                    'Conc wastewater g dw.L-1',
                    'Conc wastewater g C.L-1',
                     'Conc wastewater g S.L-1',
                    ' min_centrifugation_rate_m3_h',
                    'max_centrifugation_rate_m3_h',
                    'Total volume centrifuged m3',
                    'Twell']

    # Creates a sample and a saltelli problem with the function
    
    sampling_res = sampling_func(Tech_opdict_distributions,
                                 Biodict_distributions,
                                 Locationdict_distributions,
                                 Physicdict_distributions,
                                 size,
                                 type_sens)

    sample = sampling_res[0]

    names_param = sampling_res[1]

    names_param_op = sampling_res[2]

    names_param_bio = sampling_res[3]

    names_param_geo = sampling_res[4]

    names_param_phy = sampling_res[5]

    problem_sobol_FAST = sampling_res[6]


    # Initialize variables that will receive results
    
    results_table = np.empty((0,
                              len(names_param)+len(LCIdict) +len(names_suppl_info)+len(columns_not_float)+len(methods)),
                             dtype=float)

    # list of tables whih will contain the conribution of each process category to each impact category
    
    list_tables_contribution = [np.zeros((len(sample),
                                          len(categories_contribution)),
                                         dtype=float) for i in range(len(methods))]  
    
    # Contributions calculated by dividing by the sum of the absolute values
    # We will only keep this one eventually
    list_tables_contribution_abs=[np.zeros((len(sample),
                                          len(categories_contribution)),
                                         dtype=float) for i in range(len(methods))] 

    
    
    
    
    
    
    
    
    
    names_methods_adjusted = [a[-1] for a in methods]
        
   
    # Put the constant objects into ray.put()
    
    
    constant_inputs = ray.put([Tech_opdict,
                               Biodict,
                               Locationdict,
                               Physicdict,
                               methods,
                               list_cfs,
                               dict_mono_technosphere_lcas,
                               processes_in_categories,
                               LCIdict,
                               months_suitable_for_cultivation,
                               fraction_maxyield, 
                               fishfeed_table,
                               elemental_contents,
                               names_param_op,
                               names_param_bio,
                               names_param_geo,
                               names_param_phy,
                               names_param]) 


    
    
    #  Parallel    
    
    
    arrayresult_raw =ray.get([calculateLCA_1param_parallel.remote(constant_inputs,
                                                                  param_set) for param_set in sample])        
            

    
        

     

        
    # Names LCI
    
    # for instance
    names_LCI = arrayresult_raw[0][2]
    
    
    # Rebuild table of results
    results_table = np.array([row[0] for row in arrayresult_raw])
    
    names_for_df = names_param + names_LCI + names_suppl_info + names_methods_adjusted

        
    results_table_df = pd.DataFrame(results_table, columns=names_for_df)

    # Rebuild list of tables of contributions
        
    intermediary_list_contrib = [row[1] for row in arrayresult_raw]
        
    

    # Build list of table for contribution
        
    list_tables_contribution=[np.array([intermediary_list_contrib[i][b] for i in range(len(intermediary_list_contrib))])
                                  for b in range(len(intermediary_list_contrib[0]))]    
    
    
    
        
  

    #       
        

        

    del results_table # Save memory
    # Calculating % contribution
    
    
    
    #Calulating contribution sum
    for index_method in range(len(methods)):

        for index_row in range(len(sample)):

            
            sumrow_abs = sum([abs(a) for a in list_tables_contribution[index_method][index_row,:]])
            
            for index_col in range(len(categories_contribution)):
                
                list_tables_contribution_abs[index_method][index_row][index_col] =(
                    list_tables_contribution[index_method][index_row][index_col]
                    *100/sumrow_abs)
                




    #  Conversion to Dataframes
    list_tables_contribution_df = [pd.DataFrame(
        table, columns=categories_contribution) for table in list_tables_contribution]

    list_tables_contribution_abs_df=[pd.DataFrame(
        table, columns=categories_contribution) for table in list_tables_contribution_abs]

    
    del list_tables_contribution_abs # Save memory
    
    del list_tables_contribution # Save memory
    
    # Statistical description of contributions
    
    total_desc_stat_contri_df =pd.DataFrame()
    
    # For all methods
    for index_meth in range(len(methods)):
        
        # statisitical description of the contributions

        desc_stat_contrib = list_tables_contribution_abs_df[index_meth].describe()
   
        # add a speration row to start a new method
        separation_row = {i:'' for i in desc_stat_contrib.columns}
        
        #  Change First key of the new row to name of method
        firstkey = list(separation_row.keys())[0]

        separation_row[firstkey]=methods[index_meth][-1]
        
        # add this separation row nrow
         
        total_desc_stat_contri_df = total_desc_stat_contri_df.append(separation_row, ignore_index=True)

        # Add the statistical description
        
        total_desc_stat_contri_df=pd.concat([total_desc_stat_contri_df,
                                             desc_stat_contrib])
    
    

    # Sensitivity

    # Cleaning Dataframe, change object type to float except for qualitative variables

    columnnames_without_quali = [a for a in names_for_df if a not in columns_not_float]

    results_table_df[columnnames_without_quali] = results_table_df[columnnames_without_quali].astype(float)

    
    # Calculating sensitivity indexes for each method
    results_sobol_fast = []

    for IC in methods:  # For each column of LCIA results = for each method (Impact category)

        output = results_table_df.loc[:, IC[-1]] # Calls the result column for this impact category

        # Rearranging into a monodimensional array
        
        array_output = pd.DataFrame(output).to_numpy()

        flat_array = array_output.flat

        output_list_clean = []

        for a in flat_array:
            output_list_clean.append(a)

        output_clean = np.array(output_list_clean)

        # Performing the sensitivy analysis

        if type_sens == 'SOBOL':
            sobol_res = SALib.analyze.sobol.analyze(
                problem_sobol_FAST, output_clean, calc_second_order=False)
        elif type_sens == 'FAST':
            sobol_res = SALib.analyze.fast.analyze(
                problem_sobol_FAST, output_clean)

        results_sobol_fast.append([IC[-1], sobol_res])



     
     
    # Sensitivity
    
    # Inititalizing table for sensitivy indices for all methods
     
    name_columns=[IC[0]+'__'+ type_ind  for IC in results_sobol_fast for type_ind in list(results_sobol_fast[0][1])]                

    count=-1
    if type_sens == 'SOBOL':
        
        name_columns=[IC[0]+'__'+ type_ind  for IC in results_sobol_fast for type_ind in list(results_sobol_fast[0][1])]                

        sensi_multi = pd.DataFrame(np.zeros((len(results_sobol_fast[0][1]['ST']),
                                             len(results_sobol_fast)*4)),
                     columns=name_columns)
        
            # put values in datframes (1 column = 1 Impact category, n rows=n paramters)

        
        for IC in results_sobol_fast :
            count += 1
     
            sensi_multi.iloc[:,count*4]=IC[1]['S1']
            sensi_multi.iloc[:,count*4+1]=IC[1]['S1_conf']
            sensi_multi.iloc[:,count*4+2]=IC[1]['ST']
            sensi_multi.iloc[:,count*4+3]=IC[1]['ST_conf']     
        
    elif type_sens == 'FAST':
        
        name_columns=[IC[0]+'__'+ type_ind  for IC in results_sobol_fast for type_ind in list(results_sobol_fast[0][1])[:-1]]                

        
        sensi_multi = pd.DataFrame(np.zeros((len(results_sobol_fast[0][1]['ST']),
                                     len(results_sobol_fast)*2)),
             columns=name_columns)
        
        # put values in datframes (1 column = 1 Impact category, n rows=n paramters)

        
        for IC in results_sobol_fast :
            count += 1
     
            sensi_multi.iloc[:,count*2]=IC[1]['S1']
            sensi_multi.iloc[:,count*2+1]=IC[1]['ST']
            
    
    
  
    
    # put values in datframes (1 column = 1 Impact category, n rows=n paramters)
        

    sensi_multi['Parameter'] = names_param

    
    sensi_multi_melt=pd.melt(sensi_multi,id_vars=['Parameter'],value_vars=name_columns,var_name="Type")

     
    
    sensi_multi_melt[["Impact_Category","Type indice"]]=sensi_multi_melt["Type"].str.split("__",1,expand=True)



    # Statistical description of LCI values



    desc_stat_results = results_table_df.describe()

    toexclude = ['night_monitoring',
                 'Bio_class',
                 'market_for_substitution',
                 'Nsource',
                 ' min_centrifugation_rate_m3_h',
                 'max_centrifugation_rate_m3_h']
    
    desc_stat_results = desc_stat_results[[a for a in names_for_df if a not in toexclude ]]
    
   

    return (sample,
            problem_sobol_FAST,
            results_table_df,
            results_sobol_fast,
            sensi_multi_melt,
            desc_stat_results,
            total_desc_stat_contri_df)



