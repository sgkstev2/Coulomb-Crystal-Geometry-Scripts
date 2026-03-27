# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 15:31:54 2025

@author: Kane
"""

import pandas as pd
import numpy as np
import itertools

# Inputs:
#for each array input [start , stop , step] or input a constant value or None

def Configs_for_Simulation(Kr_arr=None, Ar_arr=None, NH3_arr=None, H2O_arr=None, Ca_arr=None, eta_arr=None, Vrf_arr=None , Vend_arr=None):
    
    df = pd.DataFrame(columns=['config#', 'Kr#', 'Ar#', 'NH3#', 'H2O#', 'Ca#', 'ETA', 'Vrf' , 'Vend'])
    
    
    def param_ranges(Kr_arr, Ar_arr, NH3_arr, H2O_arr, Ca_arr, eta_arr, Vrf_arr , Vend_arr):
        def make_range(arr):
            if isinstance(arr, (int, float)):  
                return np.array([arr])
            if arr is not None:
                return np.arange(arr[0], arr[1] + arr[2], arr[2])  # Generate range if it's a [start, stop, step] list
            
            return np.array([0])  
    
        Kr_no  = make_range(Kr_arr)
        Ar_no  = make_range(Ar_arr)
        NH3_no = make_range(NH3_arr)
        H2O_no = make_range(H2O_arr)
        Ca_no  = make_range(Ca_arr)
        eta_no = np.round(make_range(eta_arr), 3)
        Vrf_no = np.round(make_range(Vrf_arr), 3)
        Vend_no = np.round(make_range(Vend_arr), 3)
        return Kr_no, Ar_no, NH3_no, H2O_no, Ca_no, eta_no, Vrf_no , Vend_no  # Return values for unpacking

   
    Kr_no, Ar_no, NH3_no, H2O_no, Ca_no, eta_no, Vrf_no , Vend_no = param_ranges(Kr_arr, Ar_arr, NH3_arr, H2O_arr, Ca_arr, eta_arr, Vrf_arr , Vend_arr)

    
    combinations = itertools.product(Kr_no, Ar_no, NH3_no, H2O_no, Ca_no, eta_no, Vrf_no , Vend_no)    

    
    for idx, combo in enumerate(combinations, start=1):
        df.loc[len(df)] = [f'config{idx}', *combo]
            
    return df

# Define parameter ranges for the simulation we can also input constants here.
Kr_range = None 
Ar_range = None
NH3_range = None
H2O_range = None
Ca_range = [50, 1000, 50]
eta_range = 0.22
Vrf_range = [160 , 180 , 20]
Vend_range = [1.8 , 5 , 0.2]


df = Configs_for_Simulation(Kr_arr=Kr_range, Ar_arr=Ar_range, NH3_arr=NH3_range, H2O_arr=H2O_range, Ca_arr=Ca_range, eta_arr=eta_range, Vrf_arr=Vrf_range , Vend_arr=Vend_range)
df.to_csv("D:\Kane folder\\To_Run_configs_Ca/Configs for Simulation.csv", index=False)
print(len(df))

#------------------------------------------------------------------------------------------------------------------------------------------------
def Configs_for_Simulation_mixed_ratio(total_particles_arr=None, mix_species="Ar" ,  ratios=None , eta_arr=None , Vrf_arr=None, Vend_arr=None):

    def make_range(arr):
        if isinstance(arr, (int, float)):
            return np.array([arr])
        if arr is not None:
            return np.arange(arr[0], arr[1] + arr[2], arr[2])
        return np.array([0])
        
    eta_no = np.round(make_range(eta_arr), 3)
    Vrf_no = np.round(make_range(Vrf_arr), 3)
    Vend_no = np.round(make_range(Vend_arr), 3)
    
    
    def mix_numbers():
        configs = []
        idx = 1
        
        for total in total_particles_arr:
            for r in ratios:
                Ca_val = round(total * r)
                other_val = round(total - Ca_val)

                Kr_val, Ar_val, NH3_val, H2O_val = 0, 0, 0, 0
                if mix_species == "Ar": Ar_val = other_val
                elif mix_species == "Kr": Kr_val = other_val
                elif mix_species == "NH3": NH3_val = other_val
                elif mix_species == "H2O": H2O_val = other_val

                for eta, Vrf, Vend in itertools.product(eta_no, Vrf_no, Vend_no):
                    configs.append([f'config{idx}', Kr_val, Ar_val, NH3_val, H2O_val,
                                    Ca_val, eta, Vrf, Vend])
                    idx += 1
        return configs
                
    configs = mix_numbers()
                    
    df = pd.DataFrame(configs, columns=['config#', 'Kr#', 'Ar#', 'NH3#', 'H2O#', 'Ca#', 'ETA', 'Vrf', 'Vend'])
    return df

#total_particles_arr = np.arange(1300 , 1500 + 100, 100)
#ratios = [0.8]#np.arange(0.55 , 0.8 + 0.05 , 0.05)
#eta_arr = [0.2 , 0.34 , 0.02]
#Vrf_arr = [160 , 240 , 20]
#Vend_arr = [1.5, 4.5 , 0.5]

#df = Configs_for_Simulation_mixed_ratio(total_particles_arr, mix_species="Kr", ratios=ratios, eta_arr=eta_arr, Vrf_arr=Vrf_arr, Vend_arr=Vend_arr)
#df.to_csv("D:\To_Run_configs_Ca_Kr/Configs for Simulation.csv", index=False)
#print(len(df))