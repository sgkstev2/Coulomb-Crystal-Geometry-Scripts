# -*- coding: utf-8 -*-
"""
Created on Wed May 31 13:53:28 2023

@author: vr4
"""

import pandas as pd
import numpy as np
import pathlib

df = pd.read_csv(r'D:\Kane folder\To_Run_configs_Ca/Configs for Simulation.csv') # Here you need to put the location of your input excel file

configs = pd.DataFrame(df, columns=['config#'])
Xes = pd.DataFrame(df, columns=['Kr#'])
Ars = pd.DataFrame(df, columns=['Ar#'])
NH3s = pd.DataFrame(df, columns=['NH3#'])
H2Os = pd.DataFrame(df, columns=['H2O#'])
Cas = pd.DataFrame(df, columns=['Ca#'])
ETAs = pd.DataFrame(df, columns=['ETA'])
Vrfs = pd.DataFrame(df, columns=['Vrf'])
Vends = pd.DataFrame(df, columns=['Vend'])

configarr = configs.to_numpy()
Xearr = Xes.to_numpy()
Ararr = Ars.to_numpy()
NH3arr = NH3s.to_numpy()
H2Oarr = H2Os.to_numpy()
Caarr = Cas.to_numpy()
ETAarr = ETAs.to_numpy()
Vrfarr = Vrfs.to_numpy()
Vendarr = Vends.to_numpy()

rangeX = configarr.size

for i in range(0, rangeX):
    bob = np.array2string(configarr[i])
    bob = bob.replace("['", "")
    bob = bob.replace("']", "")
    #os.mkdir('{}'.format(bob))
    
    bobXe = np.array2string(Xearr[i])
    bobXe = bobXe.replace("[", "")
    bobXe = bobXe.replace("]", "")
    
    bobAr = np.array2string(Ararr[i])
    bobAr = bobAr.replace("[", "")
    bobAr = bobAr.replace("]", "") 
    
    bobNH3 = np.array2string(NH3arr[i])
    bobNH3 = bobNH3.replace("[", "")
    bobNH3 = bobNH3.replace("]", "") 
    
    bobH2O = np.array2string(H2Oarr[i])
    bobH2O = bobH2O.replace("[", "")
    bobH2O = bobH2O.replace("]", "") 
    
    bobCa = np.array2string(Caarr[i])
    bobCa = bobCa.replace("[", "")
    bobCa = bobCa.replace("]", "") 
    
    bobETA = np.array2string(ETAarr[i])
    bobETA = bobETA.replace("[", "")
    bobETA = bobETA.replace("]", "") 
    
    bobVrf = np.array2string(Vrfarr[i])
    bobVrf = bobVrf.replace("[", "")
    bobVrf = bobVrf.replace("]", "") 
    
    bobVend = np.array2string(Vendarr[i])
    bobVend = bobVend.replace("[", "")
    bobVend = bobVend.replace("]", "") 
    
    A = 'trap { \n    vrf     '
    
    A2 = '{} \n'.format(bobVrf)
    
    B = '    vend    {} \n'.format(bobVend)

    C = '    eta     {} \n'.format(bobETA)

    D = '    r0      3.5e-3 \n    z0      2.75e-3 \n    freq    3.85e6 \n    type    { \n        name cosine \n    } \n} \n\n'
    
    trap = A + A2 + B + C + D

    integrator = 'integrator {\n    stepsPerPeriod 100 \n    respasteps  500       ; Respa inner loop steps \n    coolperiods 10000 \n    histperiods   10000 \n} \n \n'

    image = 'image { \n    makeimage   true \n    scale       1; 2.5     ; Image sclaing in pixels per micron \n    blur        5.0     ; Blur radius in microns \n    dof         50.0    ; Depth of field in microns \n    nz          640     ; Number of pixels in z axis \n    nx          480     ; Number of pixels in x acis \n} \n \n'

    simulation = 'simulation {\n    threads     2 \n    seed        213 \n} \n \n'
    
    E = 'ionnumbers { \n    Ca      '
    
    F = '{} \n'.format(bobCa)
    
    F2 = '    Kr      {}\n'.format(bobXe)
    
    F3 = '    Ar      {}\n'.format(bobAr)
    
    F4 = '    NH3     {}\n'.format(bobNH3)

    F5 = '    H2O     {}\n'.format(bobH2O)
    
    Fout = ''
    
    G = '}\n\n'

    if bobXe == '0':
        Fout = F
    else:
        Fout = F + F2
    
    if bobAr == '0':
        Fout = Fout
    else:
        Fout = Fout + F3
    
    if bobNH3 == '0':
        Fout = Fout
    else:
        Fout = Fout + F4
    
    if bobH2O == '0':
        Fout = Fout
    else:
        Fout = Fout + F5

    ionnumbers = E + Fout + G

    iontype = 'iontype { \n    CaF { \n        name        CalciumFlouride \n        mass        59.0 \n        charge      1 \n    } \n    Ca { \n        name        Calcium \n        mass        40.0 \n        charge      1 \n        lasercooled true \n        beta        0.8 \n        heated      true \n        recoil      0.00001 \n        direction   1 \n                \t\tA21\t\t\t1.4e8 \n    } \n    ND3 { \n        name        Ammonia-d3 \n        mass            20.0 \n        charge          1 \n    } \n    NH3 { \n        name        Ammonia-h3 \n        mass        17.0 \n        charge      1 \n    } \n    Kr { \n        name        Krypton \n        mass        84.0 \n        charge          1 \n    } \n    Ar { \n        name        Argon \n        mass        40.0 \n        charge          1 \n    } \n    H2O { \n        name        Water-h2 \n        mass        18.0 \n        charge      1 \n    } \n}\nlaser { \n        wavelength 0.000000396908 \n        delta 7.5e8 \n        IdIsat 1 \n}'

    output = trap + integrator + image + simulation + ionnumbers + iontype

    new_dir_name = bob    
    new_dir = pathlib.Path('D:\Kane folder\To_Run_configs_Ca', new_dir_name) # This is where your simulation files will be generated 
    new_dir.mkdir(parents=True, exist_ok=True)
    new_file = new_dir / 'trap.info'
    new_file.write_text(output)    
    
